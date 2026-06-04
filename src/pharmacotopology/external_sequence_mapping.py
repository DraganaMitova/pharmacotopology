from __future__ import annotations

import csv
import gzip
import hashlib
import re
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterator, Mapping, Optional, Sequence

from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
)


SIFTS_PDB_CHAIN_UNIPROT_CSV_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/"
    "pdb_chain_uniprot.csv.gz"
)


@dataclass(frozen=True)
class ExternalSequenceMappingResult:
    row_id: str
    source_accession: str
    pdb_id: str
    chain_id: str
    sequence_length: int
    sequence_sha256: str
    mapping_status: str
    external_coupling_status: str
    failure_kind: str
    rejection_reason: str
    uniprot_accession: str
    sifts_mapping_available: bool
    mapping_source_kind: str
    mapping_sha256: str
    target_coverage: float
    focus_sequence_mapping_confidence: float
    coordinate_truth_used_to_build_constraints: bool = False
    native_truth_used_before_coupling_selection: bool = False
    structure_model_used: bool = False
    raw_sequence_exposed: bool = False

    @property
    def resolved(self) -> bool:
        return self.mapping_status == "mapping_resolved"

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _score(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _pdb_chain(source_accession: str) -> tuple[str, str]:
    pdb_id, chain_id = source_accession.split(":", 1)
    return pdb_id.upper(), chain_id


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(block)
    return f"sha256:{digest.hexdigest()}"


def download_sifts_pdb_chain_uniprot_csv(
    *,
    output_path: Path,
    url: str = SIFTS_PDB_CHAIN_UNIPROT_CSV_URL,
) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=120) as response:
        output_path.write_bytes(response.read())
    return output_path


def _open_text(path: Path) -> Iterator[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", newline="") as file:
            for line in file:
                yield line
        return
    with path.open("r", encoding="utf-8", newline="") as file:
        for line in file:
            yield line


def _non_comment_lines(path: Path) -> Iterator[str]:
    for line in _open_text(path):
        if line.startswith("#") or not line.strip():
            continue
        yield line


def _normalise_key(key: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", key.lower())


def _value(raw: Mapping[str, str], *names: str) -> str:
    by_key = {_normalise_key(key): value for key, value in raw.items()}
    for name in names:
        value = by_key.get(_normalise_key(name), "")
        if value:
            return value.strip()
    return ""


def _int_or_none(value: str) -> Optional[int]:
    match = re.search(r"-?\d+", value or "")
    if match is None:
        return None
    return int(match.group(0))


def _span_coverage(raw: Mapping[str, str], sequence_length: int) -> float:
    starts = (
        _value(raw, "res_beg", "seqres_beg", "pdb_beg", "sp_beg"),
        _value(raw, "pdb_beg", "res_beg", "seqres_beg", "sp_beg"),
    )
    ends = (
        _value(raw, "res_end", "seqres_end", "pdb_end", "sp_end"),
        _value(raw, "pdb_end", "res_end", "seqres_end", "sp_end"),
    )
    spans: list[int] = []
    for start_raw, end_raw in zip(starts, ends):
        start = _int_or_none(start_raw)
        end = _int_or_none(end_raw)
        if start is None or end is None:
            continue
        spans.append(abs(end - start) + 1)
    if not spans or sequence_length <= 0:
        return 1.0
    return _score(max(spans) / sequence_length)


def _sifts_rows_by_pdb_chain(path: Path) -> dict[tuple[str, str], list[dict[str, str]]]:
    rows: dict[tuple[str, str], list[dict[str, str]]] = {}
    reader = csv.DictReader(_non_comment_lines(path))
    for raw in reader:
        pdb_id = _value(raw, "pdb", "pdb_id", "entry_id").upper()
        chain_id = _value(raw, "chain", "chain_id", "auth_asym_id", "asym_id")
        if not pdb_id or not chain_id:
            continue
        rows.setdefault((pdb_id, chain_id), []).append(dict(raw))
    return rows


def _missing_mapping(row: RealCoordinateVisualRow, reason: str) -> ExternalSequenceMappingResult:
    pdb_id, chain_id = _pdb_chain(row.source_accession)
    failure_kind = "database_missing" if reason == "sifts_database_missing" else "mapping_failed"
    return ExternalSequenceMappingResult(
        row_id=row.row_id,
        source_accession=row.source_accession,
        pdb_id=pdb_id,
        chain_id=chain_id,
        sequence_length=row.sequence_length,
        sequence_sha256=f"sha256:{row.sequence_sha256}",
        mapping_status="database_missing" if failure_kind == "database_missing" else "mapping_failed",
        external_coupling_status="mapping_failed_no_sifts_uniprot_match",
        failure_kind=failure_kind,
        rejection_reason=reason,
        uniprot_accession="",
        sifts_mapping_available=False,
        mapping_source_kind="pdbe_sifts_pdb_chain_uniprot_csv",
        mapping_sha256="",
        target_coverage=0.0,
        focus_sequence_mapping_confidence=0.0,
    )


def map_rows_to_uniprot_via_sifts(
    rows: Sequence[RealCoordinateVisualRow],
    *,
    sifts_csv_path: Optional[Path],
) -> tuple[ExternalSequenceMappingResult, ...]:
    if sifts_csv_path is None or not sifts_csv_path.exists():
        return tuple(_missing_mapping(row, "sifts_database_missing") for row in rows)

    mapping_sha256 = _sha256_file(sifts_csv_path)
    sifts = _sifts_rows_by_pdb_chain(sifts_csv_path)
    results: list[ExternalSequenceMappingResult] = []
    for row in rows:
        pdb_id, chain_id = _pdb_chain(row.source_accession)
        candidates = sifts.get((pdb_id, chain_id), ())
        best: Optional[Mapping[str, str]] = None
        best_coverage = 0.0
        for candidate in candidates:
            accession = _value(
                candidate,
                "sp_primary",
                "uniprot",
                "uniprot_accession",
                "unp",
                "accession",
            )
            if not accession:
                continue
            coverage = _span_coverage(candidate, row.sequence_length)
            if best is None or coverage > best_coverage:
                best = candidate
                best_coverage = coverage
        if best is None:
            results.append(_missing_mapping(row, "no_sifts_uniprot_match"))
            continue
        uniprot_accession = _value(
            best,
            "sp_primary",
            "uniprot",
            "uniprot_accession",
            "unp",
            "accession",
        )
        coverage = _score(best_coverage)
        results.append(
            ExternalSequenceMappingResult(
                row_id=row.row_id,
                source_accession=row.source_accession,
                pdb_id=pdb_id,
                chain_id=chain_id,
                sequence_length=row.sequence_length,
                sequence_sha256=f"sha256:{row.sequence_sha256}",
                mapping_status="mapping_resolved",
                external_coupling_status="mapping_resolved",
                failure_kind="",
                rejection_reason="",
                uniprot_accession=uniprot_accession,
                sifts_mapping_available=True,
                mapping_source_kind="pdbe_sifts_pdb_chain_uniprot_csv",
                mapping_sha256=mapping_sha256,
                target_coverage=coverage,
                focus_sequence_mapping_confidence=coverage,
            )
        )
    return tuple(results)
