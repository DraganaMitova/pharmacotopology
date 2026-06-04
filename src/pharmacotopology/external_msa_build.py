from __future__ import annotations

import hashlib
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

from pharmacotopology.external_sequence_mapping import ExternalSequenceMappingResult
from pharmacotopology.folding_external_coupling_sources import (
    SERIOUS_EXTERNAL_COUPLING_POLICY,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
)


@dataclass(frozen=True)
class MsaBuildConfig:
    msa_cache_dir: Optional[Path] = None
    msa_database: Optional[Path] = None
    msa_tool: str = "jackhmmer"
    work_dir: Optional[Path] = None
    target_coverage_min: float = SERIOUS_EXTERNAL_COUPLING_POLICY.target_coverage_min
    effective_sequence_count_over_length_min: float = (
        SERIOUS_EXTERNAL_COUPLING_POLICY.effective_sequence_count_over_length_min
    )


@dataclass(frozen=True)
class ExternalMsaBuildResult:
    row_id: str
    source_accession: str
    uniprot_accession: str
    msa_status: str
    external_coupling_status: str
    failure_kind: str
    rejection_reason: str
    msa_source_kind: str
    msa_tool: str
    msa_database: str
    msa_sha256: str
    msa_depth: int
    effective_sequence_count: float
    effective_sequence_count_over_length: float
    target_coverage: float
    focus_sequence_mapping_confidence: float
    msa_path: Optional[Path] = None
    coordinate_truth_used_to_build_constraints: bool = False
    native_truth_used_before_coupling_selection: bool = False
    structure_model_used: bool = False
    raw_sequence_exposed: bool = False

    @property
    def usable(self) -> bool:
        return self.msa_status == "msa_attempted"

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload.pop("msa_path", None)
        return payload


def _score(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _six(value: float) -> float:
    return round(value, 6)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(block)
    return f"sha256:{digest.hexdigest()}"


def _safe_name(value: str) -> str:
    return "".join(char if char.isalnum() else "_" for char in value)


def _candidate_msa_paths(
    *,
    cache_dir: Path,
    row: RealCoordinateVisualRow,
    mapping: ExternalSequenceMappingResult,
) -> Iterable[Path]:
    stems = (
        row.row_id,
        _safe_name(row.source_accession),
        mapping.uniprot_accession,
    )
    suffixes = (".a3m", ".afa", ".aln", ".fa", ".fasta", ".sto", ".stockholm")
    for stem in stems:
        if not stem:
            continue
        for suffix in suffixes:
            yield cache_dir / f"{stem}{suffix}"


def _read_alignment_records(path: Path) -> tuple[str, ...]:
    text = path.read_text(encoding="utf-8")
    if path.suffix in {".sto", ".stockholm"}:
        records: dict[str, list[str]] = {}
        for line in text.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped == "//":
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            records.setdefault(parts[0], []).append(parts[1])
        return tuple("".join(parts) for parts in records.values())

    records_out: list[str] = []
    current: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if current:
                records_out.append("".join(current))
            current = []
            continue
        current.append(stripped)
    if current:
        records_out.append("".join(current))
    if records_out:
        return tuple(records_out)
    return tuple(line.strip() for line in text.splitlines() if line.strip())


def _alignment_coverage(records: Sequence[str], sequence_length: int) -> float:
    if not records or sequence_length <= 0:
        return 0.0
    focus = records[0]
    aligned_focus = "".join(char for char in focus if not char.islower())
    covered = sum(1 for char in aligned_focus if char not in {"-", "."})
    return _score(covered / sequence_length)


def _assess_alignment(
    *,
    row: RealCoordinateVisualRow,
    mapping: ExternalSequenceMappingResult,
    path: Path,
    config: MsaBuildConfig,
    msa_source_kind: str,
    msa_tool: str,
    msa_database: str,
) -> ExternalMsaBuildResult:
    records = _read_alignment_records(path)
    depth = len(records)
    effective_sequence_count = float(depth)
    depth_over_length = (
        effective_sequence_count / row.sequence_length
        if row.sequence_length
        else 0.0
    )
    target_coverage = _alignment_coverage(records, row.sequence_length)
    if target_coverage < config.target_coverage_min:
        status = "msa_failed_low_coverage"
        failure_kind = "msa_failed"
        reason = "target_coverage_below_threshold"
    elif depth_over_length < config.effective_sequence_count_over_length_min:
        status = "msa_failed_low_depth"
        failure_kind = "msa_failed"
        reason = "effective_sequence_count_over_length_below_threshold"
    else:
        status = "msa_attempted"
        failure_kind = ""
        reason = ""
    return ExternalMsaBuildResult(
        row_id=row.row_id,
        source_accession=row.source_accession,
        uniprot_accession=mapping.uniprot_accession,
        msa_status=status,
        external_coupling_status=status,
        failure_kind=failure_kind,
        rejection_reason=reason,
        msa_source_kind=msa_source_kind,
        msa_tool=msa_tool,
        msa_database=msa_database,
        msa_sha256=_sha256_file(path),
        msa_depth=depth,
        effective_sequence_count=_six(effective_sequence_count),
        effective_sequence_count_over_length=_six(depth_over_length),
        target_coverage=target_coverage,
        focus_sequence_mapping_confidence=mapping.focus_sequence_mapping_confidence,
        msa_path=path,
    )


def _failed_due_mapping(
    row: RealCoordinateVisualRow,
    mapping: ExternalSequenceMappingResult,
) -> ExternalMsaBuildResult:
    return ExternalMsaBuildResult(
        row_id=row.row_id,
        source_accession=row.source_accession,
        uniprot_accession=mapping.uniprot_accession,
        msa_status="mapping_failed",
        external_coupling_status=mapping.external_coupling_status,
        failure_kind=mapping.failure_kind or "mapping_failed",
        rejection_reason=mapping.rejection_reason,
        msa_source_kind="not_attempted_mapping_failed",
        msa_tool="",
        msa_database="",
        msa_sha256="",
        msa_depth=0,
        effective_sequence_count=0.0,
        effective_sequence_count_over_length=0.0,
        target_coverage=0.0,
        focus_sequence_mapping_confidence=0.0,
    )


def _missing_database(
    row: RealCoordinateVisualRow,
    mapping: ExternalSequenceMappingResult,
    config: MsaBuildConfig,
) -> ExternalMsaBuildResult:
    return ExternalMsaBuildResult(
        row_id=row.row_id,
        source_accession=row.source_accession,
        uniprot_accession=mapping.uniprot_accession,
        msa_status="database_missing",
        external_coupling_status="msa_failed",
        failure_kind="database_missing",
        rejection_reason="msa_database_missing",
        msa_source_kind="external_msa_search_not_available",
        msa_tool=config.msa_tool,
        msa_database="",
        msa_sha256="",
        msa_depth=0,
        effective_sequence_count=0.0,
        effective_sequence_count_over_length=0.0,
        target_coverage=0.0,
        focus_sequence_mapping_confidence=mapping.focus_sequence_mapping_confidence,
    )


def _missing_tool(
    row: RealCoordinateVisualRow,
    mapping: ExternalSequenceMappingResult,
    config: MsaBuildConfig,
) -> ExternalMsaBuildResult:
    return ExternalMsaBuildResult(
        row_id=row.row_id,
        source_accession=row.source_accession,
        uniprot_accession=mapping.uniprot_accession,
        msa_status="tool_missing",
        external_coupling_status="msa_failed",
        failure_kind="tool_missing",
        rejection_reason=f"msa_tool_missing:{config.msa_tool}",
        msa_source_kind="external_msa_search_not_available",
        msa_tool=config.msa_tool,
        msa_database=str(config.msa_database or ""),
        msa_sha256="",
        msa_depth=0,
        effective_sequence_count=0.0,
        effective_sequence_count_over_length=0.0,
        target_coverage=0.0,
        focus_sequence_mapping_confidence=mapping.focus_sequence_mapping_confidence,
    )


def _run_jackhmmer(
    *,
    row: RealCoordinateVisualRow,
    mapping: ExternalSequenceMappingResult,
    tool_path: str,
    config: MsaBuildConfig,
) -> Optional[Path]:
    work_dir = config.work_dir or Path(tempfile.mkdtemp(prefix="external_msa_v0_"))
    work_dir.mkdir(parents=True, exist_ok=True)
    query_path = work_dir / f"{row.row_id}.query.fasta"
    output_path = work_dir / f"{row.row_id}.sto"
    query_path.write_text(
        f">{mapping.uniprot_accession or row.row_id}\n{row.sequence}\n",
        encoding="utf-8",
    )
    command = [
        tool_path,
        "-N",
        "1",
        "-A",
        str(output_path),
        str(query_path),
        str(config.msa_database),
    ]
    completed = subprocess.run(
        command,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if completed.returncode != 0 or not output_path.exists():
        return None
    return output_path


def build_external_msas(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    mappings: Sequence[ExternalSequenceMappingResult],
    config: MsaBuildConfig = MsaBuildConfig(),
) -> tuple[ExternalMsaBuildResult, ...]:
    mapping_by_row = {mapping.row_id: mapping for mapping in mappings}
    results: list[ExternalMsaBuildResult] = []
    for row in rows:
        mapping = mapping_by_row[row.row_id]
        if not mapping.resolved:
            results.append(_failed_due_mapping(row, mapping))
            continue

        if config.msa_cache_dir is not None:
            for candidate in _candidate_msa_paths(
                cache_dir=config.msa_cache_dir,
                row=row,
                mapping=mapping,
            ):
                if candidate.exists():
                    results.append(
                        _assess_alignment(
                            row=row,
                            mapping=mapping,
                            path=candidate,
                            config=config,
                            msa_source_kind="cached_external_msa",
                            msa_tool="cache",
                            msa_database=str(config.msa_cache_dir),
                        )
                    )
                    break
            if results and results[-1].row_id == row.row_id:
                continue

        if config.msa_database is None or not config.msa_database.exists():
            results.append(_missing_database(row, mapping, config))
            continue

        tool_path = shutil.which(config.msa_tool)
        if tool_path is None:
            results.append(_missing_tool(row, mapping, config))
            continue

        output = _run_jackhmmer(
            row=row,
            mapping=mapping,
            tool_path=tool_path,
            config=config,
        )
        if output is None:
            results.append(
                ExternalMsaBuildResult(
                    row_id=row.row_id,
                    source_accession=row.source_accession,
                    uniprot_accession=mapping.uniprot_accession,
                    msa_status="msa_failed",
                    external_coupling_status="msa_failed",
                    failure_kind="msa_failed",
                    rejection_reason="msa_tool_failed",
                    msa_source_kind="external_msa_search_failed",
                    msa_tool=config.msa_tool,
                    msa_database=str(config.msa_database),
                    msa_sha256="",
                    msa_depth=0,
                    effective_sequence_count=0.0,
                    effective_sequence_count_over_length=0.0,
                    target_coverage=0.0,
                    focus_sequence_mapping_confidence=(
                        mapping.focus_sequence_mapping_confidence
                    ),
                )
            )
            continue
        results.append(
            _assess_alignment(
                row=row,
                mapping=mapping,
                path=output,
                config=config,
                msa_source_kind="jackhmmer_sequence_search",
                msa_tool=config.msa_tool,
                msa_database=str(config.msa_database),
            )
        )
    return tuple(results)
