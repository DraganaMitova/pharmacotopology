from __future__ import annotations

import gzip
import hashlib
from collections import OrderedDict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
)


@dataclass(frozen=True)
class PfamDomainMapping:
    pfam_id: str
    name: str
    description: str
    start: int
    end: int
    coverage: float

    def positions(self) -> tuple[int, ...]:
        return tuple(range(self.start, self.end + 1))


@dataclass(frozen=True)
class ExternalCovariationPair:
    i: int
    j: int
    raw_score: float
    apc_corrected_score: float
    pfam_id: str


@dataclass(frozen=True)
class PfamApcRowResult:
    row_id: str
    source_accession: str
    dca_status: str
    external_coupling_status: str
    failure_kind: str
    rejection_reason: str
    dca_tool: str
    covariation_method: str
    pfam_ids: str
    sample_depth: int
    total_depth_seen: int
    raw_pair_count: int
    accepted_pair_count: int
    msa_sha256: str
    pairs: tuple[ExternalCovariationPair, ...]

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload.pop("pairs", None)
        return payload


AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
GAP_CODE = 20


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(block)
    return f"sha256:{digest.hexdigest()}"


def read_stockholm_sample(
    path: Path,
    *,
    max_records: int = 4000,
) -> tuple[tuple[str, ...], int]:
    opener = gzip.open if path.read_bytes()[:2] == b"\x1f\x8b" else open
    records: OrderedDict[str, list[str]] = OrderedDict()
    seen_total: set[str] = set()
    with opener(path, "rt", encoding="utf-8", errors="replace") as file:
        for line in file:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped == "//":
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            name, fragment = parts[0], parts[1]
            seen_total.add(name)
            if name not in records:
                if len(records) >= max_records:
                    continue
                records[name] = []
            if name in records:
                records[name].append(fragment)
    return tuple("".join(parts) for parts in records.values()), len(seen_total)


def _encode_alignment(sequences: Sequence[str]):
    import numpy as np

    if not sequences:
        return np.zeros((0, 0), dtype=np.int16)
    encoded = {residue: index for index, residue in enumerate(AA_ALPHABET)}
    width = max(len(sequence) for sequence in sequences)
    matrix = np.full((len(sequences), width), GAP_CODE, dtype=np.int16)
    for row_index, sequence in enumerate(sequences):
        cleaned = "".join(char for char in sequence if not char.islower())
        for column_index, residue in enumerate(cleaned[:width]):
            matrix[row_index, column_index] = encoded.get(residue.upper(), GAP_CODE)
    return matrix


def _apc_pairs_from_alignment(
    *,
    sequences: Sequence[str],
    positions: Sequence[int],
    minimum_sequence_separation: int,
    column_occupancy_min: float = 0.20,
) -> tuple[ExternalCovariationPair, ...]:
    import numpy as np

    matrix = _encode_alignment(sequences)
    if matrix.size == 0 or len(positions) < 2:
        return ()
    sequence_count, _ = matrix.shape
    occupancy = (matrix != GAP_CODE).mean(axis=0)
    kept_columns = np.where(occupancy >= column_occupancy_min)[0]
    usable_count = min(len(kept_columns), len(positions))
    if usable_count < 2:
        return ()
    kept_columns = kept_columns[:usable_count]
    mapped_positions = tuple(positions[:usable_count])
    data = matrix[:, kept_columns]
    alphabet_size = len(AA_ALPHABET) + 1
    width = data.shape[1]
    mi = np.zeros((width, width), dtype=np.float64)
    for left in range(width):
        left_column = data[:, left]
        left_probs = (
            np.bincount(left_column, minlength=alphabet_size).astype(np.float64)
            / sequence_count
        )
        for right in range(left + 1, width):
            if mapped_positions[right] - mapped_positions[left] < minimum_sequence_separation:
                continue
            right_column = data[:, right]
            right_probs = (
                np.bincount(right_column, minlength=alphabet_size).astype(np.float64)
                / sequence_count
            )
            joint = (
                np.bincount(
                    left_column * alphabet_size + right_column,
                    minlength=alphabet_size * alphabet_size,
                )
                .astype(np.float64)
                .reshape((alphabet_size, alphabet_size))
                / sequence_count
            )
            nonzero = joint > 0
            denom = left_probs[:, None] * right_probs[None, :]
            value = float((joint[nonzero] * np.log(joint[nonzero] / denom[nonzero])).sum())
            mi[left, right] = value
            mi[right, left] = value
    if not (mi > 0).any():
        return ()
    row_mean = mi.mean(axis=1)
    total_mean = float(mi.mean()) or 1.0
    pairs: list[ExternalCovariationPair] = []
    for left in range(width):
        for right in range(left + 1, width):
            if mapped_positions[right] - mapped_positions[left] < minimum_sequence_separation:
                continue
            apc = mi[left, right] - (row_mean[left] * row_mean[right] / total_mean)
            pairs.append(
                ExternalCovariationPair(
                    i=mapped_positions[left],
                    j=mapped_positions[right],
                    raw_score=round(float(mi[left, right]), 6),
                    apc_corrected_score=round(float(apc), 6),
                    pfam_id="",
                )
            )
    return tuple(sorted(pairs, key=lambda pair: (-pair.apc_corrected_score, pair.i, pair.j)))


def run_pfam_apc_covariation_for_row(
    *,
    row: RealCoordinateVisualRow,
    mappings: Sequence[PfamDomainMapping],
    alignment_dir: Path,
    max_records: int = 4000,
    minimum_sequence_separation: int = 24,
) -> PfamApcRowResult:
    if not mappings:
        return PfamApcRowResult(
            row_id=row.row_id,
            source_accession=row.source_accession,
            dca_status="mapping_failed",
            external_coupling_status="mapping_failed_no_sifts_uniprot_match",
            failure_kind="mapping_failed",
            rejection_reason="no_pdbe_pfam_mapping",
            dca_tool="none",
            covariation_method="interpro_pfam_full_alignment_mi_apc",
            pfam_ids="",
            sample_depth=0,
            total_depth_seen=0,
            raw_pair_count=0,
            accepted_pair_count=0,
            msa_sha256="",
            pairs=(),
        )

    best_by_pair: dict[tuple[int, int], ExternalCovariationPair] = {}
    sample_depth = 0
    total_depth_seen = 0
    hashes: list[str] = []
    for mapping in mappings:
        alignment_path = alignment_dir / f"{mapping.pfam_id}.full.sto.gz"
        if not alignment_path.exists():
            continue
        sequences, total_depth = read_stockholm_sample(
            alignment_path,
            max_records=max_records,
        )
        sample_depth = max(sample_depth, len(sequences))
        total_depth_seen = max(total_depth_seen, total_depth)
        hashes.append(sha256_file(alignment_path))
        for pair in _apc_pairs_from_alignment(
            sequences=sequences,
            positions=mapping.positions(),
            minimum_sequence_separation=minimum_sequence_separation,
        ):
            if pair.i < 1 or pair.j > row.sequence_length or pair.j <= pair.i:
                continue
            keyed = ExternalCovariationPair(
                i=pair.i,
                j=pair.j,
                raw_score=pair.raw_score,
                apc_corrected_score=pair.apc_corrected_score,
                pfam_id=mapping.pfam_id,
            )
            key = (keyed.i, keyed.j)
            if key not in best_by_pair or (
                keyed.apc_corrected_score > best_by_pair[key].apc_corrected_score
            ):
                best_by_pair[key] = keyed
    ranked = tuple(
        sorted(
            best_by_pair.values(),
            key=lambda pair: (-pair.apc_corrected_score, pair.i, pair.j),
        )[: row.sequence_length]
    )
    status = "dca_available" if ranked else "dca_failed_tool_error"
    return PfamApcRowResult(
        row_id=row.row_id,
        source_accession=row.source_accession,
        dca_status=status,
        external_coupling_status=(
            "external_couplings_available" if ranked else "dca_failed_tool_error"
        ),
        failure_kind="" if ranked else "dca_failed",
        rejection_reason="" if ranked else "no_ranked_covariation_pairs",
        dca_tool="none",
        covariation_method="interpro_pfam_full_alignment_mi_apc",
        pfam_ids=",".join(mapping.pfam_id for mapping in mappings),
        sample_depth=sample_depth,
        total_depth_seen=total_depth_seen,
        raw_pair_count=len(best_by_pair),
        accepted_pair_count=len(ranked),
        msa_sha256=";".join(hashes),
        pairs=ranked,
    )
