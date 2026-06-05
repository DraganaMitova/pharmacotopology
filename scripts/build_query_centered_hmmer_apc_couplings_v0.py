from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from build_query_centered_hmmer_plmc_couplings_v0 import (  # noqa: E402
    DEFAULT_BASE_EXTERNAL_COUPLING_FILE,
    DEFAULT_BENCHMARK_FILE,
    DEFAULT_MINIMUM_SEQUENCE_SEPARATION,
    _alignment_column_to_row_position,
    _read_focus_fasta,
    _sha256_file,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
GAP_CODE = len(AA_ALPHABET)


@dataclass(frozen=True)
class ApcPair:
    apc_corrected_score: float
    raw_score: float
    i: int
    j: int
    alignment_column_i: int
    alignment_column_j: int


def _read_fasta_records(path: Path) -> tuple[tuple[str, str], ...]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    sequence_parts: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if name is not None:
                records.append((name, "".join(sequence_parts)))
            name = stripped[1:]
            sequence_parts = []
            continue
        sequence_parts.append(stripped)
    if name is not None:
        records.append((name, "".join(sequence_parts)))
    if not records:
        raise ValueError(f"{path} does not contain FASTA records")
    return tuple(records)


def _query_centered_apc_pairs(
    *,
    row: RealCoordinateVisualRow,
    focus_fasta_file: Path,
    minimum_sequence_separation: int,
    max_records: int,
) -> tuple[ApcPair, ...]:
    if minimum_sequence_separation < 1:
        raise ValueError("minimum_sequence_separation must be positive")
    if max_records < 2:
        raise ValueError("max_records must be at least 2")

    _, focus_alignment = _read_focus_fasta(focus_fasta_file)
    column_to_position = _alignment_column_to_row_position(
        focus_alignment=focus_alignment,
        row=row,
    )
    alignment_columns = tuple(sorted(column_to_position))
    positions = tuple(column_to_position[column] for column in alignment_columns)
    records = _read_fasta_records(focus_fasta_file)[:max_records]

    import numpy as np

    encoded = {residue: index for index, residue in enumerate(AA_ALPHABET)}
    matrix = np.full(
        (len(records), len(alignment_columns)),
        GAP_CODE,
        dtype=np.int16,
    )
    for row_index, (_, sequence) in enumerate(records):
        for column_index, alignment_column in enumerate(alignment_columns):
            if alignment_column > len(sequence):
                continue
            residue = sequence[alignment_column - 1]
            if not residue.isalpha() or residue == "-":
                continue
            matrix[row_index, column_index] = encoded.get(residue.upper(), GAP_CODE)

    sequence_count, width = matrix.shape
    if sequence_count < 2 or width < 2:
        return ()

    alphabet_size = len(AA_ALPHABET) + 1
    mi = np.zeros((width, width), dtype=np.float64)
    eligible = np.zeros((width, width), dtype=np.bool_)
    for left in range(width):
        left_column = matrix[:, left]
        left_probs = (
            np.bincount(left_column, minlength=alphabet_size).astype(np.float64)
            / sequence_count
        )
        for right in range(left + 1, width):
            if positions[right] - positions[left] < minimum_sequence_separation:
                continue
            eligible[left, right] = True
            eligible[right, left] = True
            right_column = matrix[:, right]
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
            denominator = left_probs[:, None] * right_probs[None, :]
            value = float(
                (
                    joint[nonzero]
                    * np.log(joint[nonzero] / denominator[nonzero])
                ).sum()
            )
            mi[left, right] = value
            mi[right, left] = value

    if not (mi > 0.0).any():
        return ()

    eligible_counts = eligible.sum(axis=1)
    row_mean = np.divide(
        (mi * eligible).sum(axis=1),
        eligible_counts,
        out=np.zeros(width, dtype=np.float64),
        where=eligible_counts > 0,
    )
    total_mean = float(mi[eligible].mean()) if eligible.any() else 1.0
    if total_mean == 0.0:
        total_mean = 1.0
    pairs: list[ApcPair] = []
    for left in range(width):
        for right in range(left + 1, width):
            if positions[right] - positions[left] < minimum_sequence_separation:
                continue
            apc = float(mi[left, right] - (row_mean[left] * row_mean[right] / total_mean))
            pairs.append(
                ApcPair(
                    apc_corrected_score=apc,
                    raw_score=float(mi[left, right]),
                    i=positions[left],
                    j=positions[right],
                    alignment_column_i=alignment_columns[left],
                    alignment_column_j=alignment_columns[right],
                )
            )
    return tuple(
        sorted(
            pairs,
            key=lambda pair: (
                -pair.apc_corrected_score,
                pair.i,
                pair.j,
            ),
        )
    )


def _constraint_rows(
    *,
    row: RealCoordinateVisualRow,
    pairs: Sequence[ApcPair],
    msa_sha256: str,
    msa_depth: int,
    effective_sequence_count: float,
    focus_id: str,
    minimum_sequence_separation: int,
) -> list[dict[str, object]]:
    selected = tuple(pairs[: row.sequence_length])
    max_positive_score = max(
        (max(0.0, pair.apc_corrected_score) for pair in selected),
        default=0.0,
    )
    constraints: list[dict[str, object]] = []
    for rank, pair in enumerate(selected, start=1):
        confidence = (
            max(0.01, max(0.0, pair.apc_corrected_score) / max_positive_score)
            if max_positive_score
            else 0.01
        )
        constraints.append(
            {
                "row_id": row.row_id,
                "source_accession": row.source_accession,
                "constraint_id": (
                    f"query_hmmer_apc_{row.row_id}_{pair.i}_{pair.j}_{rank}"
                ),
                "i": pair.i,
                "j": pair.j,
                "sequence_separation": pair.j - pair.i,
                "normalized_separation": round(
                    (pair.j - pair.i) / row.sequence_length,
                    6,
                ),
                "confidence": round(confidence, 6),
                "raw_score": round(pair.raw_score, 6),
                "apc_corrected_score": round(pair.apc_corrected_score, 6),
                "rank": rank,
                "rank_fraction": round(rank / row.sequence_length, 6),
                "constraint_class": "external_query_centered_hmmer_mi_apc_coupling",
                "source_kind": "external_uniref_msa_dca_v1",
                "msa_source_kind": (
                    "ebi_hmmer_jackhmmer_uniprot_focus_mi_apc"
                ),
                "msa_sha256": msa_sha256,
                "msa_depth": msa_depth,
                "effective_sequence_count": float(effective_sequence_count),
                "effective_sequence_count_over_length": round(
                    effective_sequence_count / row.sequence_length,
                    6,
                ),
                "target_coverage": 1.0,
                "focus_sequence_mapping_confidence": 1.0,
                "top_L_couplings_available": len(pairs) >= row.sequence_length,
                "coordinate_truth_used_to_build_constraint": False,
                "native_truth_used_before_coupling_selection": False,
                "structure_model_used": False,
                "raw_sequence_exposed": False,
                "apc_focus_id": focus_id,
                "apc_alignment_column_i": pair.alignment_column_i,
                "apc_alignment_column_j": pair.alignment_column_j,
                "minimum_sequence_separation": minimum_sequence_separation,
            }
        )
    return constraints


def build_query_centered_hmmer_apc_couplings_v0(
    *,
    benchmark_file: Path,
    base_external_coupling_file: Path,
    row_id: str,
    focus_fasta_file: Path,
    output: Path,
    msa_depth: int,
    effective_sequence_count: float,
    minimum_sequence_separation: int = DEFAULT_MINIMUM_SEQUENCE_SEPARATION,
    max_records: int = 2000,
) -> Path:
    rows = tuple(load_real_coordinate_visual_rows(benchmark_file))
    row_by_id = {row.row_id: row for row in rows}
    row = row_by_id[row_id]
    focus_id, _ = _read_focus_fasta(focus_fasta_file)
    pairs = _query_centered_apc_pairs(
        row=row,
        focus_fasta_file=focus_fasta_file,
        minimum_sequence_separation=minimum_sequence_separation,
        max_records=max_records,
    )
    if len(pairs) < row.sequence_length:
        raise ValueError(
            f"{row_id} has only {len(pairs)} long-range APC pairs with "
            f"sequence separation >= {minimum_sequence_separation}; "
            f"{row.sequence_length} required"
        )
    constraints = _constraint_rows(
        row=row,
        pairs=pairs,
        msa_sha256=_sha256_file(focus_fasta_file),
        msa_depth=msa_depth,
        effective_sequence_count=effective_sequence_count,
        focus_id=focus_id,
        minimum_sequence_separation=minimum_sequence_separation,
    )
    payload = json.loads(base_external_coupling_file.read_text(encoding="utf-8"))
    payload["constraints"] = [
        constraint
        for constraint in payload["constraints"]
        if constraint["row_id"] != row_id
    ] + constraints
    payload["constraints"] = sorted(
        payload["constraints"],
        key=lambda constraint: (
            str(constraint["row_id"]),
            int(constraint["rank"]),
            int(constraint["i"]),
            int(constraint["j"]),
        ),
    )
    payload["external_constraint_count"] = len(payload["constraints"])
    payload[
        f"hmmer_query_centered_apc_{row.source_accession.replace(':', '_')}_added"
    ] = True
    payload[
        f"hmmer_query_centered_apc_{row.source_accession.replace(':', '_')}_alignment_sha256"
    ] = _sha256_file(focus_fasta_file)
    payload[
        f"hmmer_query_centered_apc_{row.source_accession.replace(':', '_')}_minimum_sequence_separation"
    ] = minimum_sequence_separation
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Replace one row in the locked external-coupling JSON with "
            "query-centered EBI HMMER jackhmmer MI/APC couplings."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--base-external-coupling-file",
        default=str(DEFAULT_BASE_EXTERNAL_COUPLING_FILE),
    )
    parser.add_argument("--row-id", required=True)
    parser.add_argument("--focus-fasta-file", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--msa-depth", type=int, required=True)
    parser.add_argument("--effective-sequence-count", type=float, required=True)
    parser.add_argument(
        "--minimum-sequence-separation",
        type=int,
        default=DEFAULT_MINIMUM_SEQUENCE_SEPARATION,
    )
    parser.add_argument("--max-records", type=int, default=2000)
    args = parser.parse_args()
    output = build_query_centered_hmmer_apc_couplings_v0(
        benchmark_file=Path(args.benchmark_file),
        base_external_coupling_file=Path(args.base_external_coupling_file),
        row_id=args.row_id,
        focus_fasta_file=Path(args.focus_fasta_file),
        output=Path(args.output),
        msa_depth=args.msa_depth,
        effective_sequence_count=args.effective_sequence_count,
        minimum_sequence_separation=args.minimum_sequence_separation,
        max_records=args.max_records,
    )
    print(output)


if __name__ == "__main__":
    main()
