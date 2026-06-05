from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_BASE_EXTERNAL_COUPLING_FILE = Path(
    "data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
)


@dataclass(frozen=True)
class PlmcPair:
    score: float
    i: int
    j: int
    alignment_column_i: int
    alignment_column_j: int


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(block)
    return f"sha256:{digest.hexdigest()}"


def _read_focus_fasta(path: Path) -> tuple[str, str]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    if len(lines) < 2 or not lines[0].startswith(">"):
        raise ValueError(f"{path} must contain a focus FASTA record first")
    sequence_lines: list[str] = []
    for line in lines[1:]:
        if line.startswith(">"):
            break
        sequence_lines.append(line)
    return lines[0][1:], "".join(sequence_lines)


def _alignment_column_to_row_position(
    *,
    focus_alignment: str,
    row: RealCoordinateVisualRow,
) -> dict[int, int]:
    letters = "".join(char.upper() for char in focus_alignment if char.isalpha())
    if letters != row.sequence:
        raise ValueError(
            f"focus sequence does not match frozen row sequence for {row.row_id}"
        )
    column_to_position: dict[int, int] = {}
    row_position = 0
    for alignment_column, char in enumerate(focus_alignment, start=1):
        if not char.isalpha():
            continue
        row_position += 1
        if char.isupper():
            column_to_position[alignment_column] = row_position
    if row_position != row.sequence_length:
        raise ValueError(
            f"focus sequence length {row_position} does not match "
            f"{row.sequence_length} for {row.row_id}"
        )
    return column_to_position


def _iter_plmc_pairs(
    *,
    plmc_couplings_file: Path,
    column_to_position: dict[int, int],
) -> Iterable[PlmcPair]:
    for line in plmc_couplings_file.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        column_i = int(parts[0])
        column_j = int(parts[2])
        if column_i not in column_to_position or column_j not in column_to_position:
            continue
        i = column_to_position[column_i]
        j = column_to_position[column_j]
        if j < i:
            i, j = j, i
            column_i, column_j = column_j, column_i
        yield PlmcPair(
            score=float(parts[-1]),
            i=i,
            j=j,
            alignment_column_i=column_i,
            alignment_column_j=column_j,
        )


def _rank_unique_pairs(pairs: Iterable[PlmcPair]) -> list[PlmcPair]:
    seen: set[tuple[int, int]] = set()
    ranked: list[PlmcPair] = []
    for pair in sorted(pairs, key=lambda item: (-item.score, item.i, item.j)):
        key = (pair.i, pair.j)
        if key in seen:
            continue
        seen.add(key)
        ranked.append(pair)
    return ranked


def _constraint_rows(
    *,
    row: RealCoordinateVisualRow,
    pairs: list[PlmcPair],
    msa_sha256: str,
    msa_depth: int,
    effective_sequence_count: float,
    focus_id: str,
    hmmer_job_id: str,
    hmmer_iteration_id: str,
) -> list[dict[str, object]]:
    selected = pairs[: row.sequence_length]
    max_positive_score = max((max(0.0, pair.score) for pair in selected), default=0.0)
    constraints: list[dict[str, object]] = []
    for rank, pair in enumerate(selected, start=1):
        confidence = (
            max(0.01, max(0.0, pair.score) / max_positive_score)
            if max_positive_score
            else 0.01
        )
        constraints.append(
            {
                "row_id": row.row_id,
                "source_accession": row.source_accession,
                "constraint_id": (
                    f"query_hmmer_plmc_{row.row_id}_{pair.i}_{pair.j}_{rank}"
                ),
                "i": pair.i,
                "j": pair.j,
                "sequence_separation": pair.j - pair.i,
                "normalized_separation": round(
                    (pair.j - pair.i) / row.sequence_length,
                    6,
                ),
                "confidence": round(confidence, 6),
                "raw_score": round(pair.score, 6),
                "apc_corrected_score": round(pair.score, 6),
                "rank": rank,
                "rank_fraction": round(rank / row.sequence_length, 6),
                "constraint_class": "external_query_centered_hmmer_plmc_coupling",
                "source_kind": "external_msa_dca_plmc_v1",
                "msa_source_kind": (
                    "ebi_hmmer_jackhmmer_uniprot_iter3_stockholm_focus_plmc_"
                    "unweighted"
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
                "hmmer_job_id": hmmer_job_id,
                "hmmer_iteration_id": hmmer_iteration_id,
                "plmc_focus_id": focus_id,
                "plmc_alignment_column_i": pair.alignment_column_i,
                "plmc_alignment_column_j": pair.alignment_column_j,
            }
        )
    return constraints


def build_query_centered_hmmer_plmc_couplings_v0(
    *,
    benchmark_file: Path,
    base_external_coupling_file: Path,
    row_id: str,
    focus_fasta_file: Path,
    plmc_couplings_file: Path,
    output: Path,
    msa_depth: int,
    effective_sequence_count: float,
    hmmer_job_id: str,
    hmmer_iteration_id: str,
) -> Path:
    rows = tuple(load_real_coordinate_visual_rows(benchmark_file))
    row_by_id = {row.row_id: row for row in rows}
    row = row_by_id[row_id]
    focus_id, focus_alignment = _read_focus_fasta(focus_fasta_file)
    column_to_position = _alignment_column_to_row_position(
        focus_alignment=focus_alignment,
        row=row,
    )
    pairs = _rank_unique_pairs(
        _iter_plmc_pairs(
            plmc_couplings_file=plmc_couplings_file,
            column_to_position=column_to_position,
        )
    )
    if len(pairs) < row.sequence_length:
        raise ValueError(
            f"{row_id} has only {len(pairs)} unique mapped plmc pairs; "
            f"{row.sequence_length} required"
        )
    constraints = _constraint_rows(
        row=row,
        pairs=pairs,
        msa_sha256=_sha256_file(focus_fasta_file),
        msa_depth=msa_depth,
        effective_sequence_count=effective_sequence_count,
        focus_id=focus_id,
        hmmer_job_id=hmmer_job_id,
        hmmer_iteration_id=hmmer_iteration_id,
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
    payload[f"hmmer_query_centered_{row.source_accession.replace(':', '_')}_added"] = True
    payload[
        f"hmmer_query_centered_{row.source_accession.replace(':', '_')}_alignment_sha256"
    ] = _sha256_file(focus_fasta_file)
    payload[
        f"hmmer_query_centered_{row.source_accession.replace(':', '_')}_valid_sequence_count"
    ] = msa_depth
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
            "query-centered EBI HMMER jackhmmer + plmc couplings."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument(
        "--base-external-coupling-file",
        default=str(DEFAULT_BASE_EXTERNAL_COUPLING_FILE),
    )
    parser.add_argument("--row-id", required=True)
    parser.add_argument("--focus-fasta-file", required=True)
    parser.add_argument("--plmc-couplings-file", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--msa-depth", type=int, required=True)
    parser.add_argument("--effective-sequence-count", type=float, required=True)
    parser.add_argument("--hmmer-job-id", required=True)
    parser.add_argument("--hmmer-iteration-id", required=True)
    args = parser.parse_args()
    output = build_query_centered_hmmer_plmc_couplings_v0(
        benchmark_file=Path(args.benchmark_file),
        base_external_coupling_file=Path(args.base_external_coupling_file),
        row_id=args.row_id,
        focus_fasta_file=Path(args.focus_fasta_file),
        plmc_couplings_file=Path(args.plmc_couplings_file),
        output=Path(args.output),
        msa_depth=args.msa_depth,
        effective_sequence_count=args.effective_sequence_count,
        hmmer_job_id=args.hmmer_job_id,
        hmmer_iteration_id=args.hmmer_iteration_id,
    )
    print(output)


if __name__ == "__main__":
    main()
