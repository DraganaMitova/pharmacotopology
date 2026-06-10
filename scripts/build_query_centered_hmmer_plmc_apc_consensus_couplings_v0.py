from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from build_query_centered_hmmer_apc_couplings_v0 import (  # noqa: E402
    _query_centered_apc_pairs,
)
from build_query_centered_hmmer_plmc_couplings_v0 import (  # noqa: E402
    DEFAULT_BASE_EXTERNAL_COUPLING_FILE,
    DEFAULT_BENCHMARK_FILE,
    DEFAULT_MINIMUM_SEQUENCE_SEPARATION,
    _read_focus_fasta,
    _sha256_file,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


DEFAULT_CONSENSUS_WEIGHT = 0.75


@dataclass(frozen=True)
class ConsensusPairScore:
    score: float
    plmc_confidence: float
    apc_confidence: float
    constraint: dict[str, object]


def _apc_confidence_by_pair(
    *,
    row: RealCoordinateVisualRow,
    focus_fasta_file: Path,
    minimum_sequence_separation: int,
    max_records: int,
) -> dict[tuple[int, int], float]:
    pairs = _query_centered_apc_pairs(
        row=row,
        focus_fasta_file=focus_fasta_file,
        minimum_sequence_separation=minimum_sequence_separation,
        max_records=max_records,
    )
    selected = tuple(pairs[: row.sequence_length])
    max_positive_score = max(
        (max(0.0, pair.apc_corrected_score) for pair in selected),
        default=0.0,
    )
    if not max_positive_score:
        return {}
    return {
        (pair.i, pair.j): max(
            0.01,
            max(0.0, pair.apc_corrected_score) / max_positive_score,
        )
        for pair in selected
    }


def _consensus_pair_scores(
    *,
    base_constraints: list[dict[str, object]],
    apc_confidence_by_pair: dict[tuple[int, int], float],
    consensus_weight: float,
) -> list[ConsensusPairScore]:
    scores: list[ConsensusPairScore] = []
    for constraint in base_constraints:
        pair = (int(constraint["i"]), int(constraint["j"]))
        plmc_confidence = float(constraint["confidence"])
        apc_confidence = apc_confidence_by_pair.get(pair, 0.0)
        consensus_score = plmc_confidence + consensus_weight * min(
            plmc_confidence,
            apc_confidence,
        )
        scores.append(
            ConsensusPairScore(
                score=consensus_score,
                plmc_confidence=plmc_confidence,
                apc_confidence=apc_confidence,
                constraint=constraint,
            )
        )
    return sorted(
        scores,
        key=lambda item: (
            -item.score,
            int(item.constraint["i"]),
            int(item.constraint["j"]),
        ),
    )


def _constraint_rows(
    *,
    row: RealCoordinateVisualRow,
    pair_scores: list[ConsensusPairScore],
    focus_id: str,
    msa_sha256: str,
    consensus_weight: float,
    minimum_sequence_separation: int,
) -> list[dict[str, object]]:
    selected = pair_scores[: row.sequence_length]
    max_score = max((item.score for item in selected), default=0.0)
    if not max_score:
        raise ValueError(f"{row.row_id} has no positive PLMC/APC consensus score")
    constraints: list[dict[str, object]] = []
    for rank, item in enumerate(selected, start=1):
        raw = dict(item.constraint)
        raw.update(
            {
                "constraint_id": (
                    f"query_hmmer_plmc_apc_consensus_{row.row_id}_"
                    f"{raw['i']}_{raw['j']}_{rank}"
                ),
                "confidence": round(max(0.01, item.score / max_score), 6),
                "raw_score": round(item.score, 6),
                "apc_corrected_score": round(item.score, 6),
                "rank": rank,
                "rank_fraction": round(rank / row.sequence_length, 6),
                "constraint_class": (
                    "external_query_centered_hmmer_plmc_apc_consensus_coupling"
                ),
                "msa_source_kind": (
                    "ebi_hmmer_jackhmmer_uniprot_focus_plmc_apc_consensus"
                ),
                "coordinate_truth_used_to_build_constraint": False,
                "native_truth_used_before_coupling_selection": False,
                "structure_model_used": False,
                "raw_sequence_exposed": False,
                "plmc_apc_consensus_focus_id": focus_id,
                "plmc_apc_consensus_alignment_sha256": msa_sha256,
                "plmc_apc_consensus_weight": consensus_weight,
                "plmc_apc_consensus_method": "plmc_plus_weighted_apc_agreement",
                "fusion_plmc_confidence": round(item.plmc_confidence, 6),
                "fusion_apc_confidence": round(item.apc_confidence, 6),
                "minimum_sequence_separation": minimum_sequence_separation,
                "top_L_couplings_available": len(pair_scores) >= row.sequence_length,
            }
        )
        constraints.append(raw)
    return constraints


def build_query_centered_hmmer_plmc_apc_consensus_couplings_v0(
    *,
    benchmark_file: Path,
    base_external_coupling_file: Path,
    row_id: str,
    focus_fasta_file: Path,
    output: Path,
    consensus_weight: float = DEFAULT_CONSENSUS_WEIGHT,
    minimum_sequence_separation: int = DEFAULT_MINIMUM_SEQUENCE_SEPARATION,
    max_records: int = 2000,
) -> Path:
    if consensus_weight < 0.0:
        raise ValueError("consensus_weight must be non-negative")
    rows = tuple(load_real_coordinate_visual_rows(benchmark_file))
    row_by_id = {row.row_id: row for row in rows}
    row = row_by_id[row_id]
    focus_id, _ = _read_focus_fasta(focus_fasta_file)
    payload = json.loads(base_external_coupling_file.read_text(encoding="utf-8"))
    base_constraints = [
        constraint
        for constraint in payload["constraints"]
        if constraint["row_id"] == row_id
    ]
    if len(base_constraints) < row.sequence_length:
        raise ValueError(
            f"{row_id} has only {len(base_constraints)} base constraints; "
            f"{row.sequence_length} required"
        )
    if any(
        int(constraint["sequence_separation"]) < minimum_sequence_separation
        for constraint in base_constraints[: row.sequence_length]
    ):
        raise ValueError(
            f"{row_id} base top-L constraints are not all separated by "
            f">= {minimum_sequence_separation}"
        )
    apc_confidence = _apc_confidence_by_pair(
        row=row,
        focus_fasta_file=focus_fasta_file,
        minimum_sequence_separation=minimum_sequence_separation,
        max_records=max_records,
    )
    pair_scores = _consensus_pair_scores(
        base_constraints=base_constraints,
        apc_confidence_by_pair=apc_confidence,
        consensus_weight=consensus_weight,
    )
    constraints = _constraint_rows(
        row=row,
        pair_scores=pair_scores,
        focus_id=focus_id,
        msa_sha256=_sha256_file(focus_fasta_file),
        consensus_weight=consensus_weight,
        minimum_sequence_separation=minimum_sequence_separation,
    )
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
    source_slug = row.source_accession.replace(":", "_")
    payload["external_constraint_count"] = len(payload["constraints"])
    payload[f"hmmer_query_centered_plmc_apc_consensus_{source_slug}_added"] = True
    payload[
        f"hmmer_query_centered_plmc_apc_consensus_{source_slug}_alignment_sha256"
    ] = _sha256_file(focus_fasta_file)
    payload[
        f"hmmer_query_centered_plmc_apc_consensus_{source_slug}_weight"
    ] = consensus_weight
    payload[
        f"hmmer_query_centered_plmc_apc_consensus_{source_slug}_minimum_sequence_separation"
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
            "Reweight one locked query-centered HMMER PLMC row by independent "
            "query-centered HMMER MI/APC agreement."
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
    parser.add_argument(
        "--consensus-weight",
        type=float,
        default=DEFAULT_CONSENSUS_WEIGHT,
    )
    parser.add_argument(
        "--minimum-sequence-separation",
        type=int,
        default=DEFAULT_MINIMUM_SEQUENCE_SEPARATION,
    )
    parser.add_argument("--max-records", type=int, default=2000)
    args = parser.parse_args()
    output = build_query_centered_hmmer_plmc_apc_consensus_couplings_v0(
        benchmark_file=Path(args.benchmark_file),
        base_external_coupling_file=Path(args.base_external_coupling_file),
        row_id=args.row_id,
        focus_fasta_file=Path(args.focus_fasta_file),
        output=Path(args.output),
        consensus_weight=args.consensus_weight,
        minimum_sequence_separation=args.minimum_sequence_separation,
        max_records=args.max_records,
    )
    print(output)


if __name__ == "__main__":
    main()
