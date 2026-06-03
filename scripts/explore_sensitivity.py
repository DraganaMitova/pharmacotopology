from __future__ import annotations

import argparse
import csv
import json
import random
import statistics
import sys
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any, Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.layer import (  # noqa: E402
    DEFAULT_MECHANISM_VECTORS,
    DEFAULT_TOPOLOGY_PROFILES,
    PATHOLOGY_DIMENSIONS,
    MechanismVector,
    build_pharmacotopology_review,
    get_topology_profile,
)


DEFAULT_OUTPUT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/sensitivity_explorer_report.json"
)
DEFAULT_SAMPLES_OUTPUT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/sensitivity_explorer_samples.csv"
)


def _rounded(value: float) -> float:
    return round(value, 6)


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def _parse_range(value: str) -> tuple[float, float]:
    if ":" not in value:
        raise ValueError("Range must use lower:upper, for example 0.05:0.25")
    lower, upper = (float(part.strip()) for part in value.split(":", 1))
    if lower > upper:
        raise ValueError("Range lower bound cannot be greater than upper bound")
    return lower, upper


def _normalize_vary_field(value: str) -> tuple[str, str]:
    if value == "collapse_cost":
        return "collapse_cost", ""
    if value.startswith("delta:"):
        dimension = value.split(":", 1)[1]
    elif value.endswith("_delta"):
        dimension = value[: -len("_delta")]
    else:
        dimension = value
    if dimension not in PATHOLOGY_DIMENSIONS:
        available = ", ".join(("collapse_cost",) + PATHOLOGY_DIMENSIONS)
        raise ValueError(f"Unknown vary field {value!r}. Available: {available}")
    return "delta", dimension


def _mechanism_lookup(
    mechanisms: Sequence[MechanismVector],
) -> dict[str, MechanismVector]:
    return {mechanism.mechanism_id: mechanism for mechanism in mechanisms}


def _perturb_vector(
    vector: MechanismVector,
    rng: random.Random,
    *,
    noise: float,
    target_mechanism_id: str,
    vary_kind: str,
    vary_dimension: str,
    vary_range: tuple[float, float],
) -> MechanismVector:
    collapse_cost = _clamp(
        vector.collapse_cost + rng.uniform(-noise, noise),
        0.0,
        1.0,
    )
    deltas = {
        dimension: _rounded(float(delta) + rng.uniform(-noise, noise))
        for dimension, delta in vector.deltas.items()
    }

    if vector.mechanism_id == target_mechanism_id:
        varied_value = rng.uniform(vary_range[0], vary_range[1])
        if vary_kind == "collapse_cost":
            collapse_cost = _clamp(varied_value, 0.0, 1.0)
        else:
            deltas[vary_dimension] = _rounded(varied_value)

    return replace(vector, collapse_cost=_rounded(collapse_cost), deltas=deltas)


def _ranking_map(review: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        str(item["mechanism_id"]): item
        for item in review.get("Φ.ranking", [])
        if isinstance(item, dict)
    }


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _percentile(values: Sequence[float], percentile: float) -> float:
    if not values:
        return 0.0
    sorted_values = sorted(values)
    index = int(round((len(sorted_values) - 1) * percentile))
    return sorted_values[index]


def explore_sensitivity(
    *,
    profile_key: str,
    mechanism_id: str,
    vary_field: str,
    vary_range: tuple[float, float],
    samples: int = 100,
    noise: float = 0.03,
    seed: int = 1729,
    mechanisms: Sequence[MechanismVector] = DEFAULT_MECHANISM_VECTORS,
) -> dict[str, Any]:
    if samples <= 0:
        raise ValueError("samples must be positive")
    if noise < 0:
        raise ValueError("noise must be non-negative")

    source = get_topology_profile(profile_key)
    mechanism_ids = tuple(mechanism.mechanism_id for mechanism in mechanisms)
    if mechanism_id not in mechanism_ids:
        available = ", ".join(mechanism_ids)
        raise ValueError(f"Unknown mechanism {mechanism_id!r}. Available: {available}")
    vary_kind, vary_dimension = _normalize_vary_field(vary_field)

    baseline_review = build_pharmacotopology_review(
        source=source,
        mechanisms=mechanisms,
    )
    baseline_ranking = _ranking_map(baseline_review)
    baseline_target = baseline_ranking[mechanism_id]
    rng = random.Random(seed)
    sample_rows: list[dict[str, Any]] = []
    per_mechanism: dict[str, dict[str, list[float]]] = {
        item: {"ranks": [], "scores": []}
        for item in mechanism_ids
    }

    for sample_id in range(1, samples + 1):
        perturbed = tuple(
            _perturb_vector(
                vector,
                rng,
                noise=noise,
                target_mechanism_id=mechanism_id,
                vary_kind=vary_kind,
                vary_dimension=vary_dimension,
                vary_range=vary_range,
            )
            for vector in mechanisms
        )
        review = build_pharmacotopology_review(source=source, mechanisms=perturbed)
        ranking = _ranking_map(review)
        target_row = ranking[mechanism_id]
        varied_vector = _mechanism_lookup(perturbed)[mechanism_id]
        varied_value = (
            varied_vector.collapse_cost
            if vary_kind == "collapse_cost"
            else varied_vector.deltas[vary_dimension]
        )
        sample_rows.append(
            {
                "sample_id": sample_id,
                "profile": profile_key,
                "mechanism_id": mechanism_id,
                "vary_field": vary_field,
                "varied_value": _rounded(float(varied_value)),
                "rank": int(target_row["rank"]),
                "net_topology_health_score": float(
                    target_row["net_topology_health_score"]
                ),
                "pathology_reduction_score": float(
                    target_row["pathology_reduction_score"]
                ),
                "collapse_cost_score": float(target_row["collapse_cost_score"]),
                "top_mechanism_id": str(review["Φ.ranking"][0]["mechanism_id"]),
            }
        )
        for item_id, item in ranking.items():
            per_mechanism[item_id]["ranks"].append(float(item["rank"]))
            per_mechanism[item_id]["scores"].append(
                float(item["net_topology_health_score"])
            )

    target_scores = [
        float(row["net_topology_health_score"])
        for row in sample_rows
    ]
    target_ranks = [int(row["rank"]) for row in sample_rows]
    robustness_rows = []
    for item_id, values in per_mechanism.items():
        ranks = values["ranks"]
        scores = values["scores"]
        rank_std = statistics.pstdev(ranks) if len(ranks) > 1 else 0.0
        score_std = statistics.pstdev(scores) if len(scores) > 1 else 0.0
        baseline_rank = int(baseline_ranking[item_id]["rank"])
        changed = sum(1 for rank in ranks if int(rank) != baseline_rank)
        robustness_score = _clamp(
            1.0
            - (rank_std / max(len(mechanism_ids) - 1, 1))
            - min(score_std, 1.0),
            0.0,
            1.0,
        )
        robustness_rows.append(
            {
                "mechanism_id": item_id,
                "baseline_rank": baseline_rank,
                "mean_rank": _rounded(_mean(ranks)),
                "rank_std": _rounded(rank_std),
                "rank_change_rate": _rounded(changed / len(ranks)),
                "mean_net_topology_health_score": _rounded(_mean(scores)),
                "score_std": _rounded(score_std),
                "robustness_score": _rounded(robustness_score),
            }
        )

    return {
        "explorer_kind": "monte_carlo_assumption_sensitivity",
        "profile": profile_key,
        "mechanism_id": mechanism_id,
        "vary_field": vary_field,
        "vary_range": {
            "lower": vary_range[0],
            "upper": vary_range[1],
        },
        "samples": samples,
        "noise": noise,
        "seed": seed,
        "clinical_use_allowed": False,
        "practical_use": "assumption_sensitivity_exploration",
        "baseline": {
            "rank": int(baseline_target["rank"]),
            "net_topology_health_score": float(
                baseline_target["net_topology_health_score"]
            ),
            "pathology_reduction_score": float(
                baseline_target["pathology_reduction_score"]
            ),
            "collapse_cost_score": float(baseline_target["collapse_cost_score"]),
        },
        "target_distribution": {
            "mean_net_topology_health_score": _rounded(_mean(target_scores)),
            "net_score_std": _rounded(
                statistics.pstdev(target_scores) if len(target_scores) > 1 else 0.0
            ),
            "net_score_p05": _rounded(_percentile(target_scores, 0.05)),
            "net_score_p95": _rounded(_percentile(target_scores, 0.95)),
            "mean_rank": _rounded(_mean([float(rank) for rank in target_ranks])),
            "rank_change_rate": _rounded(
                sum(
                    1
                    for rank in target_ranks
                    if rank != int(baseline_target["rank"])
                )
                / len(target_ranks)
            ),
        },
        "robustness": sorted(
            robustness_rows,
            key=lambda row: float(row["robustness_score"]),
            reverse=True,
        ),
        "samples_table": sample_rows,
    }


def write_explorer_outputs(
    report: Mapping[str, Any],
    output_path: Path = DEFAULT_OUTPUT_PATH,
    samples_output_path: Path = DEFAULT_SAMPLES_OUTPUT_PATH,
) -> tuple[Path, Path]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    samples_output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    sample_rows = list(report.get("samples_table", []))
    with samples_output_path.open("w", encoding="utf-8", newline="") as file:
        if sample_rows:
            writer = csv.DictWriter(file, fieldnames=list(sample_rows[0]))
            writer.writeheader()
            writer.writerows(sample_rows)
    return output_path, samples_output_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Explore pharmacotopology sensitivity with deterministic Monte Carlo sampling."
    )
    parser.add_argument(
        "--profile",
        choices=sorted(DEFAULT_TOPOLOGY_PROFILES),
        default="schizophrenia_like",
        help="Synthetic topology profile to explore.",
    )
    parser.add_argument(
        "--mechanism",
        default="nmda_support_like",
        help="Mechanism vector to vary.",
    )
    parser.add_argument(
        "--vary",
        nargs=2,
        metavar=("FIELD", "LOWER:UPPER"),
        default=("collapse_cost", "0.05:0.25"),
        help=(
            "Field and range to vary. Use collapse_cost, delta:<dimension>, "
            "or <dimension>_delta."
        ),
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=100,
        help="Number of deterministic Monte Carlo samples.",
    )
    parser.add_argument(
        "--noise",
        type=float,
        default=0.03,
        help="Uniform noise applied to all mechanism deltas and collapse costs.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1729,
        help="Random seed for reproducible exploration.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT_PATH),
        help="Path for the sensitivity explorer JSON report.",
    )
    parser.add_argument(
        "--samples-output",
        default=str(DEFAULT_SAMPLES_OUTPUT_PATH),
        help="Path for the sensitivity explorer sample CSV.",
    )
    args = parser.parse_args()

    report = explore_sensitivity(
        profile_key=args.profile,
        mechanism_id=args.mechanism,
        vary_field=args.vary[0],
        vary_range=_parse_range(args.vary[1]),
        samples=args.samples,
        noise=args.noise,
        seed=args.seed,
    )
    output_path, samples_output_path = write_explorer_outputs(
        report,
        Path(args.output),
        Path(args.samples_output),
    )
    print(output_path)
    print(samples_output_path)


if __name__ == "__main__":
    main()
