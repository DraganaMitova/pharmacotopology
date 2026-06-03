from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.layer import (  # noqa: E402
    DEFAULT_MECHANISM_VECTORS,
    DEFAULT_NORMAL_BOUNDED_PROFILE,
    DEFAULT_TOPOLOGY_PROFILES,
    PATHOLOGY_DIMENSIONS,
    TOPOLOGY_PRESSURE_MAX,
    TOPOLOGY_PRESSURE_MIN,
    MechanismVector,
    TopologyProfile,
    build_pharmacotopology_review,
    get_topology_profile,
)


DEFAULT_RUN_DIR = Path("first_contact_clean_pharmacotopology_layer_run")
DEFAULT_OUTPUT_PATH = DEFAULT_RUN_DIR / "sensitivity_analysis_report.json"
DEFAULT_CSV_OUTPUT_PATH = DEFAULT_RUN_DIR / "sensitivity_rankings.csv"


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def _shift_profile(
    profile: TopologyProfile,
    dimension: str,
    delta: float,
) -> TopologyProfile:
    dimensions = dict(profile.dimensions)
    dimensions[dimension] = round(
        _clamp(
            dimensions[dimension] + delta,
            TOPOLOGY_PRESSURE_MIN,
            TOPOLOGY_PRESSURE_MAX,
        ),
        6,
    )
    return TopologyProfile(
        profile_id=f"{profile.profile_id}.{dimension}.{delta:+.3f}",
        description=(
            f"Sensitivity variant of {profile.profile_id}; synthetic pressure "
            "shift only, not a diagnosis or patient model."
        ),
        dimensions=dimensions,
    )


def _ranking_map(review: Mapping[str, Any]) -> dict[str, int]:
    return {
        str(item["mechanism_id"]): int(item["rank"])
        for item in review.get("Φ.ranking", [])
        if isinstance(item, dict)
    }


def _ranking_order(review: Mapping[str, Any]) -> tuple[str, ...]:
    return tuple(
        str(item["mechanism_id"])
        for item in review.get("Φ.ranking", [])
        if isinstance(item, dict)
    )


def _top_row(review: Mapping[str, Any]) -> Mapping[str, Any]:
    ranking = review.get("Φ.ranking", [])
    if not ranking:
        return {}
    top = ranking[0]
    return top if isinstance(top, dict) else {}


def _dimension_list(raw_dimensions: Optional[Sequence[str]]) -> tuple[str, ...]:
    if raw_dimensions:
        dimensions = tuple(raw_dimensions)
    else:
        dimensions = PATHOLOGY_DIMENSIONS
    unknown = sorted(set(dimensions).difference(PATHOLOGY_DIMENSIONS))
    if unknown:
        available = ", ".join(PATHOLOGY_DIMENSIONS)
        raise ValueError(f"Unknown dimensions {unknown}. Available: {available}")
    return dimensions


def run_sensitivity_analysis(
    *,
    profile: TopologyProfile,
    target: TopologyProfile = DEFAULT_NORMAL_BOUNDED_PROFILE,
    mechanisms: Sequence[MechanismVector] = DEFAULT_MECHANISM_VECTORS,
    pressure_step: float = 0.08,
    dimensions: Optional[Sequence[str]] = None,
) -> dict[str, Any]:
    selected_dimensions = _dimension_list(dimensions)
    baseline_review = build_pharmacotopology_review(
        source=profile,
        target=target,
        mechanisms=mechanisms,
    )
    baseline_rank = _ranking_map(baseline_review)
    baseline_top = _top_row(baseline_review)
    baseline_top_id = str(baseline_top.get("mechanism_id", ""))
    variants: list[dict[str, Any]] = []

    for dimension in selected_dimensions:
        baseline_value = float(profile.dimensions[dimension])
        for direction, delta in (("down", -pressure_step), ("up", pressure_step)):
            shifted_profile = _shift_profile(profile, dimension, delta)
            review = build_pharmacotopology_review(
                source=shifted_profile,
                target=target,
                mechanisms=mechanisms,
            )
            variant_rank = _ranking_map(review)
            rank_shifts = {
                mechanism_id: variant_rank[mechanism_id] - rank
                for mechanism_id, rank in baseline_rank.items()
            }
            max_abs_rank_shift = max(
                abs(shift) for shift in rank_shifts.values()
            ) if rank_shifts else 0
            top = _top_row(review)
            top_mechanism_id = str(top.get("mechanism_id", ""))
            variants.append(
                {
                    "dimension": dimension,
                    "direction": direction,
                    "delta": round(delta, 6),
                    "baseline_value": baseline_value,
                    "adjusted_value": float(shifted_profile.dimensions[dimension]),
                    "top_mechanism_id": top_mechanism_id,
                    "top_net_topology_health_score": float(
                        top.get("net_topology_health_score", 0.0)
                    ),
                    "top_changed": top_mechanism_id != baseline_top_id,
                    "max_abs_rank_shift": max_abs_rank_shift,
                    "rank_shifts": rank_shifts,
                    "ranking_order": _ranking_order(review),
                }
            )

    return {
        "sensitivity_kind": "local_pressure_sweep",
        "profile": asdict(profile),
        "target_profile": asdict(target),
        "pressure_step": pressure_step,
        "dimensions_reviewed": selected_dimensions,
        "baseline_top_mechanism_id": baseline_top_id,
        "baseline_net_topology_health_score": float(
            baseline_top.get("net_topology_health_score", 0.0)
        ),
        "clinical_use_allowed": False,
        "practical_use": "ranking_robustness_review",
        "variants": variants,
    }


def write_sensitivity_outputs(
    report: Mapping[str, Any],
    output_path: Path = DEFAULT_OUTPUT_PATH,
    csv_output_path: Path = DEFAULT_CSV_OUTPUT_PATH,
) -> tuple[Path, Path]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    csv_output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    csv_rows = [
        {
            "dimension": variant["dimension"],
            "direction": variant["direction"],
            "delta": variant["delta"],
            "baseline_value": variant["baseline_value"],
            "adjusted_value": variant["adjusted_value"],
            "top_mechanism_id": variant["top_mechanism_id"],
            "top_net_topology_health_score": variant[
                "top_net_topology_health_score"
            ],
            "top_changed": variant["top_changed"],
            "max_abs_rank_shift": variant["max_abs_rank_shift"],
            "ranking_order": ";".join(variant["ranking_order"]),
        }
        for variant in report.get("variants", [])
    ]
    with csv_output_path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)

    return output_path, csv_output_path


def _parse_dimensions(value: str) -> tuple[str, ...]:
    if not value.strip():
        return ()
    return tuple(part.strip() for part in value.split(",") if part.strip())


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run a local pressure sensitivity analysis for rankings."
    )
    parser.add_argument(
        "--profile",
        choices=sorted(DEFAULT_TOPOLOGY_PROFILES),
        default="schizophrenia_like",
        help="Synthetic topology profile to sweep.",
    )
    parser.add_argument(
        "--pressure-step",
        type=float,
        default=0.08,
        help="Amount to move each selected topology pressure up and down.",
    )
    parser.add_argument(
        "--dimensions",
        default="",
        help="Comma-separated dimensions to sweep. Defaults to all dimensions.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT_PATH),
        help="Path for the sensitivity JSON report.",
    )
    parser.add_argument(
        "--csv-output",
        default=str(DEFAULT_CSV_OUTPUT_PATH),
        help="Path for the sensitivity CSV summary.",
    )
    args = parser.parse_args()

    report = run_sensitivity_analysis(
        profile=get_topology_profile(args.profile),
        pressure_step=args.pressure_step,
        dimensions=_parse_dimensions(args.dimensions),
    )
    output_path, csv_output_path = write_sensitivity_outputs(
        report,
        Path(args.output),
        Path(args.csv_output),
    )
    print(output_path)
    print(csv_output_path)


if __name__ == "__main__":
    main()
