from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


DEFAULT_INPUT_PATH = Path("first_contact_clean_pharmacotopology_layer_run/memory.jsonl")
DEFAULT_RANKINGS_OUTPUT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/pharmacotopology_rankings.csv"
)
DEFAULT_DELTAS_OUTPUT_PATH = Path(
    "first_contact_clean_pharmacotopology_layer_run/pharmacotopology_deltas.csv"
)
REQUIRED_REVIEW_KEY = "Φ.review"


def _json_rows(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parsed = json.loads(line)
        if isinstance(parsed, dict):
            rows.append(parsed)
    return rows


def read_phi_packet(memory_path: Path) -> dict[str, Any]:
    for row in _json_rows(memory_path):
        content = row.get("content")
        if not isinstance(content, str):
            continue
        packet = json.loads(content)
        if isinstance(packet, dict) and REQUIRED_REVIEW_KEY in packet:
            return packet
    raise ValueError(f"No {REQUIRED_REVIEW_KEY} packet found in {memory_path}")


def _join(values: object) -> str:
    if isinstance(values, (list, tuple)):
        return ";".join(str(value) for value in values)
    return str(values)


def _interval_value(interval: Mapping[str, Any], key: str) -> float:
    return float(interval.get(key, 0.0))


def _direction(delta: float) -> str:
    if delta < 0.0:
        return "reduce"
    if delta > 0.0:
        return "worsen"
    return "neutral"


def _result_lookup(review: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        str(result["mechanism_id"]): result
        for result in review.get("Φ.results", [])
        if isinstance(result, dict)
    }


def _ranking_rows(review: Mapping[str, Any]) -> list[dict[str, Any]]:
    results = _result_lookup(review)
    rows: list[dict[str, Any]] = []
    for item in review.get("Φ.ranking", []):
        if not isinstance(item, dict):
            continue
        mechanism_id = str(item["mechanism_id"])
        result = results[mechanism_id]
        pathology_interval = result.get("pathology_reduction_interval", {})
        collapse_interval = result.get("collapse_cost_interval", {})
        net_interval = result.get("net_topology_health_interval", {})
        rows.append(
            {
                "rank": int(item["rank"]),
                "mechanism_id": mechanism_id,
                "mechanism_family": result["mechanism_family"],
                "fit_label": item["fit_label"],
                "pathology_reduction_score": item["pathology_reduction_score"],
                "pathology_reduction_lower": _interval_value(
                    pathology_interval,
                    "lower",
                ),
                "pathology_reduction_upper": _interval_value(
                    pathology_interval,
                    "upper",
                ),
                "collapse_cost_score": item["collapse_cost_score"],
                "collapse_cost_lower": _interval_value(collapse_interval, "lower"),
                "collapse_cost_upper": _interval_value(collapse_interval, "upper"),
                "net_topology_health_score": item["net_topology_health_score"],
                "net_topology_health_lower": _interval_value(net_interval, "lower"),
                "net_topology_health_upper": _interval_value(net_interval, "upper"),
                "evidence_stage": result["evidence_stage"],
                "evidence_weight": item["evidence_weight"],
                "uncertainty_radius": item["uncertainty_radius"],
                "evidence_readiness_label": item["evidence_readiness_label"],
                "primary_evidence_sources": _join(
                    result.get("primary_evidence_sources", ())
                ),
                "evidence_refs": _join(result.get("evidence_refs", ())),
                "calibration_blockers": _join(
                    result.get("calibration_blockers", ())
                ),
                "confidence_interval_kind": result.get(
                    "confidence_interval_kind",
                    "model_uncertainty_interval",
                ),
                "improved_dimensions": _join(result["improved_dimensions"]),
                "worsened_dimensions": _join(result["worsened_dimensions"]),
            }
        )
    return rows


def _delta_rows(review: Mapping[str, Any]) -> list[dict[str, Any]]:
    results = _result_lookup(review)
    source_dimensions = review["Φ.source_profile"]["dimensions"]
    target_dimensions = review["Φ.target_profile"]["dimensions"]
    dimensions = list(source_dimensions)
    rows: list[dict[str, Any]] = []
    rank_lookup = {
        str(item["mechanism_id"]): int(item["rank"])
        for item in review.get("Φ.ranking", [])
        if isinstance(item, dict)
    }

    for mechanism_id, result in results.items():
        topology_delta = result.get("topology_delta", {})
        resulting_state = result.get("resulting_state", {})
        for dimension in dimensions:
            delta = float(topology_delta.get(dimension, 0.0))
            rows.append(
                {
                    "rank": rank_lookup[mechanism_id],
                    "mechanism_id": mechanism_id,
                    "dimension": dimension,
                    "baseline_value": float(source_dimensions[dimension]),
                    "target_value": float(target_dimensions[dimension]),
                    "resulting_value": float(resulting_state.get(dimension, 0.0)),
                    "delta": delta,
                    "direction": _direction(delta),
                }
            )
        collapse_cost = float(result["collapse_cost_score"])
        rows.append(
            {
                "rank": rank_lookup[mechanism_id],
                "mechanism_id": mechanism_id,
                "dimension": "collapse_cost",
                "baseline_value": 0.0,
                "target_value": "",
                "resulting_value": collapse_cost,
                "delta": collapse_cost,
                "direction": "cost",
            }
        )

    return sorted(rows, key=lambda row: (int(row["rank"]), str(row["dimension"])))


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return path
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def export_csv(
    review: Mapping[str, Any],
    rankings_output_path: Path = DEFAULT_RANKINGS_OUTPUT_PATH,
    deltas_output_path: Path = DEFAULT_DELTAS_OUTPUT_PATH,
) -> tuple[Path, Path]:
    rankings_path = _write_csv(rankings_output_path, _ranking_rows(review))
    deltas_path = _write_csv(deltas_output_path, _delta_rows(review))
    return rankings_path, deltas_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export pharmacotopology rankings and deltas to CSV."
    )
    parser.add_argument(
        "--input",
        default=str(DEFAULT_INPUT_PATH),
        help="Path to the pharmacotopology memory.jsonl file.",
    )
    parser.add_argument(
        "--rankings-output",
        default=str(DEFAULT_RANKINGS_OUTPUT_PATH),
        help="Path for the mechanism ranking CSV.",
    )
    parser.add_argument(
        "--deltas-output",
        default=str(DEFAULT_DELTAS_OUTPUT_PATH),
        help="Path for the topology delta CSV.",
    )
    args = parser.parse_args()

    packet = read_phi_packet(Path(args.input))
    rankings_path, deltas_path = export_csv(
        packet[REQUIRED_REVIEW_KEY],
        Path(args.rankings_output),
        Path(args.deltas_output),
    )
    print(rankings_path)
    print(deltas_path)


if __name__ == "__main__":
    main()
