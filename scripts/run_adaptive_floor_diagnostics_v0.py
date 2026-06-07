from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import mean
from typing import Mapping, Sequence

from pharmacotopology.folding_coupling_nucleus_selector import (
    COUPLING_ADAPTIVE_LOW_SIGNAL_FUTURE_MIN,
    COUPLING_FUTURE_PRESERVATION_MIN,
    COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR_DEFAULT,
    adaptive_coupling_floor_report,
    build_coupling_nucleus_context,
    coupling_nucleus_score,
    _adaptive_coupling_floor_profile,
    _adaptive_future_preservation_threshold,
    _passes_coupling_future_gate,
    _passes_low_signal_adaptive_rescue,
)
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
MANIFEST = REPO_ROOT / "data" / "all_locked_real_external_coupling_holdout_manifest_v0.locked.json"


def _rounded(value: float) -> float:
    return round(value, 6)


def _target_manifest_entries() -> tuple[Mapping[str, object], ...]:
    parsed = json.loads(MANIFEST.read_text(encoding="utf-8"))
    targets = parsed.get("targets", [])
    if not isinstance(targets, list):
        raise ValueError("locked holdout manifest targets must be a list")
    return tuple(target for target in targets if isinstance(target, Mapping))


def _load_rows(path: str) -> tuple[RealCoordinateVisualRow, ...]:
    return tuple(load_real_coordinate_visual_rows(REPO_ROOT / path))


def _row_diagnostics(
    *,
    target_id: str,
    benchmark_file: str,
    external_coupling_file: str,
) -> list[dict[str, object]]:
    rows = _load_rows(benchmark_file)
    dataset = load_coupling_dataset(REPO_ROOT / external_coupling_file)
    context = build_coupling_nucleus_context(rows=rows, coupling_dataset=dataset)
    output: list[dict[str, object]] = []
    for row in rows:
        row_events = tuple(
            event for event in context.competitive_events if event.row_id == row.row_id
        )
        if not row_events:
            continue
        profile = _adaptive_coupling_floor_profile(row_events[0], context)
        row_assessments = [
            context.assessment_by_event_id[event.event_id] for event in row_events
        ]
        adaptive_threshold = _adaptive_future_preservation_threshold(
            row_assessments,
            default=COUPLING_FUTURE_PRESERVATION_MIN,
            floor=profile.future_preservation_floor,
        )
        strict_pass = [
            event
            for event in row_events
            if context.assessment_by_event_id[event.event_id].direct_support_score >= 0.22
            and context.assessment_by_event_id[event.event_id].future_preservation_score
            >= COUPLING_FUTURE_PRESERVATION_MIN
            and context.assessment_by_event_id[event.event_id].blocked_future_pressure
            <= 0.16
        ]
        adaptive_pass = [
            event for event in row_events if _passes_coupling_future_gate(event, context)
        ]
        low_signal_pass = [
            event
            for event in row_events
            if _passes_low_signal_adaptive_rescue(event, context, profile)
        ]
        future_scores = [
            context.assessment_by_event_id[event.event_id].future_preservation_score
            for event in row_events
        ]
        output.append(
            {
                "target_id": target_id,
                "benchmark_file": benchmark_file,
                "external_coupling_file": external_coupling_file,
                "row_id": row.row_id,
                "source_accession": row.source_accession,
                "architecture_axis": row.truth_axes.get("architecture_axis", "unknown"),
                "sequence_length": row.sequence_length,
                "constraint_count": len(dataset.constraints_by_row_id().get(row.row_id, ())),
                **adaptive_coupling_floor_report(row_events[0], context),
                "row_event_count": len(row_events),
                "strict_future_gate_pass_count": len(strict_pass),
                "adaptive_future_gate_pass_count": len(adaptive_pass),
                "low_signal_rescue_pass_count": len(low_signal_pass),
                "adaptive_gap_threshold": _rounded(adaptive_threshold),
                "future_score_min": _rounded(min(future_scores, default=0.0)),
                "future_score_mean": _rounded(mean(future_scores)) if future_scores else 0.0,
                "future_score_max": _rounded(max(future_scores, default=0.0)),
                "low_signal_future_min": COUPLING_ADAPTIVE_LOW_SIGNAL_FUTURE_MIN,
                "old_static_floor": COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR_DEFAULT,
                "best_low_signal_event_score": _rounded(
                    max((coupling_nucleus_score(event, context) for event in low_signal_pass), default=0.0)
                ),
            }
        )
    return output


def _csv_write(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _frontier_1cll_contact_map_precision() -> dict[str, object]:
    benchmark_file = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
    frontier_file = OUTPUT_DIR / "external_coupling_trace_loop_frontier.csv"
    if not frontier_file.exists():
        return {"status": "missing_frontier_file", "perfect_contact_map": False}
    rows = load_real_coordinate_visual_rows(benchmark_file)
    row = next(candidate for candidate in rows if candidate.source_accession == "1CLL:A")
    native_pairs = set(row.native_contact_pairs())
    frontier_rows: list[Mapping[str, str]] = []
    with frontier_file.open(newline="", encoding="utf-8") as handle:
        for item in csv.DictReader(handle):
            if item.get("row_id") == row.row_id:
                frontier_rows.append(item)
    predicted_pairs: set[tuple[int, int]] = set()
    event_pair_counts: list[dict[str, object]] = []
    for item in frontier_rows:
        event_pairs = {
            (left, right)
            for left in range(int(item["segment_a_start"]), int(item["segment_a_end"]) + 1)
            for right in range(int(item["segment_b_start"]), int(item["segment_b_end"]) + 1)
            if right - left >= 3
        }
        predicted_pairs.update(event_pairs)
        true_event_pairs = event_pairs & native_pairs
        event_pair_counts.append(
            {
                "event_id": item.get("event_id", ""),
                "candidate_region_pair_count": len(event_pairs),
                "true_positive_region_pair_count": len(true_event_pairs),
                "region_pair_precision": _rounded(
                    len(true_event_pairs) / len(event_pairs)
                )
                if event_pairs
                else 0.0,
            }
        )
    true_positive_pairs = predicted_pairs & native_pairs
    precision = _rounded(len(true_positive_pairs) / len(predicted_pairs)) if predicted_pairs else 0.0
    recall = _rounded(len(true_positive_pairs) / len(native_pairs)) if native_pairs else 0.0
    return {
        "status": "ok",
        "source_accession": row.source_accession,
        "row_id": row.row_id,
        "frontier_event_count": len(frontier_rows),
        "native_contact_pair_count": len(native_pairs),
        "predicted_region_pair_count": len(predicted_pairs),
        "true_positive_region_pair_count": len(true_positive_pairs),
        "false_positive_region_pair_count": len(predicted_pairs - native_pairs),
        "false_negative_native_pair_count": len(native_pairs - predicted_pairs),
        "region_contact_precision": precision,
        "region_contact_recall": recall,
        "perfect_contact_map": bool(precision == 1.0 and recall == 1.0),
        "interpretation": (
            "7 frontier events is not a perfect contact map unless precision and recall are both 1.0"
        ),
        "event_region_precision_rows": event_pair_counts,
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    row_reports: list[dict[str, object]] = []
    for target in _target_manifest_entries():
        benchmark_file = str(target["benchmark_file"])
        external_coupling_file = str(target["external_coupling_file"])
        row_reports.extend(
            _row_diagnostics(
                target_id=str(target["target_id"]),
                benchmark_file=benchmark_file,
                external_coupling_file=external_coupling_file,
            )
        )

    expansion_summary = {
        "locked_target_count": len(_target_manifest_entries()),
        "locked_row_count": len(row_reports),
        "available_multidomain_or_segmented_rows": [
            row["source_accession"]
            for row in row_reports
            if row.get("architecture_axis") == "multidomain_or_segmented"
        ],
        "unknown_architecture_rows": [
            row["source_accession"]
            for row in row_reports
            if row.get("architecture_axis") == "unknown"
        ],
        "note": (
            "This scans every complete locked external-coupling target present in the archive. "
            "It does not fetch new proteins from the internet."
        ),
    }
    contact_precision = _frontier_1cll_contact_map_precision()
    report = {
        "report_kind": "adaptive_coupling_floor_diagnostics_v0",
        "coordinate_truth_used_to_build_constraints": False,
        "native_truth_used_before_coupling_selection": False,
        "raw_sequence_exposed": False,
        "old_static_future_floor": COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR_DEFAULT,
        "adaptive_floor_inputs": [
            "phase_mode",
            "sequence_complexity",
            "effective_sequence_count_over_length",
            "target_coverage",
            "row_future_preservation_ceiling_noise_guard",
        ],
        "dataset_expansion_summary": expansion_summary,
        "one_cll_contact_map_precision_check": contact_precision,
        "row_reports": row_reports,
    }
    (OUTPUT_DIR / "adaptive_floor_diagnostics_v0.json").write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _csv_write(OUTPUT_DIR / "adaptive_floor_diagnostics_rows_v0.csv", row_reports)
    (OUTPUT_DIR / "1cll_contact_map_precision_v0.json").write_text(
        json.dumps(contact_precision, indent=2, sort_keys=True),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
