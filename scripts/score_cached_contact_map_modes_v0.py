from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
for root in (REPO_ROOT, SRC_ROOT):
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))

from pharmacotopology.artifact_io import write_csv_rows  # noqa: E402
from pharmacotopology.folding_coupling_nucleus_selector import (  # noqa: E402
    CouplingSelectorMetric,
)
from pharmacotopology.folding_external_coupling_trace_loop import (  # noqa: E402
    TraceLoopRun,
)
from pharmacotopology.folding_nucleus_closure_search import (  # noqa: E402
    NucleusClosureEvent,
)

from scripts.run_blind_external_holdout_battery_v0 import (  # noqa: E402
    BATTERY_REPORT_KIND,
    DEFAULT_COUPLING_DIR,
    REPO_ROOT as BATTERY_REPO_ROOT,
    SATURATION_BOUNDARY_SELECTOR,
    SELF_CRITICAL_QUALITY_SWITCH_SELECTOR,
    _apply_self_critical_quality_switch,
    _exact_contact_metrics,
    _external_coupling_quality_summary,
    _load_json,
    _mean_field,
    _resolve_path,
    _scaffold_contact_metrics,
    _target_coupling_file,
    import_external_coupling_dataset,
    load_real_coordinate_visual_rows,
)


CACHED_MODE_SCORE_REPORT_KIND = "cached_contact_map_mode_score_v0"
DEFAULT_TARGET_MANIFEST = Path(
    "data/all_locked_real_external_coupling_holdout_manifest_v0.locked.json"
)
DEFAULT_SELECTED_EVENTS = Path(
    "first_contact_clean_pharmacotopology_layer_run/"
    "all_locked_real_external_holdout_v0/blind_external_holdout_selected_events.csv"
)
DEFAULT_TARGET_ROWS = Path(
    "first_contact_clean_pharmacotopology_layer_run/"
    "all_locked_real_external_holdout_v0/blind_external_holdout_rows.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    "first_contact_clean_pharmacotopology_layer_run/"
    "cached_contact_map_mode_score_v0"
)


def _read_csv_rows(path: Path) -> tuple[dict[str, str], ...]:
    if not path.exists() or path.stat().st_size == 0:
        return ()
    with path.open(newline="", encoding="utf-8") as file:
        return tuple(dict(row) for row in csv.DictReader(file))


def _float(raw: object, default: float = 0.0) -> float:
    try:
        return float(raw)
    except (TypeError, ValueError):
        return default


def _int(raw: object, default: int = 0) -> int:
    try:
        return int(float(raw))
    except (TypeError, ValueError):
        return default


def _cached_event(
    raw: Mapping[str, str],
    *,
    row_by_id: Mapping[str, object],
) -> NucleusClosureEvent:
    row = row_by_id[str(raw["row_id"])]
    segment_a_start = _int(raw["segment_a_start"])
    segment_a_end = _int(raw["segment_a_end"])
    segment_b_start = _int(raw["segment_b_start"])
    segment_b_end = _int(raw["segment_b_end"])
    candidate_contact_count = sum(
        1
        for left in range(segment_a_start, segment_a_end + 1)
        for right in range(segment_b_start, segment_b_end + 1)
        if right - left >= 3
    )
    return NucleusClosureEvent(
        row_id=str(raw["row_id"]),
        source_accession=str(raw.get("source_accession", row.source_accession)),
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        event_id=str(raw["event_id"]),
        segment_a_start=segment_a_start,
        segment_a_end=segment_a_end,
        segment_b_start=segment_b_start,
        segment_b_end=segment_b_end,
        sequence_span=segment_b_end - segment_a_start + 1,
        normalized_span=_float(raw.get("normalized_span")),
        candidate_contact_count=candidate_contact_count,
        contact_cluster_gain=_float(raw.get("coupling_selectivity_score")),
        secondary_structure_compatibility=0.0,
        hydrophobic_burial_gain=_float(raw.get("burial_gain")),
        registry_support=_float(raw.get("direct_support_score")),
        loop_entropy_cost=0.0,
        geometry_violation_cost=0.0,
        frustration_cost=_float(raw.get("unsatisfied_polar_penalty")),
        isolation_penalty=0.0,
        nucleus_score=_float(raw.get("coupling_nucleus_score")),
        closure_event_stability=_float(
            raw.get("multiscale_critical_coherence_score"),
            _float(raw.get("coupling_nucleus_score")),
        ),
        native_contact_count_after_scoring=_int(
            raw.get("native_contact_count_after_scoring")
        ),
        native_long_range_contact_count_after_scoring=_int(
            raw.get("native_contact_count_after_scoring")
        ),
        native_label_attached_after_event_generation=True,
        native_truth_used_before_event_generation=False,
        raw_sequence_exposed=str(raw.get("raw_sequence_exposed", "False")) == "True",
    )


def _cached_metric(
    *,
    target_id: str,
    target_row: Mapping[str, str] | None,
    selected_event_count: int,
) -> CouplingSelectorMetric:
    source = target_row or {}
    return CouplingSelectorMetric(
        selector_name=SATURATION_BOUNDARY_SELECTOR,
        selected_event_count=selected_event_count,
        false_nucleus_rate=_float(source.get("false_nucleus_rate")),
        contact_cluster_precision=_float(source.get("cluster_precision")),
        long_range_contact_recall=_float(source.get("long_range_contact_recall")),
        coupling_constraint_recall=_float(source.get("coupling_constraint_recall")),
        real_vs_decoy_coupling_enrichment_ratio=_float(
            source.get("real_vs_decoy_coupling_enrichment_ratio")
        ),
        real_beats_decoy_coupling_score_rate=0.0,
        mean_selected_coupling_selectivity_score=0.0,
        mean_decoy_coupling_selectivity_score=0.0,
        mean_coupling_decoy_selectivity_margin=0.0,
        mean_coupling_nucleus_score=0.0,
        mean_decoy_coupling_nucleus_score=0.0,
        mean_coupling_nucleus_decoy_margin=0.0,
        real_vs_decoy_coupling_nucleus_enrichment_ratio=0.0,
        real_beats_decoy_coupling_nucleus_score_rate=0.0,
        survives_targets=selected_event_count > 0,
        coordinate_truth_used_to_build_constraints=False,
        native_truth_used_before_coupling_selection=False,
        raw_sequence_exposed=False,
    )


def _events_for_target(
    selected_rows: Sequence[Mapping[str, str]],
    *,
    target_id: str,
    row_by_id: Mapping[str, object],
    selector_name: str,
) -> tuple[NucleusClosureEvent, ...]:
    seen: set[tuple[str, str]] = set()
    events: list[NucleusClosureEvent] = []
    for raw in selected_rows:
        if str(raw.get("target_id", "")) != target_id:
            continue
        if str(raw.get("selector_name", "")) != selector_name:
            continue
        key = (str(raw.get("row_id", "")), str(raw.get("event_id", "")))
        if key in seen:
            continue
        seen.add(key)
        events.append(_cached_event(raw, row_by_id=row_by_id))
    return tuple(events)


def score_cached_contact_map_modes_v0(
    *,
    target_manifest: Path,
    selected_events_csv: Path,
    target_rows_csv: Path,
    coupling_dir: Path,
    output_dir: Path,
    selector_name: str = SATURATION_BOUNDARY_SELECTOR,
) -> tuple[Path, Path]:
    manifest = _load_json(target_manifest)
    selected_rows = _read_csv_rows(selected_events_csv)
    target_rows_by_id = {
        str(row.get("target_id", "")): row for row in _read_csv_rows(target_rows_csv)
    }
    output_rows: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []

    for raw_target in manifest.get("targets", []):
        if not isinstance(raw_target, Mapping):
            continue
        target_id = str(raw_target.get("target_id", "")).strip()
        if not target_id:
            continue
        try:
            benchmark_file = _resolve_path(
                raw_target.get("benchmark_file", ""),
                base_dir=BATTERY_REPO_ROOT,
            )
            external_coupling_file = _target_coupling_file(
                raw_target,
                coupling_dir=coupling_dir,
            )
            rows = load_real_coordinate_visual_rows(benchmark_file)
            row_by_id = {row.row_id: row for row in rows}
            import_result = import_external_coupling_dataset(
                rows=rows,
                external_coupling_file=external_coupling_file,
            )
            events = _events_for_target(
                selected_rows,
                target_id=target_id,
                row_by_id=row_by_id,
                selector_name=selector_name,
            )
            metric = _cached_metric(
                target_id=target_id,
                target_row=target_rows_by_id.get(target_id),
                selected_event_count=len(events),
            )
            run = TraceLoopRun(
                selector_name=selector_name,
                dataset=import_result.dataset,
                metric=metric,
                selected_events=events,
                selected_rows=(),
                constraint_count=len(import_result.dataset.constraints),
                control_kind="cached_external_real_selected_events",
            )
            exact_metrics = _exact_contact_metrics(
                rows=rows,
                dataset=import_result.dataset,
                run=run,
            )
            scaffold_metrics = _scaffold_contact_metrics(
                rows=rows,
                dataset=import_result.dataset,
                run=run,
                exact_metrics=exact_metrics,
            )
            row_out: dict[str, object] = {
                "target_id": target_id,
                "selector_name": selector_name,
                "selected_event_count": len(events),
                "cached_selected_event_count": len(events),
                "cached_source_selected_events_csv": str(selected_events_csv),
                "cached_source_target_rows_csv": str(target_rows_csv),
                "cache_claim_allowed": False,
            }
            row_out.update(exact_metrics.to_dict())
            row_out.update(scaffold_metrics.to_dict())
            row_out.update(_external_coupling_quality_summary(external_coupling_file))
            output_rows.append(row_out)
        except Exception as exc:  # noqa: BLE001
            failures.append(
                {
                    "target_id": target_id,
                    "failure_kind": exc.__class__.__name__,
                    "failure_message": str(exc),
                }
            )

    switch_audit = _apply_self_critical_quality_switch(output_rows)
    report = {
        "report_kind": CACHED_MODE_SCORE_REPORT_KIND,
        "source_report_kind": BATTERY_REPORT_KIND,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "target_manifest": str(target_manifest),
        "selected_events_csv": str(selected_events_csv),
        "target_rows_csv": str(target_rows_csv),
        "selector_name": selector_name,
        "target_count": len(output_rows),
        "failure_count": len(failures),
        "cache_claim_allowed": False,
        "classification": "cached_mode_score_inconclusive",
        "self_critical_quality_switch": switch_audit,
        "mean_phase_coverage_scaffold_precision_recall_f1": _mean_field(
            output_rows,
            "phase_coverage_scaffold_precision_recall_f1",
        ),
        "mean_region_density_top_l_scaffold_precision_recall_f1": _mean_field(
            output_rows,
            "region_density_top_l_scaffold_precision_recall_f1",
        ),
        "mean_phase_density_conflict_shell_scaffold_precision_recall_f1": _mean_field(
            output_rows,
            "phase_density_conflict_shell_scaffold_precision_recall_f1",
        ),
        "mean_phase_density_conflict_phase_confidence_scaffold_precision_recall_f1": _mean_field(
            output_rows,
            "phase_density_conflict_phase_confidence_scaffold_precision_recall_f1",
        ),
        "mean_self_critical_quality_switch_precision_recall_f1": _mean_field(
            output_rows,
            "self_critical_quality_switch_precision_recall_f1",
        ),
        "self_critical_quality_switch_contact_map_perfect_rate": _mean_field(
            output_rows,
            "self_critical_quality_switch_contact_map_perfect",
        ),
        "folding_problem_solved": False,
        "failures": failures,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "cached_contact_map_mode_score_report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    rows_path = write_csv_rows(
        output_rows,
        output_dir / "cached_contact_map_mode_score_rows.csv",
    )
    return report_path, rows_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Score contact-map completion modes from a frozen selected-events "
            "CSV without rerunning external-coupling frontier selection."
        )
    )
    parser.add_argument("--target-manifest", default=str(DEFAULT_TARGET_MANIFEST))
    parser.add_argument("--selected-events", default=str(DEFAULT_SELECTED_EVENTS))
    parser.add_argument("--target-rows", default=str(DEFAULT_TARGET_ROWS))
    parser.add_argument("--coupling-dir", default=str(DEFAULT_COUPLING_DIR))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--selector-name", default=SATURATION_BOUNDARY_SELECTOR)
    args = parser.parse_args()
    for output in score_cached_contact_map_modes_v0(
        target_manifest=Path(args.target_manifest),
        selected_events_csv=Path(args.selected_events),
        target_rows_csv=Path(args.target_rows),
        coupling_dir=Path(args.coupling_dir),
        output_dir=Path(args.output_dir),
        selector_name=args.selector_name,
    ):
        print(output)


if __name__ == "__main__":
    main()
