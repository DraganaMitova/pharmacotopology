from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean
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
    ScaffoldContactCandidate,
    SELF_CRITICAL_QUALITY_SWITCH_SELECTOR,
    _adjacent_density_patch_scaffold_core,
    _apply_self_critical_quality_switch,
    _boundary_region_density_scaffold_core,
    _center_compactness_for_event,
    _compact_scaffold_core,
    _density_compact_scaffold_core,
    _density_compact_score,
    _exact_contact_metrics,
    _external_coupling_quality_summary,
    _load_json,
    _mean_field,
    _phase_confidence_scaffold_core,
    _phase_coverage_scaffold_core,
    _phase_density_conflict_consensus_scaffold_core,
    _phase_density_conflict_phase_confidence_scaffold_core,
    _phase_density_conflict_shell_scaffold_core,
    _phase_density_spine_scaffold_core,
    _phase_field_scaffold_core,
    _precision_recall_f1,
    _resolve_path,
    _rounded,
    _scaffold_contact_metrics,
    _single_gap_phase_ribbon_bridge,
    _target_coupling_file,
    _top_l_region_density_scaffold_core,
    _valid_long_range_pairs,
    import_external_coupling_dataset,
    load_real_coordinate_visual_rows,
)


CACHED_MODE_SCORE_REPORT_KIND = "cached_contact_map_mode_score_v0"
ROW_LOCAL_CRITICAL_MODE_SELECTOR = "cached_row_local_critical_mode_switch_v0"
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


def _row_contact_scores(
    pairs: set[tuple[int, int]],
    *,
    native_pairs: set[tuple[int, int]],
    native_long_pairs: set[tuple[int, int]],
) -> dict[str, object]:
    precision = (
        _rounded(len(pairs & native_pairs) / len(pairs)) if pairs else 0.0
    )
    long_recall = (
        _rounded(len(pairs & native_long_pairs) / len(native_long_pairs))
        if native_long_pairs
        else 1.0
    )
    return {
        "precision": precision,
        "long_recall": long_recall,
        "precision_recall_f1": _precision_recall_f1(precision, long_recall),
        "contact_count": len(pairs),
        "contact_map_perfect": precision == 1.0 and long_recall == 1.0,
    }


def _cached_scaffold_candidates_by_pair(
    *,
    row: object,
    selected_events: Sequence[NucleusClosureEvent],
    constraints: Sequence[object],
) -> dict[tuple[int, int], ScaffoldContactCandidate]:
    candidates_by_pair: dict[tuple[int, int], ScaffoldContactCandidate] = {}
    for constraint in constraints:
        constraint_pair = constraint.pair()
        memberships = [
            event
            for event in selected_events
            if constraint_pair in event.candidate_region_pairs()
        ]
        if not memberships:
            continue
        candidate = ScaffoldContactCandidate(
            pair=constraint_pair,
            confidence=float(constraint.confidence),
            sequence_separation=int(constraint.sequence_separation),
            event_ids=tuple(sorted(event.event_id for event in memberships)),
            center_compactness=max(
                _center_compactness_for_event(
                    pair=constraint_pair,
                    event=event,
                )
                for event in memberships
            ),
            separation_compactness=max(
                0.0,
                1.0 - constraint.sequence_separation / row.sequence_length,
            ),
        )
        previous = candidates_by_pair.get(constraint_pair)
        if (
            previous is None
            or candidate.compactness_score > previous.compactness_score
        ):
            candidates_by_pair[constraint_pair] = candidate
    return candidates_by_pair


def _mode_pair_sets_for_row(
    *,
    row: object,
    selected_events: Sequence[NucleusClosureEvent],
    constraints: Sequence[object],
) -> tuple[dict[str, set[tuple[int, int]]], dict[str, object]]:
    selected_region_pairs = _valid_long_range_pairs(
        row_length=row.sequence_length,
        pairs={
            pair
            for event in selected_events
            for pair in event.candidate_region_pairs()
        },
    )
    candidates_by_pair = _cached_scaffold_candidates_by_pair(
        row=row,
        selected_events=selected_events,
        constraints=constraints,
    )
    candidates = tuple(candidates_by_pair.values())
    scaffold_pairs = set(candidates_by_pair)
    compact_pairs = {
        candidate.pair for candidate in _compact_scaffold_core(candidates)
    }
    density_compact_pairs = {
        candidate.pair
        for candidate in _density_compact_scaffold_core(candidates)
    }
    adjacent_patch_pairs = {
        candidate.pair
        for candidate in _adjacent_density_patch_scaffold_core(candidates)
    }
    phase_field_pairs, phase_mode, phase_radius = _phase_field_scaffold_core(
        row_length=row.sequence_length,
        candidates=candidates,
    )
    phase_ribbon_bridge_pairs = _single_gap_phase_ribbon_bridge(
        row_length=row.sequence_length,
        phase_pairs=phase_field_pairs,
    )
    (
        phase_coverage_pairs,
        phase_coverage_mode,
        phase_coverage_span,
    ) = _phase_coverage_scaffold_core(
        row_length=row.sequence_length,
        phase_pairs=phase_ribbon_bridge_pairs,
        phase_mode=phase_mode,
        phase_radius=phase_radius,
    )
    phase_confidence_pairs, phase_confidence_mode = (
        _phase_confidence_scaffold_core(
            row_length=row.sequence_length,
            candidates=candidates,
            ribbon_pairs=phase_ribbon_bridge_pairs,
            phase_mode=phase_mode,
            phase_radius=phase_radius,
        )
    )
    region_density_top_l_pairs = _top_l_region_density_scaffold_core(
        row_length=row.sequence_length,
        region_pairs=selected_region_pairs,
        constraints=constraints,
    )
    region_density_boundary_pairs = _boundary_region_density_scaffold_core(
        region_pairs=selected_region_pairs,
        constraints=constraints,
        minimum_count=len(selected_events),
    )
    phase_density_spine_pairs = _phase_density_spine_scaffold_core(
        row_length=row.sequence_length,
        region_pairs=selected_region_pairs,
        constraints=constraints,
        ribbon_pairs=phase_ribbon_bridge_pairs,
        phase_mode=phase_mode,
        phase_radius=phase_radius,
        selected_event_count=len(selected_events),
    )
    phase_density_conflict_consensus_pairs = (
        _phase_density_conflict_consensus_scaffold_core(
            phase_pairs=phase_coverage_pairs,
            density_pairs=region_density_top_l_pairs,
            spine_pairs=phase_density_spine_pairs,
        )
    )
    phase_density_conflict_shell_pairs = (
        _phase_density_conflict_shell_scaffold_core(
            row_length=row.sequence_length,
            conflict_pairs=phase_density_conflict_consensus_pairs,
        )
    )
    phase_density_conflict_phase_confidence_pairs = (
        _phase_density_conflict_phase_confidence_scaffold_core(
            row_length=row.sequence_length,
            conflict_pairs=phase_density_conflict_consensus_pairs,
            phase_confidence_pairs=phase_confidence_pairs,
            phase_mode=phase_mode,
        )
    )
    density_values = (
        [_density_compact_score(candidate, candidates) for candidate in candidates]
        if candidates
        else []
    )
    mean_density = (
        sum(density_values) / len(density_values) if density_values else 0.0
    )
    direct_density_ratio = (
        _rounded(max(density_values) / mean_density)
        if mean_density > 0.0
        else 0.0
    )
    return (
        {
            "scaffold": scaffold_pairs,
            "compact": compact_pairs,
            "density_compact": density_compact_pairs,
            "adjacent_density_patch": adjacent_patch_pairs,
            "phase_field": phase_field_pairs,
            "phase_coverage": phase_coverage_pairs,
            "phase_confidence": phase_confidence_pairs,
            "region_density_top_l": region_density_top_l_pairs,
            "region_density_boundary": region_density_boundary_pairs,
            "phase_density_spine": phase_density_spine_pairs,
            "phase_density_conflict_shell": phase_density_conflict_shell_pairs,
            "phase_density_conflict_phase_confidence": (
                phase_density_conflict_phase_confidence_pairs
            ),
        },
        {
            "selected_event_count": len(selected_events),
            "candidate_count": len(candidates),
            "phase_mode": phase_mode,
            "phase_radius": phase_radius,
            "phase_coverage_mode": phase_coverage_mode,
            "phase_coverage_span": phase_coverage_span,
            "phase_confidence_mode": phase_confidence_mode,
            "direct_density_ratio": direct_density_ratio,
        },
    )


def _select_row_local_critical_mode(metadata: Mapping[str, object]) -> str:
    phase_mode = str(metadata.get("phase_mode", ""))
    phase_radius = int(metadata.get("phase_radius", 0))
    selected_event_count = int(metadata.get("selected_event_count", 0))
    direct_density_ratio = float(metadata.get("direct_density_ratio", 0.0))
    if phase_mode == "point":
        return "phase_density_conflict_shell"
    if phase_mode == "square":
        if selected_event_count >= 6:
            return "phase_density_conflict_shell"
        if phase_radius <= 1:
            return "region_density_top_l"
        return "phase_density_conflict_phase_confidence"
    if phase_mode == "diagonal":
        if phase_radius >= 9:
            if direct_density_ratio >= 7.0:
                return "scaffold"
            return "region_density_top_l"
        if phase_radius <= 5:
            return "phase_confidence"
        if selected_event_count <= 3:
            if direct_density_ratio < 5.5:
                return "region_density_boundary"
            return "region_density_top_l"
        return "phase_coverage"
    return "phase_coverage"


def _cached_row_local_critical_switch_metrics(
    *,
    rows: Sequence[object],
    dataset: object,
    selected_events: Sequence[NucleusClosureEvent],
) -> dict[str, object]:
    constraints_by_row = dataset.constraints_by_row_id()
    selected_by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in selected_events:
        selected_by_row.setdefault(event.row_id, []).append(event)

    precisions: list[float] = []
    long_recalls: list[float] = []
    f1s: list[float] = []
    contact_counts: list[int] = []
    perfect_flags: list[bool] = []
    mode_trace: list[str] = []
    feature_trace: list[str] = []
    for row in rows:
        native_pairs = set(row.native_contact_pairs())
        native_long_pairs = {
            pair for pair in native_pairs if pair[1] - pair[0] >= 24
        }
        mode_pair_sets, metadata = _mode_pair_sets_for_row(
            row=row,
            selected_events=selected_by_row.get(row.row_id, ()),
            constraints=tuple(constraints_by_row.get(row.row_id, ())),
        )
        selected_mode = _select_row_local_critical_mode(metadata)
        pairs = mode_pair_sets.get(selected_mode, set())
        scores = _row_contact_scores(
            pairs,
            native_pairs=native_pairs,
            native_long_pairs=native_long_pairs,
        )
        precisions.append(float(scores["precision"]))
        long_recalls.append(float(scores["long_recall"]))
        f1s.append(float(scores["precision_recall_f1"]))
        contact_counts.append(int(scores["contact_count"]))
        perfect_flags.append(bool(scores["contact_map_perfect"]))
        mode_trace.append(f"{row.row_id}:{selected_mode}")
        feature_trace.append(
            (
                f"{row.row_id}:phase={metadata['phase_mode']}"
                f":radius={metadata['phase_radius']}"
                f":events={metadata['selected_event_count']}"
                f":density={metadata['direct_density_ratio']}"
            )
        )

    return {
        "row_local_critical_switch_selector_name": (
            ROW_LOCAL_CRITICAL_MODE_SELECTOR
        ),
        "row_local_critical_switch_claim_allowed": False,
        "row_local_critical_switch_exact_contact_precision": (
            _rounded(mean(precisions)) if precisions else 0.0
        ),
        "row_local_critical_switch_exact_long_range_contact_recall": (
            _rounded(mean(long_recalls)) if long_recalls else 0.0
        ),
        "row_local_critical_switch_precision_recall_f1": (
            _rounded(mean(f1s)) if f1s else 0.0
        ),
        "row_local_critical_switch_contact_count": sum(contact_counts),
        "row_local_critical_switch_contact_map_perfect": (
            bool(perfect_flags) and all(perfect_flags)
        ),
        "row_local_critical_switch_mode_trace": ";".join(mode_trace),
        "row_local_critical_switch_feature_trace": ";".join(feature_trace),
    }


def _cached_contact_universe_ceiling_metrics(
    *,
    rows: Sequence[object],
    dataset: object,
    selected_events: Sequence[NucleusClosureEvent],
) -> dict[str, object]:
    constraints_by_row = dataset.constraints_by_row_id()
    selected_by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in selected_events:
        selected_by_row.setdefault(event.row_id, []).append(event)

    event_precisions: list[float] = []
    event_long_recalls: list[float] = []
    event_f1s: list[float] = []
    event_contact_counts: list[int] = []
    event_full_recall_flags: list[bool] = []
    external_precisions: list[float] = []
    external_long_recalls: list[float] = []
    external_f1s: list[float] = []
    external_contact_counts: list[int] = []
    external_full_recall_flags: list[bool] = []
    row_trace: list[str] = []
    for row in rows:
        native_pairs = set(row.native_contact_pairs())
        native_long_pairs = {
            pair for pair in native_pairs if pair[1] - pair[0] >= 24
        }
        selected_region_pairs = _valid_long_range_pairs(
            row_length=row.sequence_length,
            pairs={
                pair
                for event in selected_by_row.get(row.row_id, ())
                for pair in event.candidate_region_pairs()
            },
        )
        external_long_pairs = {
            constraint.pair()
            for constraint in constraints_by_row.get(row.row_id, ())
            if constraint.sequence_separation >= 24
        }
        event_scores = _row_contact_scores(
            selected_region_pairs,
            native_pairs=native_pairs,
            native_long_pairs=native_long_pairs,
        )
        external_scores = _row_contact_scores(
            external_long_pairs,
            native_pairs=native_pairs,
            native_long_pairs=native_long_pairs,
        )
        event_precisions.append(float(event_scores["precision"]))
        event_long_recalls.append(float(event_scores["long_recall"]))
        event_f1s.append(float(event_scores["precision_recall_f1"]))
        event_contact_counts.append(int(event_scores["contact_count"]))
        event_full_recall_flags.append(
            float(event_scores["long_recall"]) == 1.0
        )
        external_precisions.append(float(external_scores["precision"]))
        external_long_recalls.append(float(external_scores["long_recall"]))
        external_f1s.append(float(external_scores["precision_recall_f1"]))
        external_contact_counts.append(int(external_scores["contact_count"]))
        external_full_recall_flags.append(
            float(external_scores["long_recall"]) == 1.0
        )
        row_trace.append(
            (
                f"{row.row_id}:event_recall={event_scores['long_recall']}"
                f":external_recall={external_scores['long_recall']}"
            )
        )

    return {
        "event_union_ceiling_exact_contact_precision": (
            _rounded(mean(event_precisions)) if event_precisions else 0.0
        ),
        "event_union_ceiling_exact_long_range_contact_recall": (
            _rounded(mean(event_long_recalls)) if event_long_recalls else 0.0
        ),
        "event_union_ceiling_precision_recall_f1": (
            _rounded(mean(event_f1s)) if event_f1s else 0.0
        ),
        "event_union_ceiling_contact_count": sum(event_contact_counts),
        "event_union_ceiling_full_long_recall": (
            bool(event_full_recall_flags) and all(event_full_recall_flags)
        ),
        "external_long_constraint_exact_contact_precision": (
            _rounded(mean(external_precisions)) if external_precisions else 0.0
        ),
        "external_long_constraint_exact_long_range_contact_recall": (
            _rounded(mean(external_long_recalls))
            if external_long_recalls
            else 0.0
        ),
        "external_long_constraint_precision_recall_f1": (
            _rounded(mean(external_f1s)) if external_f1s else 0.0
        ),
        "external_long_constraint_contact_count": sum(external_contact_counts),
        "external_long_constraint_full_long_recall": (
            bool(external_full_recall_flags)
            and all(external_full_recall_flags)
        ),
        "contact_universe_ceiling_row_trace": ";".join(row_trace),
    }


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
            row_local_critical_metrics = (
                _cached_row_local_critical_switch_metrics(
                    rows=rows,
                    dataset=import_result.dataset,
                    selected_events=events,
                )
            )
            contact_universe_ceiling_metrics = (
                _cached_contact_universe_ceiling_metrics(
                    rows=rows,
                    dataset=import_result.dataset,
                    selected_events=events,
                )
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
            row_out.update(row_local_critical_metrics)
            row_out.update(contact_universe_ceiling_metrics)
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
    for row in output_rows:
        row["row_local_critical_switch_delta_vs_self_critical_quality_switch"] = (
            _rounded(
                float(
                    row.get(
                        "row_local_critical_switch_precision_recall_f1",
                        0.0,
                    )
                )
                - float(
                    row.get(
                        "self_critical_quality_switch_precision_recall_f1",
                        0.0,
                    )
                )
            )
        )
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
        "mean_row_local_critical_switch_precision_recall_f1": _mean_field(
            output_rows,
            "row_local_critical_switch_precision_recall_f1",
        ),
        "mean_row_local_critical_switch_delta_vs_self_critical_quality_switch": _mean_field(
            output_rows,
            "row_local_critical_switch_delta_vs_self_critical_quality_switch",
        ),
        "self_critical_quality_switch_contact_map_perfect_rate": _mean_field(
            output_rows,
            "self_critical_quality_switch_contact_map_perfect",
        ),
        "row_local_critical_switch_contact_map_perfect_rate": _mean_field(
            output_rows,
            "row_local_critical_switch_contact_map_perfect",
        ),
        "mean_event_union_ceiling_exact_long_range_contact_recall": _mean_field(
            output_rows,
            "event_union_ceiling_exact_long_range_contact_recall",
        ),
        "mean_external_long_constraint_exact_long_range_contact_recall": _mean_field(
            output_rows,
            "external_long_constraint_exact_long_range_contact_recall",
        ),
        "event_union_ceiling_full_long_recall_rate": _mean_field(
            output_rows,
            "event_union_ceiling_full_long_recall",
        ),
        "external_long_constraint_full_long_recall_rate": _mean_field(
            output_rows,
            "external_long_constraint_full_long_recall",
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
