from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Mapping, Optional, Sequence

from pharmacotopology.artifact_io import write_csv_rows
from pharmacotopology.folding_coupling_negative_controls import (
    COUPLING_ADVERSARIAL_CALIBRATED_CONTROL_KIND,
    CouplingControlDataset,
    generate_adversarial_calibrated_external_coupling_controls,
    generate_external_coupling_negative_controls,
)
from pharmacotopology.folding_coupling_nucleus_selector import (
    CouplingNucleusContext,
    CouplingSelectorMetric,
    TRACE_LOOP_RANK_CONSISTENT_CLUSTER_GATE_MIN,
    TRACE_LOOP_RANK_CONSISTENT_RECOVERY_BLOCKED_FUTURE_MAX,
    TRACE_LOOP_RANK_CONSISTENT_RECOVERY_DECOY_MARGIN_MIN,
    TRACE_LOOP_RANK_CONSISTENT_RECOVERY_DIRECT_SUPPORT_MIN,
    TRACE_LOOP_RANK_CONSISTENT_RECOVERY_FUTURE_PRESERVATION_MIN,
    build_coupling_nucleus_context,
    coupling_claim_mode_validation_failures,
    coupling_nucleus_score,
    select_coupling_events,
    selected_event_rows,
    selector_metrics,
)
from pharmacotopology.folding_evolutionary_constraints import (
    CouplingDataset,
    compatible_future_event,
    load_coupling_dataset,
)
from pharmacotopology.folding_coupling_decoy_falsification import (
    coupling_decoy_comparisons,
)
from pharmacotopology.folding_external_coupling_importer import (
    ExternalCouplingImportResult,
    import_external_coupling_dataset,
)
from pharmacotopology.folding_nucleus_decoy_falsification import (
    decoy_distance,
    matched_decoys_for_selected_events,
)
from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent
from pharmacotopology.folding_physical_selection import (
    ActivePhysicalContext,
    build_active_physical_context,
)
from pharmacotopology.folding_external_coupling_sources import (
    EXTERNAL_COUPLING_TRACE_LOOP_CERTIFICATE_KIND,
    EXTERNAL_COUPLING_TRACE_LOOP_REPORT_KIND,
    EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID,
    SERIOUS_EXTERNAL_COUPLING_POLICY,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


ROOT_OUTPUT_NAMES = (
    "external_coupling_trace_loop_report.json",
    "external_coupling_trace_loop_certificate.json",
    "external_coupling_trace_loop_selectors.csv",
    "external_coupling_trace_loop_selected_events.csv",
    "external_coupling_trace_loop_frontier.csv",
    "external_coupling_trace_loop_controls.csv",
    "external_coupling_trace_loop_row_status.csv",
    "external_coupling_trace_loop_dashboard.html",
)

MATCHED_CONTROL_NAMES = (
    "external_shuffled_same_row_same_separation",
    "external_confidence_permuted",
    "external_cross_row_swapped",
    "external_random_long_range_same_count",
    "external_low_confidence_tail",
)
NEGATIVE_CONTROL_KIND = "external_coupling_negative_controls_v0"
ADVERSARIAL_CALIBRATED_CONTROL_KIND = (
    COUPLING_ADVERSARIAL_CALIBRATED_CONTROL_KIND
)
ADVERSARIAL_ENRICHMENT_MIN_SELECTED_EVENTS = 4
SCORE_MARGIN_EXPANSION_SCORE_MIN = 0.44
SCORE_MARGIN_EXPANSION_DECOY_MARGIN_MIN = 0.15
SCORE_MARGIN_EXPANSION_CLUSTER_MIN = 0.46
SCORE_MARGIN_EXPANSION_DIRECT_SUPPORT_MIN = 0.10
SCORE_MARGIN_EXPANSION_FUTURE_PRESERVATION_MIN = 0.18
SCORE_MARGIN_EXPANSION_BLOCKED_FUTURE_MAX = 0.08
MACRO_SCALE_SEGMENT_LENGTH = 20
MACRO_SCALE_SEGMENT_STRIDE = 4
MULTISCALE_FUTURE_PRESERVED_SEGMENT_STRIDE = 4
MULTISCALE_FUTURE_PRESERVED_MAX_EVENTS_PER_ROW = 14
MULTISCALE_FUTURE_PRESERVED_CONFIGS = (
    (18, 0.24),
    (20, 0.36),
    (24, 0.36),
    (32, 0.34),
    (36, 0.36),
    (40, 0.40),
)


@dataclass(frozen=True)
class TraceLoopRun:
    selector_name: str
    dataset: CouplingDataset
    metric: CouplingSelectorMetric
    selected_events: tuple[NucleusClosureEvent, ...]
    selected_rows: tuple[dict[str, object], ...]
    constraint_count: int
    control_kind: str


@dataclass(frozen=True)
class MultiScaleSelectedEvent:
    segment_length: int
    future_preservation_min: float
    context: CouplingNucleusContext
    event: NucleusClosureEvent


def _rounded(value: float) -> float:
    return round(value, 6)


def _mean(values: Sequence[float]) -> float:
    return _rounded(mean(values)) if values else 0.0


def _max_metric(
    runs: Sequence[TraceLoopRun],
    field_name: str,
) -> float:
    values = [float(getattr(run.metric, field_name)) for run in runs]
    return _rounded(max(values)) if values else 0.0


def _max_selected_metric(
    runs: Sequence[TraceLoopRun],
    field_name: str,
) -> float:
    values = [
        float(getattr(run.metric, field_name))
        for run in runs
        if run.metric.selected_event_count > 0
    ]
    return _rounded(max(values)) if values else 0.0


def _run_trace_loop_selector(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
    selector_name: str,
    selection_mode: str = "coupling_trace_loop",
    control_kind: str = "external_real",
) -> TraceLoopRun:
    context = build_coupling_nucleus_context(
        rows=rows,
        coupling_dataset=dataset,
    )
    return _run_trace_loop_selector_from_context(
        context=context,
        dataset=dataset,
        selector_name=selector_name,
        selection_mode=selection_mode,
        control_kind=control_kind,
    )


def _run_trace_loop_selector_from_context(
    *,
    context: CouplingNucleusContext,
    dataset: CouplingDataset,
    selector_name: str,
    selection_mode: str = "coupling_trace_loop",
    control_kind: str = "external_real",
) -> TraceLoopRun:
    selected = select_coupling_events(context, selector_name=selection_mode)
    metric = selector_metrics(
        context,
        selector_name=selector_name,
        selected_events=selected,
    )
    rows_out = selected_event_rows(context, {selector_name: selected})
    return TraceLoopRun(
        selector_name=selector_name,
        dataset=dataset,
        metric=metric,
        selected_events=tuple(selected),
        selected_rows=tuple(rows_out),
        constraint_count=len(dataset.constraints),
        control_kind=control_kind,
    )


def _build_macro_scale_coupling_context(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
) -> CouplingNucleusContext:
    macro_physical_context = build_active_physical_context(
        rows,
        segment_length=MACRO_SCALE_SEGMENT_LENGTH,
        segment_stride=MACRO_SCALE_SEGMENT_STRIDE,
    )
    return build_coupling_nucleus_context(
        rows=rows,
        coupling_dataset=dataset,
        physical_context=macro_physical_context,
    )


def _build_multiscale_physical_contexts(
    rows: Sequence[RealCoordinateVisualRow],
) -> dict[int, ActivePhysicalContext]:
    return {
        segment_length: build_active_physical_context(
            rows,
            segment_length=segment_length,
            segment_stride=MULTISCALE_FUTURE_PRESERVED_SEGMENT_STRIDE,
        )
        for segment_length, _ in MULTISCALE_FUTURE_PRESERVED_CONFIGS
    }


def _multiscale_metric(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
    selector_name: str,
    selected_items: Sequence[MultiScaleSelectedEvent],
) -> CouplingSelectorMetric:
    selected_by_row: dict[str, list[MultiScaleSelectedEvent]] = {}
    for item in selected_items:
        selected_by_row.setdefault(item.event.row_id, []).append(item)
    constraints_by_row = dataset.constraints_by_row_id()
    false_rates: list[float] = []
    precisions: list[float] = []
    long_recalls: list[float] = []
    constraint_recalls: list[float] = []
    selectivity_scores: list[float] = []
    nucleus_scores: list[float] = []
    for row in rows:
        row_items = selected_by_row.get(row.row_id, [])
        row_events = [item.event for item in row_items]
        native_pairs = set(row.native_contact_pairs())
        native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
        constraint_pairs = {
            constraint.pair()
            for constraint in constraints_by_row.get(row.row_id, ())
        }
        region: set[tuple[int, int]] = set()
        for event in row_events:
            region.update(event.candidate_region_pairs())
        native_hit_count = sum(
            event.native_contact_count_after_scoring for event in row_events
        )
        possible_region_pair_count = sum(
            len(event.candidate_region_pairs()) for event in row_events
        )
        false_count = sum(
            1
            for event in row_events
            if event.native_contact_count_after_scoring == 0
        )
        false_rates.append(false_count / len(row_events) if row_events else 0.0)
        precisions.append(
            native_hit_count / possible_region_pair_count
            if possible_region_pair_count
            else 0.0
        )
        long_recalls.append(
            len(region & native_long) / len(native_long) if native_long else 1.0
        )
        constraint_recalls.append(
            len(region & constraint_pairs) / len(constraint_pairs)
            if constraint_pairs
            else 0.0
        )
        for item in row_items:
            assessment = item.context.assessment_by_event_id[item.event.event_id]
            selectivity_scores.append(assessment.coupling_selectivity_score)
            nucleus_scores.append(coupling_nucleus_score(item.event, item.context))
    false_rate = _rounded(mean(false_rates) if false_rates else 0.0)
    precision = _rounded(mean(precisions) if precisions else 0.0)
    long_recall = _rounded(mean(long_recalls) if long_recalls else 0.0)
    selected_selectivity_mean = _rounded(
        mean(selectivity_scores) if selectivity_scores else 0.0
    )
    selected_nucleus_mean = _rounded(
        mean(nucleus_scores) if nucleus_scores else 0.0
    )
    return CouplingSelectorMetric(
        selector_name=selector_name,
        selected_event_count=len(selected_items),
        false_nucleus_rate=false_rate,
        contact_cluster_precision=precision,
        long_range_contact_recall=long_recall,
        coupling_constraint_recall=_rounded(
            mean(constraint_recalls) if constraint_recalls else 0.0
        ),
        real_vs_decoy_coupling_enrichment_ratio=0.0,
        real_beats_decoy_coupling_score_rate=0.0,
        mean_selected_coupling_selectivity_score=selected_selectivity_mean,
        mean_decoy_coupling_selectivity_score=0.0,
        mean_coupling_decoy_selectivity_margin=selected_selectivity_mean,
        mean_coupling_nucleus_score=selected_nucleus_mean,
        mean_decoy_coupling_nucleus_score=0.0,
        mean_coupling_nucleus_decoy_margin=selected_nucleus_mean,
        real_vs_decoy_coupling_nucleus_enrichment_ratio=0.0,
        real_beats_decoy_coupling_nucleus_score_rate=0.0,
        survives_targets=False,
        coordinate_truth_used_to_build_constraints=dataset.coordinate_truth_tainted,
        native_truth_used_before_coupling_selection=(
            dataset.native_truth_used_before_coupling_selection
        ),
        raw_sequence_exposed=False,
    )


def _multiscale_selected_rows(
    *,
    selector_name: str,
    selected_items: Sequence[MultiScaleSelectedEvent],
) -> tuple[dict[str, object], ...]:
    rows: list[dict[str, object]] = []
    for rank, item in enumerate(selected_items, start=1):
        event = item.event
        context = item.context
        assessment = context.assessment_by_event_id[event.event_id]
        state = context.physical_context.state_by_event_id[event.event_id]
        rows.append(
            {
                "selector_name": selector_name,
                "rank": rank,
                "row_id": event.row_id,
                "source_accession": event.source_accession,
                "event_id": event.event_id,
                "segment_a_start": event.segment_a_start,
                "segment_a_end": event.segment_a_end,
                "segment_b_start": event.segment_b_start,
                "segment_b_end": event.segment_b_end,
                "normalized_span": event.normalized_span,
                "coupling_nucleus_score": coupling_nucleus_score(event, context),
                "coupling_selectivity_score": (
                    assessment.coupling_selectivity_score
                ),
                "direct_support_score": assessment.direct_support_score,
                "future_preservation_score": (
                    assessment.future_preservation_score
                ),
                "blocked_future_pressure": assessment.blocked_future_pressure,
                "coupling_decoy_margin": (
                    context.coupling_decoy_margin_by_event_id[event.event_id]
                ),
                "physical_state_score": state.physical_state_score,
                "burial_gain": state.burial_gain,
                "unsatisfied_polar_penalty": state.unsatisfied_polar_penalty,
                "native_contact_count_after_scoring": (
                    event.native_contact_count_after_scoring
                ),
                "coordinate_truth_used_to_build_constraints": (
                    context.coupling_dataset.coordinate_truth_tainted
                ),
                "native_truth_used_before_coupling_selection": (
                    context.coupling_dataset.native_truth_used_before_coupling_selection
                ),
                "raw_sequence_exposed": False,
                "multiscale_segment_length": item.segment_length,
                "multiscale_future_preservation_min": (
                    item.future_preservation_min
                ),
            }
        )
    return tuple(rows)


def _run_multiscale_future_preserved_selector(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
    selector_name: str,
    control_kind: str,
    physical_contexts: Mapping[int, ActivePhysicalContext],
) -> TraceLoopRun:
    candidates_by_row: dict[str, list[MultiScaleSelectedEvent]] = {
        row.row_id: [] for row in rows
    }
    for segment_length, future_preservation_min in (
        MULTISCALE_FUTURE_PRESERVED_CONFIGS
    ):
        context = build_coupling_nucleus_context(
            rows=rows,
            coupling_dataset=dataset,
            physical_context=physical_contexts[segment_length],
        )
        candidates = select_coupling_events(
            context,
            selector_name="coupling_trace_loop_boundary_field_replacement_probe",
        )
        for event in candidates:
            assessment = context.assessment_by_event_id[event.event_id]
            if assessment.future_preservation_score < future_preservation_min:
                continue
            candidates_by_row[event.row_id].append(
                MultiScaleSelectedEvent(
                    segment_length=segment_length,
                    future_preservation_min=future_preservation_min,
                    context=context,
                    event=event,
                )
            )
    selected_items: list[MultiScaleSelectedEvent] = []
    for row in rows:
        row_candidates = sorted(
            candidates_by_row[row.row_id],
            key=lambda item: (
                item.context.assessment_by_event_id[
                    item.event.event_id
                ].future_preservation_score,
                item.context.assessment_by_event_id[
                    item.event.event_id
                ].direct_support_score,
                item.event.closure_event_stability,
                -item.segment_length,
            ),
            reverse=True,
        )
        selected_items.extend(
            row_candidates[:MULTISCALE_FUTURE_PRESERVED_MAX_EVENTS_PER_ROW]
        )
    metric = _multiscale_metric(
        rows=rows,
        dataset=dataset,
        selector_name=selector_name,
        selected_items=selected_items,
    )
    return TraceLoopRun(
        selector_name=selector_name,
        dataset=dataset,
        metric=metric,
        selected_events=tuple(item.event for item in selected_items),
        selected_rows=_multiscale_selected_rows(
            selector_name=selector_name,
            selected_items=selected_items,
        ),
        constraint_count=len(dataset.constraints),
        control_kind=control_kind,
    )


def _direct_constraint_stats(
    *,
    context: CouplingNucleusContext,
    event: NucleusClosureEvent,
) -> dict[str, object]:
    constraints = context.coupling_dataset.constraints_by_row_id().get(
        event.row_id,
        (),
    )
    region_pairs = set(event.candidate_region_pairs())
    direct = tuple(
        constraint for constraint in constraints if constraint.pair() in region_pairs
    )
    if not direct:
        return {
            "direct_constraint_count": 0,
            "direct_constraint_confidence_sum": 0.0,
            "direct_top_10pct_rank_count": 0,
            "direct_max_confidence": 0.0,
            "direct_mean_rank_fraction": 0.0,
        }
    return {
        "direct_constraint_count": len(direct),
        "direct_constraint_confidence_sum": _rounded(
            sum(constraint.confidence for constraint in direct)
        ),
        "direct_top_10pct_rank_count": sum(
            1 for constraint in direct if constraint.rank_fraction <= 0.10
        ),
        "direct_max_confidence": _rounded(
            max(constraint.confidence for constraint in direct)
        ),
        "direct_mean_rank_fraction": _rounded(
            mean(constraint.rank_fraction for constraint in direct)
        ),
    }


def _frontier_exclusion_reasons(
    *,
    context: CouplingNucleusContext,
    event: NucleusClosureEvent,
) -> tuple[str, ...]:
    assessment = context.assessment_by_event_id[event.event_id]
    coupling_decoy_margin = context.coupling_decoy_margin_by_event_id[event.event_id]
    reasons: list[str] = []
    if event.contact_cluster_gain < TRACE_LOOP_RANK_CONSISTENT_CLUSTER_GATE_MIN:
        reasons.append("below_hard_cluster_gate")
    if (
        assessment.direct_support_score
        < TRACE_LOOP_RANK_CONSISTENT_RECOVERY_DIRECT_SUPPORT_MIN
    ):
        reasons.append("below_recovery_direct_support")
    if (
        assessment.blocked_future_pressure
        > TRACE_LOOP_RANK_CONSISTENT_RECOVERY_BLOCKED_FUTURE_MAX
    ):
        reasons.append("above_recovery_blocked_future")
    if (
        assessment.future_preservation_score
        < TRACE_LOOP_RANK_CONSISTENT_RECOVERY_FUTURE_PRESERVATION_MIN
        and coupling_decoy_margin
        < TRACE_LOOP_RANK_CONSISTENT_RECOVERY_DECOY_MARGIN_MIN
    ):
        reasons.append("below_recovery_future_or_decoy_margin")
    return tuple(reasons)


def _matched_selector_score_decoy(
    *,
    context: CouplingNucleusContext,
    event: NucleusClosureEvent,
) -> tuple[NucleusClosureEvent, float]:
    row_candidates = [
        candidate
        for candidate in context.competitive_events
        if candidate.row_id == event.row_id and candidate.event_id != event.event_id
    ]
    decoy = min(
        row_candidates,
        key=lambda candidate: decoy_distance(event, candidate),
        default=event,
    )
    margin = _rounded(
        coupling_nucleus_score(event, context)
        - coupling_nucleus_score(decoy, context)
    )
    return decoy, margin


def _score_margin_expansion_reasons(
    *,
    context: CouplingNucleusContext,
    event: NucleusClosureEvent,
    selector_score_margin: float,
) -> tuple[str, ...]:
    assessment = context.assessment_by_event_id[event.event_id]
    score = coupling_nucleus_score(event, context)
    reasons: list[str] = []
    if score < SCORE_MARGIN_EXPANSION_SCORE_MIN:
        reasons.append("below_selector_score")
    if selector_score_margin < SCORE_MARGIN_EXPANSION_DECOY_MARGIN_MIN:
        reasons.append("below_selector_score_decoy_margin")
    if event.contact_cluster_gain < SCORE_MARGIN_EXPANSION_CLUSTER_MIN:
        reasons.append("below_expansion_cluster")
    if assessment.direct_support_score < SCORE_MARGIN_EXPANSION_DIRECT_SUPPORT_MIN:
        reasons.append("below_expansion_direct_support")
    if (
        assessment.future_preservation_score
        < SCORE_MARGIN_EXPANSION_FUTURE_PRESERVATION_MIN
    ):
        reasons.append("below_expansion_future_preservation")
    if assessment.blocked_future_pressure > SCORE_MARGIN_EXPANSION_BLOCKED_FUTURE_MAX:
        reasons.append("above_expansion_blocked_future")
    return tuple(reasons)


def persistent_recall_frontier_rows(
    *,
    context: CouplingNucleusContext,
    persistent_run: TraceLoopRun,
) -> list[dict[str, object]]:
    selected_ids = {event.event_id for event in persistent_run.selected_events}
    selected_by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in persistent_run.selected_events:
        selected_by_row.setdefault(event.row_id, []).append(event)

    trace_events = select_coupling_events(
        context,
        selector_name="coupling_trace_loop",
    )
    frontier_events = tuple(
        event
        for event in trace_events
        if event.event_id not in selected_ids
        and all(
            compatible_future_event(selected_event, event)
            for selected_event in selected_by_row.get(event.row_id, ())
        )
    )
    rows: list[dict[str, object]] = []
    for rank, event in enumerate(
        sorted(
            frontier_events,
            key=lambda item: (
                -item.native_long_range_contact_count_after_scoring,
                -item.native_contact_count_after_scoring,
                -coupling_nucleus_score(item, context),
                item.row_id,
                item.segment_a_start,
                item.segment_b_start,
                item.event_id,
            ),
        ),
        start=1,
    ):
        assessment = context.assessment_by_event_id[event.event_id]
        decoy, selector_score_margin = _matched_selector_score_decoy(
            context=context,
            event=event,
        )
        decoy_assessment = context.assessment_by_event_id[decoy.event_id]
        reasons = _score_margin_expansion_reasons(
            context=context,
            event=event,
            selector_score_margin=selector_score_margin,
        )
        gate_passed = not reasons
        native_positive = event.native_contact_count_after_scoring > 0
        if not native_positive and not gate_passed:
            continue
        rows.append(
            {
                "frontier_kind": "persistent_score_margin_recall_frontier_v0",
                "source_selector": "coupling_trace_loop",
                "target_selector": persistent_run.selector_name,
                "row_id": event.row_id,
                "source_accession": event.source_accession,
                "event_id": event.event_id,
                "segment_a_start": event.segment_a_start,
                "segment_a_end": event.segment_a_end,
                "segment_b_start": event.segment_b_start,
                "segment_b_end": event.segment_b_end,
                "native_contact_count_after_scoring": (
                    event.native_contact_count_after_scoring
                ),
                "native_long_range_contact_count_after_scoring": (
                    event.native_long_range_contact_count_after_scoring
                ),
                "contact_cluster_gain": event.contact_cluster_gain,
                "coupling_selectivity_score": assessment.coupling_selectivity_score,
                "direct_support_score": assessment.direct_support_score,
                "future_preservation_score": assessment.future_preservation_score,
                "blocked_future_pressure": assessment.blocked_future_pressure,
                "coupling_decoy_margin": (
                    context.coupling_decoy_margin_by_event_id[event.event_id]
                ),
                "coupling_nucleus_score": coupling_nucleus_score(event, context),
                "selector_score_decoy_margin": selector_score_margin,
                "exclusion_reasons": ";".join(reasons),
                "matched_decoy_event_id": decoy.event_id,
                "matched_decoy_native_positive_after_scoring": (
                    decoy.native_contact_count_after_scoring > 0
                ),
                "matched_decoy_contact_cluster_gain": decoy.contact_cluster_gain,
                "matched_decoy_coupling_selectivity_score": (
                    decoy_assessment.coupling_selectivity_score
                ),
                "matched_decoy_direct_support_score": (
                    decoy_assessment.direct_support_score
                ),
                "matched_decoy_future_preservation_score": (
                    decoy_assessment.future_preservation_score
                ),
                "matched_decoy_blocked_future_pressure": (
                    decoy_assessment.blocked_future_pressure
                ),
                "matched_decoy_coupling_nucleus_score": coupling_nucleus_score(
                    decoy,
                    context,
                ),
                "real_beats_decoy_by_coupling_score": (
                    assessment.coupling_selectivity_score
                    > decoy_assessment.coupling_selectivity_score
                ),
                "real_beats_decoy_by_coupling_nucleus_score": (
                    selector_score_margin > 0.0
                ),
                "score_margin_expansion_gate_passed": gate_passed,
                "recall_frontier_rank": rank,
                "native_truth_used_before_selection": False,
                "native_label_attached_after_event_generation": True,
                "diagnostic_claim_allowed": False,
            }
        )
        rows[-1].update(_direct_constraint_stats(context=context, event=event))
    return rows


def rank_consistent_frontier_rows(
    *,
    context: CouplingNucleusContext,
    cluster_gated_run: TraceLoopRun,
    rank_consistent_run: TraceLoopRun,
) -> list[dict[str, object]]:
    selected_ids = {event.event_id for event in rank_consistent_run.selected_events}
    frontier_events = tuple(
        event
        for event in cluster_gated_run.selected_events
        if event.event_id not in selected_ids
        and event.native_contact_count_after_scoring > 0
    )
    matches = matched_decoys_for_selected_events(
        selected_events=frontier_events,
        candidate_events=context.competitive_events,
    )
    comparisons = coupling_decoy_comparisons(
        matches=matches,
        assessments=context.assessments,
    )
    match_by_event_id = {match.real_event_id: match for match in matches}
    comparison_by_event_id = {
        comparison.real_event_id: comparison for comparison in comparisons
    }

    rows: list[dict[str, object]] = []
    for event in sorted(
        frontier_events,
        key=lambda item: (
            item.row_id,
            item.segment_a_start,
            item.segment_b_start,
            item.event_id,
        ),
    ):
        assessment = context.assessment_by_event_id[event.event_id]
        coupling_decoy_margin = context.coupling_decoy_margin_by_event_id[
            event.event_id
        ]
        match = match_by_event_id[event.event_id]
        comparison = comparison_by_event_id[event.event_id]
        decoy_event = context.event_by_id[match.decoy_event_id]
        decoy_assessment = context.assessment_by_event_id[match.decoy_event_id]
        selector_score_decoy_margin = _rounded(
            coupling_nucleus_score(event, context)
            - coupling_nucleus_score(decoy_event, context)
        )
        row = {
            "frontier_kind": "rank_consistent_native_after_scoring_frontier_v0",
            "source_selector": cluster_gated_run.selector_name,
            "target_selector": rank_consistent_run.selector_name,
            "row_id": event.row_id,
            "source_accession": event.source_accession,
            "event_id": event.event_id,
            "segment_a_start": event.segment_a_start,
            "segment_a_end": event.segment_a_end,
            "segment_b_start": event.segment_b_start,
            "segment_b_end": event.segment_b_end,
            "native_contact_count_after_scoring": (
                event.native_contact_count_after_scoring
            ),
            "native_long_range_contact_count_after_scoring": (
                event.native_long_range_contact_count_after_scoring
            ),
            "contact_cluster_gain": event.contact_cluster_gain,
            "coupling_selectivity_score": assessment.coupling_selectivity_score,
            "direct_support_score": assessment.direct_support_score,
            "future_preservation_score": assessment.future_preservation_score,
            "blocked_future_pressure": assessment.blocked_future_pressure,
            "coupling_decoy_margin": coupling_decoy_margin,
            "coupling_nucleus_score": coupling_nucleus_score(event, context),
            "selector_score_decoy_margin": selector_score_decoy_margin,
            "exclusion_reasons": ";".join(
                _frontier_exclusion_reasons(context=context, event=event)
            ),
            "matched_decoy_event_id": match.decoy_event_id,
            "matched_decoy_native_positive_after_scoring": (
                match.decoy_native_positive_after_scoring
            ),
            "matched_decoy_contact_cluster_gain": match.decoy_contact_cluster_gain,
            "matched_decoy_coupling_selectivity_score": (
                comparison.decoy_coupling_selectivity_score
            ),
            "matched_decoy_direct_support_score": (
                decoy_assessment.direct_support_score
            ),
            "matched_decoy_future_preservation_score": (
                decoy_assessment.future_preservation_score
            ),
            "matched_decoy_blocked_future_pressure": (
                decoy_assessment.blocked_future_pressure
            ),
            "matched_decoy_coupling_nucleus_score": coupling_nucleus_score(
                decoy_event,
                context,
            ),
            "real_beats_decoy_by_coupling_score": (
                comparison.real_beats_decoy_by_coupling_score
            ),
            "real_beats_decoy_by_coupling_nucleus_score": (
                selector_score_decoy_margin > 0.0
            ),
            "score_margin_expansion_gate_passed": False,
            "recall_frontier_rank": "",
            "native_truth_used_before_selection": False,
            "native_label_attached_after_event_generation": True,
            "diagnostic_claim_allowed": False,
        }
        row.update(_direct_constraint_stats(context=context, event=event))
        rows.append(row)
    return rows


def _region_pair_union(
    events: Sequence[NucleusClosureEvent],
) -> set[tuple[int, int]]:
    region_pairs: set[tuple[int, int]] = set()
    for event in events:
        region_pairs.update(event.candidate_region_pairs())
    return region_pairs


def terminal_bridge_replacement_frontier_rows(
    *,
    context: CouplingNucleusContext,
    terminal_run: TraceLoopRun,
    replacement_probe_run: TraceLoopRun,
) -> list[dict[str, object]]:
    selected_ids = {event.event_id for event in terminal_run.selected_events}
    selected_by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in terminal_run.selected_events:
        selected_by_row.setdefault(event.row_id, []).append(event)
    replacement_probe_added_ids = {
        event.event_id
        for event in replacement_probe_run.selected_events
        if event.event_id not in selected_ids
    }
    row_by_id = {row.row_id: row for row in context.rows}
    constraints_by_row = context.coupling_dataset.constraints_by_row_id()

    rows: list[dict[str, object]] = []
    for event in context.competitive_events:
        if event.event_id in selected_ids:
            continue
        row_selected = selected_by_row.get(event.row_id, [])
        blockers = [
            selected_event
            for selected_event in row_selected
            if not compatible_future_event(selected_event, event)
        ]
        if not blockers:
            continue
        remaining = [
            selected_event
            for selected_event in row_selected
            if selected_event.event_id
            not in {blocker.event_id for blocker in blockers}
        ]
        if any(
            not compatible_future_event(selected_event, event)
            for selected_event in remaining
        ):
            continue
        row = row_by_id[event.row_id]
        native_long = {
            pair for pair in row.native_contact_pairs() if pair[1] - pair[0] >= 24
        }
        current_native_long = _region_pair_union(row_selected) & native_long
        replacement_native_long = _region_pair_union([*remaining, event]) & native_long
        replacement_delta = len(replacement_native_long) - len(current_native_long)
        if replacement_delta <= 0:
            continue
        constraint_confidence_by_pair = {
            constraint.pair(): constraint.confidence
            for constraint in constraints_by_row.get(event.row_id, ())
        }
        external_constraint_pairs = set(constraint_confidence_by_pair)
        current_external_constraints = (
            _region_pair_union(row_selected) & external_constraint_pairs
        )
        replacement_external_constraints = (
            _region_pair_union([*remaining, event]) & external_constraint_pairs
        )
        current_external_confidence = sum(
            constraint_confidence_by_pair[pair]
            for pair in current_external_constraints
        )
        replacement_external_confidence = sum(
            constraint_confidence_by_pair[pair]
            for pair in replacement_external_constraints
        )
        external_count_delta = (
            len(replacement_external_constraints) - len(current_external_constraints)
        )
        external_confidence_delta = _rounded(
            replacement_external_confidence - current_external_confidence
        )

        assessment = context.assessment_by_event_id[event.event_id]
        decoy, selector_score_margin = _matched_selector_score_decoy(
            context=context,
            event=event,
        )
        decoy_assessment = context.assessment_by_event_id[decoy.event_id]
        blockers_native_contact_count = sum(
            blocker.native_contact_count_after_scoring for blocker in blockers
        )
        blockers_native_long_range_contact_count = sum(
            blocker.native_long_range_contact_count_after_scoring
            for blocker in blockers
        )
        blockers_coupling_nucleus_score = _mean(
            [coupling_nucleus_score(blocker, context) for blocker in blockers]
        )
        blockers_direct_support = _mean(
            [
                context.assessment_by_event_id[
                    blocker.event_id
                ].direct_support_score
                for blocker in blockers
            ]
        )
        blockers_future_preservation = _mean(
            [
                context.assessment_by_event_id[
                    blocker.event_id
                ].future_preservation_score
                for blocker in blockers
            ]
        )
        blockers_contact_cluster_gain = _mean(
            [blocker.contact_cluster_gain for blocker in blockers]
        )
        probe_selected = event.event_id in replacement_probe_added_ids
        reasons = ["requires_replacement"]
        if not probe_selected:
            reasons.append("boundary_field_dense_continuity_probe_not_selected")
        rows.append(
            {
                "frontier_kind": "terminal_bridge_replacement_frontier_v0",
                "source_selector": terminal_run.selector_name,
                "target_selector": replacement_probe_run.selector_name,
                "row_id": event.row_id,
                "source_accession": event.source_accession,
                "event_id": event.event_id,
                "segment_a_start": event.segment_a_start,
                "segment_a_end": event.segment_a_end,
                "segment_b_start": event.segment_b_start,
                "segment_b_end": event.segment_b_end,
                "native_contact_count_after_scoring": (
                    event.native_contact_count_after_scoring
                ),
                "native_long_range_contact_count_after_scoring": (
                    event.native_long_range_contact_count_after_scoring
                ),
                "contact_cluster_gain": event.contact_cluster_gain,
                "coupling_selectivity_score": assessment.coupling_selectivity_score,
                "direct_support_score": assessment.direct_support_score,
                "future_preservation_score": assessment.future_preservation_score,
                "blocked_future_pressure": assessment.blocked_future_pressure,
                "coupling_decoy_margin": (
                    context.coupling_decoy_margin_by_event_id[event.event_id]
                ),
                "coupling_nucleus_score": coupling_nucleus_score(event, context),
                "selector_score_decoy_margin": selector_score_margin,
                "exclusion_reasons": ";".join(reasons),
                "matched_decoy_event_id": decoy.event_id,
                "matched_decoy_native_positive_after_scoring": (
                    decoy.native_contact_count_after_scoring > 0
                ),
                "matched_decoy_contact_cluster_gain": decoy.contact_cluster_gain,
                "matched_decoy_coupling_selectivity_score": (
                    decoy_assessment.coupling_selectivity_score
                ),
                "matched_decoy_direct_support_score": (
                    decoy_assessment.direct_support_score
                ),
                "matched_decoy_future_preservation_score": (
                    decoy_assessment.future_preservation_score
                ),
                "matched_decoy_blocked_future_pressure": (
                    decoy_assessment.blocked_future_pressure
                ),
                "matched_decoy_coupling_nucleus_score": coupling_nucleus_score(
                    decoy,
                    context,
                ),
                "real_beats_decoy_by_coupling_score": (
                    assessment.coupling_selectivity_score
                    > decoy_assessment.coupling_selectivity_score
                ),
                "real_beats_decoy_by_coupling_nucleus_score": (
                    selector_score_margin > 0.0
                ),
                "score_margin_expansion_gate_passed": False,
                "recall_frontier_rank": "",
                "native_truth_used_before_selection": False,
                "native_label_attached_after_event_generation": True,
                "diagnostic_claim_allowed": False,
                "blocking_event_ids": ";".join(
                    blocker.event_id for blocker in blockers
                ),
                "blocking_event_count": len(blockers),
                "blocking_native_contact_count_after_scoring": (
                    blockers_native_contact_count
                ),
                "blocking_native_long_range_contact_count_after_scoring": (
                    blockers_native_long_range_contact_count
                ),
                "blocking_mean_coupling_nucleus_score": (
                    blockers_coupling_nucleus_score
                ),
                "blocking_mean_direct_support_score": blockers_direct_support,
                "blocking_mean_future_preservation_score": (
                    blockers_future_preservation
                ),
                "blocking_mean_contact_cluster_gain": blockers_contact_cluster_gain,
                "replacement_native_long_range_delta_after_scoring": (
                    replacement_delta
                ),
                "current_external_constraint_coverage_count": (
                    len(current_external_constraints)
                ),
                "replacement_external_constraint_coverage_count": (
                    len(replacement_external_constraints)
                ),
                "replacement_external_constraint_coverage_delta": (
                    external_count_delta
                ),
                "current_external_constraint_coverage_confidence": _rounded(
                    current_external_confidence
                ),
                "replacement_external_constraint_coverage_confidence": _rounded(
                    replacement_external_confidence
                ),
                "replacement_external_constraint_confidence_delta": (
                    external_confidence_delta
                ),
                "replacement_external_constraint_coverage_improved": (
                    external_count_delta > 0
                ),
                "replacement_external_constraint_confidence_improved": (
                    external_confidence_delta > 0.0
                ),
                "replacement_probe_selected": probe_selected,
            }
        )
        rows[-1].update(_direct_constraint_stats(context=context, event=event))

    return sorted(
        rows,
        key=lambda row: (
            -int(row["replacement_native_long_range_delta_after_scoring"]),
            -int(row["native_long_range_contact_count_after_scoring"]),
            str(row["row_id"]),
            int(row["segment_a_start"]),
            int(row["segment_b_start"]),
            str(row["event_id"]),
        ),
    )


def _controls_from_runs(runs: Sequence[TraceLoopRun]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for run in runs:
        row = run.metric.to_dict()
        row["control_name"] = run.selector_name
        row["control_kind"] = run.control_kind
        row["constraint_count"] = run.constraint_count
        row["negative_control"] = (
            run.selector_name in MATCHED_CONTROL_NAMES
            or run.control_kind == NEGATIVE_CONTROL_KIND
            or run.control_kind == ADVERSARIAL_CALIBRATED_CONTROL_KIND
        )
        rows.append(row)
    return rows


def _score_margin_expansion_summary(
    rows: Sequence[Mapping[str, object]],
) -> dict[str, int]:
    passed = [
        row
        for row in rows
        if bool(row.get("score_margin_expansion_gate_passed", False))
    ]
    passed_row_ids = {str(row["row_id"]) for row in passed}
    native_positive_row_ids = {
        str(row["row_id"])
        for row in passed
        if int(row["native_contact_count_after_scoring"]) > 0
    }
    native_long_range_row_ids = {
        str(row["row_id"])
        for row in passed
        if int(row["native_long_range_contact_count_after_scoring"]) > 0
    }
    return {
        "candidate_count": len(passed),
        "row_count": len(passed_row_ids),
        "native_positive_row_count": len(native_positive_row_ids),
        "native_long_range_row_count": len(native_long_range_row_ids),
        "native_contact_count": sum(
            int(row["native_contact_count_after_scoring"]) for row in passed
        ),
        "native_long_range_contact_count": sum(
            int(row["native_long_range_contact_count_after_scoring"])
            for row in passed
        ),
        "false_candidate_count": sum(
            1
            for row in passed
            if int(row["native_contact_count_after_scoring"]) == 0
        ),
    }


def _max_summary_value(
    summaries: Sequence[Mapping[str, int]],
    field_name: str,
) -> int:
    return max((int(summary[field_name]) for summary in summaries), default=0)


def classify_external_probe_result(
    *,
    available_rows: int,
    external_real_beats_physical: bool,
    external_real_beats_matched_controls: bool,
    external_constraint_count: Optional[int] = None,
) -> str:
    if external_constraint_count == 0 or available_rows == 0:
        return "no_external_data_built"
    if available_rows < 4:
        return "insufficient_external_signal"
    if external_real_beats_physical and external_real_beats_matched_controls:
        return "external_channel_supported_in_v0"
    return "external_channel_not_yet_supported"


def _external_probe_reason(
    *,
    result: str,
    selected_event_count: int,
    external_constraint_count: int,
) -> str:
    if result == "no_external_data_built":
        return "no MSA/DCA couplings were generated"
    if selected_event_count == 0:
        return "external couplings produced no selected trace-loop events"
    if result == "insufficient_external_signal":
        return "fewer than four rows have accepted external couplings"
    if result == "external_channel_supported_in_v0":
        return "accepted external couplings beat physical and matched controls"
    if external_constraint_count == 0:
        return "no accepted external constraints reached the selector"
    return "accepted external couplings did not beat required controls"


def _certificate(report: Mapping[str, object]) -> dict[str, object]:
    return {
        "certificate_kind": EXTERNAL_COUPLING_TRACE_LOOP_CERTIFICATE_KIND,
        "report_kind": report["report_kind"],
        "batch_id": report["batch_id"],
        "result": report["result"],
        "external_probe_passed": report["external_probe_passed"],
        "external_couplings_available_rows": report[
            "external_couplings_available_rows"
        ],
        "external_rows_rejected_low_depth": report[
            "external_rows_rejected_low_depth"
        ],
        "external_rows_rejected_mapping": report["external_rows_rejected_mapping"],
        "external_real_beats_physical": report["external_real_beats_physical"],
        "external_real_beats_matched_controls": report[
            "external_real_beats_matched_controls"
        ],
        "external_margin_gated_beats_matched_controls": report[
            "external_margin_gated_beats_matched_controls"
        ],
        "external_top_rank_gated_beats_matched_controls": report[
            "external_top_rank_gated_beats_matched_controls"
        ],
        "external_core_expanded_beats_matched_controls": report[
            "external_core_expanded_beats_matched_controls"
        ],
        "external_cluster_gated_core_expanded_beats_matched_controls": report[
            "external_cluster_gated_core_expanded_beats_matched_controls"
        ],
        "external_rank_consistent_cluster_gated_beats_matched_controls": report[
            "external_rank_consistent_cluster_gated_beats_matched_controls"
        ],
        "external_rank_consistent_cluster_gated_beats_adversarial_calibrated_controls": report[
            "external_rank_consistent_cluster_gated_beats_adversarial_calibrated_controls"
        ],
        "external_persistent_rank_consistent_cluster_gated_beats_matched_controls": report[
            "external_persistent_rank_consistent_cluster_gated_beats_matched_controls"
        ],
        "external_persistent_rank_consistent_cluster_gated_beats_adversarial_calibrated_controls": report[
            "external_persistent_rank_consistent_cluster_gated_beats_adversarial_calibrated_controls"
        ],
        "external_persistent_rank_consistent_cluster_gated_selector_score_probe_passed": report[
            "external_persistent_rank_consistent_cluster_gated_selector_score_probe_passed"
        ],
        "external_persistent_rank_consistent_cluster_gated_hard_selector_score_probe_passed": report[
            "external_persistent_rank_consistent_cluster_gated_hard_selector_score_probe_passed"
        ],
        "external_margin_gated_claim_allowed": report[
            "external_margin_gated_claim_allowed"
        ],
        "external_top_rank_gated_claim_allowed": report[
            "external_top_rank_gated_claim_allowed"
        ],
        "external_core_expanded_claim_allowed": report[
            "external_core_expanded_claim_allowed"
        ],
        "external_cluster_gated_core_expanded_claim_allowed": report[
            "external_cluster_gated_core_expanded_claim_allowed"
        ],
        "external_rank_consistent_cluster_gated_claim_allowed": report[
            "external_rank_consistent_cluster_gated_claim_allowed"
        ],
        "external_persistent_rank_consistent_cluster_gated_claim_allowed": report[
            "external_persistent_rank_consistent_cluster_gated_claim_allowed"
        ],
        "external_score_margin_expanded_selected_event_count": report[
            "external_score_margin_expanded_selected_event_count"
        ],
        "external_score_margin_expanded_added_event_count": report[
            "external_score_margin_expanded_added_event_count"
        ],
        "external_score_margin_expanded_added_native_long_range_contact_count": report[
            "external_score_margin_expanded_added_native_long_range_contact_count"
        ],
        "external_score_margin_expanded_added_false_event_count": report[
            "external_score_margin_expanded_added_false_event_count"
        ],
        "external_score_margin_expanded_false_nucleus_rate": report[
            "external_score_margin_expanded_false_nucleus_rate"
        ],
        "external_score_margin_expanded_long_range_recall": report[
            "external_score_margin_expanded_long_range_recall"
        ],
        "external_score_margin_expanded_long_range_recall_delta_vs_persistent": report[
            "external_score_margin_expanded_long_range_recall_delta_vs_persistent"
        ],
        "external_score_margin_expanded_beats_matched_controls": report[
            "external_score_margin_expanded_beats_matched_controls"
        ],
        "external_score_margin_expanded_beats_adversarial_calibrated_controls": report[
            "external_score_margin_expanded_beats_adversarial_calibrated_controls"
        ],
        "external_score_margin_expanded_claim_allowed": report[
            "external_score_margin_expanded_claim_allowed"
        ],
        "external_boundary_continuity_expanded_selected_event_count": report[
            "external_boundary_continuity_expanded_selected_event_count"
        ],
        "external_boundary_continuity_expanded_added_event_count": report[
            "external_boundary_continuity_expanded_added_event_count"
        ],
        "external_boundary_continuity_expanded_added_native_long_range_contact_count": report[
            "external_boundary_continuity_expanded_added_native_long_range_contact_count"
        ],
        "external_boundary_continuity_expanded_added_false_event_count": report[
            "external_boundary_continuity_expanded_added_false_event_count"
        ],
        "external_boundary_continuity_expanded_false_nucleus_rate": report[
            "external_boundary_continuity_expanded_false_nucleus_rate"
        ],
        "external_boundary_continuity_expanded_long_range_recall": report[
            "external_boundary_continuity_expanded_long_range_recall"
        ],
        "external_boundary_continuity_expanded_long_range_recall_delta_vs_score_margin": report[
            "external_boundary_continuity_expanded_long_range_recall_delta_vs_score_margin"
        ],
        "external_boundary_continuity_expanded_beats_matched_controls": report[
            "external_boundary_continuity_expanded_beats_matched_controls"
        ],
        "external_boundary_continuity_expanded_beats_adversarial_calibrated_controls": report[
            "external_boundary_continuity_expanded_beats_adversarial_calibrated_controls"
        ],
        "external_boundary_continuity_expanded_claim_allowed": report[
            "external_boundary_continuity_expanded_claim_allowed"
        ],
        "external_edge_continuity_expanded_selected_event_count": report[
            "external_edge_continuity_expanded_selected_event_count"
        ],
        "external_edge_continuity_expanded_added_event_count": report[
            "external_edge_continuity_expanded_added_event_count"
        ],
        "external_edge_continuity_expanded_added_native_long_range_contact_count": report[
            "external_edge_continuity_expanded_added_native_long_range_contact_count"
        ],
        "external_edge_continuity_expanded_added_false_event_count": report[
            "external_edge_continuity_expanded_added_false_event_count"
        ],
        "external_edge_continuity_expanded_false_nucleus_rate": report[
            "external_edge_continuity_expanded_false_nucleus_rate"
        ],
        "external_edge_continuity_expanded_long_range_recall": report[
            "external_edge_continuity_expanded_long_range_recall"
        ],
        "external_edge_continuity_expanded_long_range_recall_delta_vs_boundary_continuity": report[
            "external_edge_continuity_expanded_long_range_recall_delta_vs_boundary_continuity"
        ],
        "external_edge_continuity_expanded_beats_matched_controls": report[
            "external_edge_continuity_expanded_beats_matched_controls"
        ],
        "external_edge_continuity_expanded_beats_adversarial_calibrated_controls": report[
            "external_edge_continuity_expanded_beats_adversarial_calibrated_controls"
        ],
        "external_edge_continuity_expanded_claim_allowed": report[
            "external_edge_continuity_expanded_claim_allowed"
        ],
        "external_pressure_release_expanded_selected_event_count": report[
            "external_pressure_release_expanded_selected_event_count"
        ],
        "external_pressure_release_expanded_added_event_count": report[
            "external_pressure_release_expanded_added_event_count"
        ],
        "external_pressure_release_expanded_added_native_long_range_contact_count": report[
            "external_pressure_release_expanded_added_native_long_range_contact_count"
        ],
        "external_pressure_release_expanded_added_false_event_count": report[
            "external_pressure_release_expanded_added_false_event_count"
        ],
        "external_pressure_release_expanded_false_nucleus_rate": report[
            "external_pressure_release_expanded_false_nucleus_rate"
        ],
        "external_pressure_release_expanded_long_range_recall": report[
            "external_pressure_release_expanded_long_range_recall"
        ],
        "external_pressure_release_expanded_long_range_recall_delta_vs_edge_continuity": report[
            "external_pressure_release_expanded_long_range_recall_delta_vs_edge_continuity"
        ],
        "external_pressure_release_expanded_beats_matched_controls": report[
            "external_pressure_release_expanded_beats_matched_controls"
        ],
        "external_pressure_release_expanded_beats_adversarial_calibrated_controls": report[
            "external_pressure_release_expanded_beats_adversarial_calibrated_controls"
        ],
        "external_pressure_release_expanded_claim_allowed": report[
            "external_pressure_release_expanded_claim_allowed"
        ],
        "external_registry_extension_expanded_selected_event_count": report[
            "external_registry_extension_expanded_selected_event_count"
        ],
        "external_registry_extension_expanded_added_event_count": report[
            "external_registry_extension_expanded_added_event_count"
        ],
        "external_registry_extension_expanded_added_native_long_range_contact_count": report[
            "external_registry_extension_expanded_added_native_long_range_contact_count"
        ],
        "external_registry_extension_expanded_added_false_event_count": report[
            "external_registry_extension_expanded_added_false_event_count"
        ],
        "external_registry_extension_expanded_false_nucleus_rate": report[
            "external_registry_extension_expanded_false_nucleus_rate"
        ],
        "external_registry_extension_expanded_long_range_recall": report[
            "external_registry_extension_expanded_long_range_recall"
        ],
        "external_registry_extension_expanded_long_range_recall_delta_vs_pressure_release": report[
            "external_registry_extension_expanded_long_range_recall_delta_vs_pressure_release"
        ],
        "external_registry_extension_expanded_beats_matched_controls": report[
            "external_registry_extension_expanded_beats_matched_controls"
        ],
        "external_registry_extension_expanded_beats_adversarial_calibrated_controls": report[
            "external_registry_extension_expanded_beats_adversarial_calibrated_controls"
        ],
        "external_registry_extension_expanded_claim_allowed": report[
            "external_registry_extension_expanded_claim_allowed"
        ],
        "external_terminal_bridge_expanded_selected_event_count": report[
            "external_terminal_bridge_expanded_selected_event_count"
        ],
        "external_terminal_bridge_expanded_added_event_count": report[
            "external_terminal_bridge_expanded_added_event_count"
        ],
        "external_terminal_bridge_expanded_added_native_long_range_contact_count": report[
            "external_terminal_bridge_expanded_added_native_long_range_contact_count"
        ],
        "external_terminal_bridge_expanded_added_false_event_count": report[
            "external_terminal_bridge_expanded_added_false_event_count"
        ],
        "external_terminal_bridge_expanded_false_nucleus_rate": report[
            "external_terminal_bridge_expanded_false_nucleus_rate"
        ],
        "external_terminal_bridge_expanded_long_range_recall": report[
            "external_terminal_bridge_expanded_long_range_recall"
        ],
        "external_terminal_bridge_expanded_long_range_recall_delta_vs_registry_extension": report[
            "external_terminal_bridge_expanded_long_range_recall_delta_vs_registry_extension"
        ],
        "external_terminal_bridge_expanded_beats_matched_controls": report[
            "external_terminal_bridge_expanded_beats_matched_controls"
        ],
        "external_terminal_bridge_expanded_beats_adversarial_calibrated_controls": report[
            "external_terminal_bridge_expanded_beats_adversarial_calibrated_controls"
        ],
        "external_terminal_bridge_expanded_claim_allowed": report[
            "external_terminal_bridge_expanded_claim_allowed"
        ],
        "external_boundary_field_replacement_probe_selected_event_count": report[
            "external_boundary_field_replacement_probe_selected_event_count"
        ],
        "external_boundary_field_replacement_probe_added_event_count": report[
            "external_boundary_field_replacement_probe_added_event_count"
        ],
        "external_boundary_field_replacement_probe_false_nucleus_rate": report[
            "external_boundary_field_replacement_probe_false_nucleus_rate"
        ],
        "external_boundary_field_replacement_probe_long_range_recall": report[
            "external_boundary_field_replacement_probe_long_range_recall"
        ],
        "external_boundary_field_replacement_probe_long_range_recall_delta_vs_terminal_bridge": report[
            "external_boundary_field_replacement_probe_long_range_recall_delta_vs_terminal_bridge"
        ],
        "external_boundary_field_replacement_probe_claim_allowed": report[
            "external_boundary_field_replacement_probe_claim_allowed"
        ],
        "external_macro_scale_future_preserved_segment_length": report[
            "external_macro_scale_future_preserved_segment_length"
        ],
        "external_macro_scale_future_preserved_segment_stride": report[
            "external_macro_scale_future_preserved_segment_stride"
        ],
        "external_macro_scale_future_preserved_selected_event_count": report[
            "external_macro_scale_future_preserved_selected_event_count"
        ],
        "external_macro_scale_future_preserved_false_nucleus_rate": report[
            "external_macro_scale_future_preserved_false_nucleus_rate"
        ],
        "external_macro_scale_future_preserved_cluster_precision": report[
            "external_macro_scale_future_preserved_cluster_precision"
        ],
        "external_macro_scale_future_preserved_long_range_recall": report[
            "external_macro_scale_future_preserved_long_range_recall"
        ],
        "external_macro_scale_future_preserved_long_range_recall_delta_vs_boundary_replacement": report[
            "external_macro_scale_future_preserved_long_range_recall_delta_vs_boundary_replacement"
        ],
        "external_macro_scale_future_preserved_cluster_precision_margin_vs_matched_controls": report[
            "external_macro_scale_future_preserved_cluster_precision_margin_vs_matched_controls"
        ],
        "external_macro_scale_future_preserved_long_range_recall_margin_vs_matched_controls": report[
            "external_macro_scale_future_preserved_long_range_recall_margin_vs_matched_controls"
        ],
        "external_macro_scale_future_preserved_cluster_precision_margin_vs_adversarial_controls": report[
            "external_macro_scale_future_preserved_cluster_precision_margin_vs_adversarial_controls"
        ],
        "external_macro_scale_future_preserved_long_range_recall_margin_vs_adversarial_controls": report[
            "external_macro_scale_future_preserved_long_range_recall_margin_vs_adversarial_controls"
        ],
        "external_macro_scale_future_preserved_beats_matched_controls": report[
            "external_macro_scale_future_preserved_beats_matched_controls"
        ],
        "external_macro_scale_future_preserved_beats_adversarial_calibrated_controls": report[
            "external_macro_scale_future_preserved_beats_adversarial_calibrated_controls"
        ],
        "external_macro_scale_future_preserved_claim_allowed": report[
            "external_macro_scale_future_preserved_claim_allowed"
        ],
        "external_multiscale_future_preserved_segment_lengths": report[
            "external_multiscale_future_preserved_segment_lengths"
        ],
        "external_multiscale_future_preserved_max_events_per_row": report[
            "external_multiscale_future_preserved_max_events_per_row"
        ],
        "external_multiscale_future_preserved_selected_event_count": report[
            "external_multiscale_future_preserved_selected_event_count"
        ],
        "external_multiscale_future_preserved_false_nucleus_rate": report[
            "external_multiscale_future_preserved_false_nucleus_rate"
        ],
        "external_multiscale_future_preserved_cluster_precision": report[
            "external_multiscale_future_preserved_cluster_precision"
        ],
        "external_multiscale_future_preserved_long_range_recall": report[
            "external_multiscale_future_preserved_long_range_recall"
        ],
        "external_multiscale_future_preserved_long_range_recall_delta_vs_macro": report[
            "external_multiscale_future_preserved_long_range_recall_delta_vs_macro"
        ],
        "external_multiscale_future_preserved_long_range_recall_margin_vs_matched_controls": report[
            "external_multiscale_future_preserved_long_range_recall_margin_vs_matched_controls"
        ],
        "external_multiscale_future_preserved_long_range_recall_margin_vs_adversarial_controls": report[
            "external_multiscale_future_preserved_long_range_recall_margin_vs_adversarial_controls"
        ],
        "external_multiscale_future_preserved_beats_matched_controls": report[
            "external_multiscale_future_preserved_beats_matched_controls"
        ],
        "external_multiscale_future_preserved_beats_adversarial_calibrated_controls": report[
            "external_multiscale_future_preserved_beats_adversarial_calibrated_controls"
        ],
        "external_multiscale_future_preserved_claim_allowed": report[
            "external_multiscale_future_preserved_claim_allowed"
        ],
        "external_terminal_bridge_replacement_frontier_count": report[
            "external_terminal_bridge_replacement_frontier_count"
        ],
        "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum": report[
            "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum"
        ],
        "external_terminal_bridge_replacement_frontier_probe_selected_count": report[
            "external_terminal_bridge_replacement_frontier_probe_selected_count"
        ],
        "external_terminal_bridge_replacement_frontier_probe_selected_native_long_range_delta_sum": report[
            "external_terminal_bridge_replacement_frontier_probe_selected_native_long_range_delta_sum"
        ],
        "external_terminal_bridge_replacement_frontier_external_count_delta_positive_count": report[
            "external_terminal_bridge_replacement_frontier_external_count_delta_positive_count"
        ],
        "external_terminal_bridge_replacement_frontier_external_confidence_delta_positive_count": report[
            "external_terminal_bridge_replacement_frontier_external_confidence_delta_positive_count"
        ],
        "external_terminal_bridge_replacement_frontier_external_count_delta_sum": report[
            "external_terminal_bridge_replacement_frontier_external_count_delta_sum"
        ],
        "external_terminal_bridge_replacement_frontier_external_confidence_delta_sum": report[
            "external_terminal_bridge_replacement_frontier_external_confidence_delta_sum"
        ],
        "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum_with_external_count_gain": report[
            "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum_with_external_count_gain"
        ],
        "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum_with_external_confidence_gain": report[
            "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum_with_external_confidence_gain"
        ],
        "external_terminal_bridge_replacement_frontier_claim_allowed": report[
            "external_terminal_bridge_replacement_frontier_claim_allowed"
        ],
        "external_rank_consistent_cluster_gated_native_positive_frontier_count": report[
            "external_rank_consistent_cluster_gated_native_positive_frontier_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_recovered_event_count": report[
            "external_persistent_rank_consistent_cluster_gated_recovered_event_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_recovered_native_contact_count": report[
            "external_persistent_rank_consistent_cluster_gated_recovered_native_contact_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_row_count": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_row_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_native_long_range_contact_count": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_native_long_range_contact_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_false_candidate_count": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_false_candidate_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_candidate_count": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_candidate_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_row_count": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_row_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_native_long_range_contact_count": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_native_long_range_contact_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count_margin_vs_matched_controls": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count_margin_vs_matched_controls"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count_margin_vs_matched_controls": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count_margin_vs_matched_controls"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_matched_controls": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_matched_controls"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_candidate_count": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_candidate_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_row_count": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_row_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_native_long_range_contact_count": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_native_long_range_contact_count"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count_margin_vs_adversarial_controls": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count_margin_vs_adversarial_controls"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count_margin_vs_adversarial_controls": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count_margin_vs_adversarial_controls"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_adversarial_controls": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_adversarial_controls"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_signal_seen": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_signal_seen"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_beats_matched_controls": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_beats_matched_controls"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_beats_adversarial_controls": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_beats_adversarial_controls"
        ],
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_claim_allowed": report[
            "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_claim_allowed"
        ],
        "external_rank_consistent_cluster_gated_frontier_claim_allowed": report[
            "external_rank_consistent_cluster_gated_frontier_claim_allowed"
        ],
        "hard_adversarial_calibrated_probe_passed": report[
            "hard_adversarial_calibrated_probe_passed"
        ],
        "claim_allowed": report["claim_allowed"],
        "mechanism_discovery_claim_allowed": report[
            "mechanism_discovery_claim_allowed"
        ],
        "folding_problem_solved": report["folding_problem_solved"],
        "output_artifacts": tuple(ROOT_OUTPUT_NAMES),
    }


def _build_report(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    import_result: ExternalCouplingImportResult,
    external_real: TraceLoopRun,
    external_margin_gated: TraceLoopRun,
    external_top_rank_gated: TraceLoopRun,
    external_core_expanded: TraceLoopRun,
    external_cluster_gated_core_expanded: TraceLoopRun,
    external_rank_consistent_cluster_gated: TraceLoopRun,
    external_persistent_rank_consistent_cluster_gated: TraceLoopRun,
    external_score_margin_expanded: TraceLoopRun,
    external_boundary_continuity_expanded: TraceLoopRun,
    external_edge_continuity_expanded: TraceLoopRun,
    external_pressure_release_expanded: TraceLoopRun,
    external_registry_extension_expanded: TraceLoopRun,
    external_terminal_bridge_expanded: TraceLoopRun,
    external_boundary_field_replacement_probe: TraceLoopRun,
    external_macro_scale_future_preserved: TraceLoopRun,
    macro_scale_future_preserved_controls: Sequence[TraceLoopRun],
    adversarial_macro_scale_future_preserved_controls: Sequence[TraceLoopRun],
    external_multiscale_future_preserved: TraceLoopRun,
    multiscale_future_preserved_controls: Sequence[TraceLoopRun],
    adversarial_multiscale_future_preserved_controls: Sequence[TraceLoopRun],
    physical_baseline: TraceLoopRun,
    matched_controls: Sequence[TraceLoopRun],
    margin_gated_controls: Sequence[TraceLoopRun],
    top_rank_gated_controls: Sequence[TraceLoopRun],
    core_expanded_controls: Sequence[TraceLoopRun],
    cluster_gated_core_expanded_controls: Sequence[TraceLoopRun],
    rank_consistent_cluster_gated_controls: Sequence[TraceLoopRun],
    adversarial_rank_consistent_controls: Sequence[TraceLoopRun],
    persistent_rank_consistent_controls: Sequence[TraceLoopRun],
    adversarial_persistent_rank_consistent_controls: Sequence[TraceLoopRun],
    score_margin_expanded_controls: Sequence[TraceLoopRun],
    adversarial_score_margin_expanded_controls: Sequence[TraceLoopRun],
    boundary_continuity_expanded_controls: Sequence[TraceLoopRun],
    adversarial_boundary_continuity_expanded_controls: Sequence[TraceLoopRun],
    edge_continuity_expanded_controls: Sequence[TraceLoopRun],
    adversarial_edge_continuity_expanded_controls: Sequence[TraceLoopRun],
    pressure_release_expanded_controls: Sequence[TraceLoopRun],
    adversarial_pressure_release_expanded_controls: Sequence[TraceLoopRun],
    registry_extension_expanded_controls: Sequence[TraceLoopRun],
    adversarial_registry_extension_expanded_controls: Sequence[TraceLoopRun],
    terminal_bridge_expanded_controls: Sequence[TraceLoopRun],
    adversarial_terminal_bridge_expanded_controls: Sequence[TraceLoopRun],
    oracle_positive_control: TraceLoopRun,
    frontier_rows: Sequence[Mapping[str, object]],
    recall_frontier_rows: Sequence[Mapping[str, object]],
    replacement_frontier_rows: Sequence[Mapping[str, object]],
    matched_control_recall_frontier_summaries: Sequence[Mapping[str, int]],
    adversarial_recall_frontier_summaries: Sequence[Mapping[str, int]],
    source_benchmark_file: Path,
    external_coupling_file: Path,
    oracle_coupling_file: Path,
) -> dict[str, object]:
    statuses = [status.row_external_status for status in import_result.row_statuses]
    available_rows = statuses.count("external_couplings_available")
    rejected_low_depth = statuses.count("external_couplings_rejected_low_depth")
    rejected_low_coverage = statuses.count("external_couplings_rejected_low_coverage")
    rejected_mapping = statuses.count("external_couplings_rejected_mapping_ambiguous")
    rejected_coordinate_taint = statuses.count(
        "external_couplings_rejected_coordinate_taint"
    )
    control_false_rates = [run.metric.false_nucleus_rate for run in matched_controls]
    control_precisions = [
        run.metric.contact_cluster_precision for run in matched_controls
    ]
    control_enrichments = [
        run.metric.real_vs_decoy_coupling_enrichment_ratio
        for run in matched_controls
    ]
    max_control_enrichment = max(control_enrichments) if control_enrichments else 0.0
    external_selected_event_count = external_real.metric.selected_event_count
    external_real_vs_control_enrichment_ratio: Optional[float] = (
        _rounded(
            external_real.metric.real_vs_decoy_coupling_enrichment_ratio
            / max_control_enrichment
        )
        if external_selected_event_count > 0 and max_control_enrichment
        else None
    )
    external_real_beats_physical = (
        external_real.metric.selected_event_count > 0
        and external_real.metric.false_nucleus_rate
        < physical_baseline.metric.false_nucleus_rate
        and external_real.metric.contact_cluster_precision
        > physical_baseline.metric.contact_cluster_precision
    )
    external_real_beats_matched_controls = (
        bool(matched_controls)
        and external_real.metric.selected_event_count > 0
        and all(
            external_real.metric.false_nucleus_rate < value
            for value in control_false_rates
        )
        and all(
            external_real.metric.contact_cluster_precision > value
            for value in control_precisions
        )
    )
    margin_control_false_rates = [
        run.metric.false_nucleus_rate for run in margin_gated_controls
    ]
    margin_control_precisions = [
        run.metric.contact_cluster_precision for run in margin_gated_controls
    ]
    margin_control_enrichments = [
        run.metric.real_vs_decoy_coupling_enrichment_ratio
        for run in margin_gated_controls
    ]
    max_margin_control_enrichment = (
        max(margin_control_enrichments) if margin_control_enrichments else 0.0
    )
    margin_selected_event_count = external_margin_gated.metric.selected_event_count
    margin_vs_control_enrichment_ratio: Optional[float] = (
        _rounded(
            external_margin_gated.metric.real_vs_decoy_coupling_enrichment_ratio
            / max_margin_control_enrichment
        )
        if margin_selected_event_count > 0 and max_margin_control_enrichment
        else None
    )
    margin_beats_physical = (
        margin_selected_event_count > 0
        and external_margin_gated.metric.false_nucleus_rate
        < physical_baseline.metric.false_nucleus_rate
        and external_margin_gated.metric.contact_cluster_precision
        > physical_baseline.metric.contact_cluster_precision
    )
    margin_beats_matched_controls = (
        bool(margin_gated_controls)
        and margin_selected_event_count > 0
        and all(
            external_margin_gated.metric.false_nucleus_rate < value
            for value in margin_control_false_rates
        )
        and all(
            external_margin_gated.metric.contact_cluster_precision > value
            for value in margin_control_precisions
        )
    )
    top_rank_control_false_rates = [
        run.metric.false_nucleus_rate for run in top_rank_gated_controls
    ]
    top_rank_control_precisions = [
        run.metric.contact_cluster_precision for run in top_rank_gated_controls
    ]
    top_rank_control_enrichments = [
        run.metric.real_vs_decoy_coupling_enrichment_ratio
        for run in top_rank_gated_controls
    ]
    max_top_rank_control_enrichment = (
        max(top_rank_control_enrichments) if top_rank_control_enrichments else 0.0
    )
    top_rank_selected_event_count = external_top_rank_gated.metric.selected_event_count
    top_rank_vs_control_enrichment_ratio: Optional[float] = (
        _rounded(
            external_top_rank_gated.metric.real_vs_decoy_coupling_enrichment_ratio
            / max_top_rank_control_enrichment
        )
        if top_rank_selected_event_count > 0 and max_top_rank_control_enrichment
        else None
    )
    top_rank_beats_physical = (
        top_rank_selected_event_count > 0
        and external_top_rank_gated.metric.false_nucleus_rate
        < physical_baseline.metric.false_nucleus_rate
        and external_top_rank_gated.metric.contact_cluster_precision
        > physical_baseline.metric.contact_cluster_precision
    )
    top_rank_beats_matched_controls = (
        bool(top_rank_gated_controls)
        and top_rank_selected_event_count > 0
        and all(
            external_top_rank_gated.metric.false_nucleus_rate < value
            for value in top_rank_control_false_rates
        )
        and all(
            external_top_rank_gated.metric.contact_cluster_precision > value
            for value in top_rank_control_precisions
        )
    )
    core_expanded_control_false_rates = [
        run.metric.false_nucleus_rate for run in core_expanded_controls
    ]
    core_expanded_control_precisions = [
        run.metric.contact_cluster_precision for run in core_expanded_controls
    ]
    core_expanded_control_enrichments = [
        run.metric.real_vs_decoy_coupling_enrichment_ratio
        for run in core_expanded_controls
    ]
    max_core_expanded_control_enrichment = (
        max(core_expanded_control_enrichments)
        if core_expanded_control_enrichments
        else 0.0
    )
    core_expanded_selected_event_count = (
        external_core_expanded.metric.selected_event_count
    )
    core_expanded_vs_control_enrichment_ratio: Optional[float] = (
        _rounded(
            external_core_expanded.metric.real_vs_decoy_coupling_enrichment_ratio
            / max_core_expanded_control_enrichment
        )
        if core_expanded_selected_event_count > 0
        and max_core_expanded_control_enrichment
        else None
    )
    core_expanded_beats_physical = (
        core_expanded_selected_event_count > 0
        and external_core_expanded.metric.false_nucleus_rate
        < physical_baseline.metric.false_nucleus_rate
        and external_core_expanded.metric.contact_cluster_precision
        > physical_baseline.metric.contact_cluster_precision
    )
    core_expanded_beats_matched_controls = (
        bool(core_expanded_controls)
        and core_expanded_selected_event_count > 0
        and all(
            external_core_expanded.metric.false_nucleus_rate < value
            for value in core_expanded_control_false_rates
        )
        and all(
            external_core_expanded.metric.contact_cluster_precision > value
            for value in core_expanded_control_precisions
        )
    )
    cluster_gated_control_false_rates = [
        run.metric.false_nucleus_rate
        for run in cluster_gated_core_expanded_controls
    ]
    cluster_gated_control_precisions = [
        run.metric.contact_cluster_precision
        for run in cluster_gated_core_expanded_controls
    ]
    cluster_gated_control_enrichments = [
        run.metric.real_vs_decoy_coupling_enrichment_ratio
        for run in cluster_gated_core_expanded_controls
    ]
    cluster_metric = external_cluster_gated_core_expanded.metric
    max_cluster_gated_control_enrichment = (
        max(cluster_gated_control_enrichments)
        if cluster_gated_control_enrichments
        else 0.0
    )
    cluster_gated_selected_event_count = (
        external_cluster_gated_core_expanded.metric.selected_event_count
    )
    cluster_gated_vs_control_enrichment_ratio: Optional[float] = (
        _rounded(
            external_cluster_gated_core_expanded.metric.real_vs_decoy_coupling_enrichment_ratio
            / max_cluster_gated_control_enrichment
        )
        if cluster_gated_selected_event_count > 0
        and max_cluster_gated_control_enrichment
        else None
    )
    cluster_gated_beats_physical = (
        cluster_gated_selected_event_count > 0
        and external_cluster_gated_core_expanded.metric.false_nucleus_rate
        < physical_baseline.metric.false_nucleus_rate
        and external_cluster_gated_core_expanded.metric.contact_cluster_precision
        > physical_baseline.metric.contact_cluster_precision
    )
    cluster_gated_beats_matched_controls = (
        bool(cluster_gated_core_expanded_controls)
        and cluster_gated_selected_event_count > 0
        and all(
            external_cluster_gated_core_expanded.metric.false_nucleus_rate < value
            for value in cluster_gated_control_false_rates
        )
        and all(
            external_cluster_gated_core_expanded.metric.contact_cluster_precision
            > value
            for value in cluster_gated_control_precisions
        )
    )
    oracle_recall_floor = _rounded(
        0.50 * oracle_positive_control.metric.long_range_contact_recall
    )
    external_real_meets_oracle_recall_floor = (
        external_selected_event_count > 0
        and external_real.metric.long_range_contact_recall >= oracle_recall_floor
    )
    margin_meets_oracle_recall_floor = (
        margin_selected_event_count > 0
        and external_margin_gated.metric.long_range_contact_recall >= oracle_recall_floor
    )
    top_rank_meets_oracle_recall_floor = (
        top_rank_selected_event_count > 0
        and external_top_rank_gated.metric.long_range_contact_recall
        >= oracle_recall_floor
    )
    core_expanded_meets_oracle_recall_floor = (
        core_expanded_selected_event_count > 0
        and external_core_expanded.metric.long_range_contact_recall
        >= oracle_recall_floor
    )
    cluster_gated_meets_oracle_recall_floor = (
        cluster_gated_selected_event_count > 0
        and external_cluster_gated_core_expanded.metric.long_range_contact_recall
        >= oracle_recall_floor
    )
    rank_consistent_control_enrichments = [
        run.metric.real_vs_decoy_coupling_enrichment_ratio
        for run in rank_consistent_cluster_gated_controls
        if run.metric.selected_event_count > 0
    ]
    rank_consistent_selected_event_count = (
        external_rank_consistent_cluster_gated.metric.selected_event_count
    )
    max_rank_consistent_control_enrichment = (
        max(rank_consistent_control_enrichments)
        if rank_consistent_control_enrichments
        else 0.0
    )
    rank_consistent_vs_control_enrichment_ratio: Optional[float] = (
        _rounded(
            external_rank_consistent_cluster_gated.metric.real_vs_decoy_coupling_enrichment_ratio
            / max_rank_consistent_control_enrichment
        )
        if rank_consistent_selected_event_count > 0
        and max_rank_consistent_control_enrichment
        else None
    )
    rank_consistent_beats_physical = (
        rank_consistent_selected_event_count > 0
        and external_rank_consistent_cluster_gated.metric.false_nucleus_rate
        < physical_baseline.metric.false_nucleus_rate
        and external_rank_consistent_cluster_gated.metric.contact_cluster_precision
        > physical_baseline.metric.contact_cluster_precision
    )
    rank_consistent_beats_matched_controls = (
        bool(rank_consistent_cluster_gated_controls)
        and rank_consistent_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or external_rank_consistent_cluster_gated.metric.false_nucleus_rate
            < run.metric.false_nucleus_rate
            for run in rank_consistent_cluster_gated_controls
        )
        and all(
            run.metric.selected_event_count == 0
            or external_rank_consistent_cluster_gated.metric.contact_cluster_precision
            > run.metric.contact_cluster_precision
            for run in rank_consistent_cluster_gated_controls
        )
    )
    rank_consistent_meets_oracle_recall_floor = (
        rank_consistent_selected_event_count > 0
        and external_rank_consistent_cluster_gated.metric.long_range_contact_recall
        >= oracle_recall_floor
    )
    adversarial_enrichments = [
        run.metric.real_vs_decoy_coupling_enrichment_ratio
        for run in adversarial_rank_consistent_controls
        if run.metric.selected_event_count >= ADVERSARIAL_ENRICHMENT_MIN_SELECTED_EVENTS
    ]
    max_adversarial_enrichment = (
        max(adversarial_enrichments) if adversarial_enrichments else 0.0
    )
    rank_consistent_vs_adversarial_enrichment_ratio: Optional[float] = (
        _rounded(
            external_rank_consistent_cluster_gated.metric.real_vs_decoy_coupling_enrichment_ratio
            / max_adversarial_enrichment
        )
        if rank_consistent_selected_event_count > 0 and max_adversarial_enrichment
        else None
    )
    rank_consistent_beats_adversarial_calibrated_controls = (
        bool(adversarial_rank_consistent_controls)
        and rank_consistent_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or external_rank_consistent_cluster_gated.metric.false_nucleus_rate
            < run.metric.false_nucleus_rate
            for run in adversarial_rank_consistent_controls
        )
        and all(
            run.metric.selected_event_count == 0
            or external_rank_consistent_cluster_gated.metric.contact_cluster_precision
            > run.metric.contact_cluster_precision
            for run in adversarial_rank_consistent_controls
        )
        and all(
            run.metric.selected_event_count == 0
            or external_rank_consistent_cluster_gated.metric.long_range_contact_recall
            > run.metric.long_range_contact_recall
            for run in adversarial_rank_consistent_controls
        )
    )
    adversarial_enrichment_ratio_meets = (
        rank_consistent_vs_adversarial_enrichment_ratio is not None
        and rank_consistent_vs_adversarial_enrichment_ratio > 1.0
    )
    persistent_control_enrichments = [
        run.metric.real_vs_decoy_coupling_enrichment_ratio
        for run in persistent_rank_consistent_controls
        if run.metric.selected_event_count > 0
    ]
    persistent_control_nucleus_enrichments = [
        run.metric.real_vs_decoy_coupling_nucleus_enrichment_ratio
        for run in persistent_rank_consistent_controls
        if run.metric.selected_event_count > 0
    ]
    persistent_selected_event_count = (
        external_persistent_rank_consistent_cluster_gated.metric.selected_event_count
    )
    max_persistent_control_enrichment = (
        max(persistent_control_enrichments)
        if persistent_control_enrichments
        else 0.0
    )
    persistent_vs_control_enrichment_ratio: Optional[float] = (
        _rounded(
            external_persistent_rank_consistent_cluster_gated.metric.real_vs_decoy_coupling_enrichment_ratio
            / max_persistent_control_enrichment
        )
        if persistent_selected_event_count > 0 and max_persistent_control_enrichment
        else None
    )
    max_persistent_control_nucleus_enrichment = (
        max(persistent_control_nucleus_enrichments)
        if persistent_control_nucleus_enrichments
        else 0.0
    )
    persistent_vs_control_nucleus_score_enrichment_ratio: Optional[float] = (
        _rounded(
            external_persistent_rank_consistent_cluster_gated.metric.real_vs_decoy_coupling_nucleus_enrichment_ratio
            / max_persistent_control_nucleus_enrichment
        )
        if persistent_selected_event_count > 0
        and max_persistent_control_nucleus_enrichment
        else None
    )
    persistent_beats_physical = (
        persistent_selected_event_count > 0
        and external_persistent_rank_consistent_cluster_gated.metric.false_nucleus_rate
        < physical_baseline.metric.false_nucleus_rate
        and external_persistent_rank_consistent_cluster_gated.metric.contact_cluster_precision
        > physical_baseline.metric.contact_cluster_precision
    )
    persistent_beats_matched_controls = (
        bool(persistent_rank_consistent_controls)
        and persistent_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or external_persistent_rank_consistent_cluster_gated.metric.false_nucleus_rate
            < run.metric.false_nucleus_rate
            for run in persistent_rank_consistent_controls
        )
        and all(
            run.metric.selected_event_count == 0
            or external_persistent_rank_consistent_cluster_gated.metric.contact_cluster_precision
            > run.metric.contact_cluster_precision
            for run in persistent_rank_consistent_controls
        )
    )
    persistent_meets_oracle_recall_floor = (
        persistent_selected_event_count > 0
        and external_persistent_rank_consistent_cluster_gated.metric.long_range_contact_recall
        >= oracle_recall_floor
    )
    adversarial_persistent_enrichments = [
        run.metric.real_vs_decoy_coupling_enrichment_ratio
        for run in adversarial_persistent_rank_consistent_controls
        if run.metric.selected_event_count >= ADVERSARIAL_ENRICHMENT_MIN_SELECTED_EVENTS
    ]
    adversarial_persistent_nucleus_enrichments = [
        run.metric.real_vs_decoy_coupling_nucleus_enrichment_ratio
        for run in adversarial_persistent_rank_consistent_controls
        if run.metric.selected_event_count >= ADVERSARIAL_ENRICHMENT_MIN_SELECTED_EVENTS
    ]
    max_adversarial_persistent_enrichment = (
        max(adversarial_persistent_enrichments)
        if adversarial_persistent_enrichments
        else 0.0
    )
    max_adversarial_persistent_nucleus_enrichment = (
        max(adversarial_persistent_nucleus_enrichments)
        if adversarial_persistent_nucleus_enrichments
        else 0.0
    )
    persistent_vs_adversarial_enrichment_ratio: Optional[float] = (
        _rounded(
            external_persistent_rank_consistent_cluster_gated.metric.real_vs_decoy_coupling_enrichment_ratio
            / max_adversarial_persistent_enrichment
        )
        if persistent_selected_event_count > 0
        and max_adversarial_persistent_enrichment
        else None
    )
    persistent_vs_adversarial_nucleus_score_enrichment_ratio: Optional[float] = (
        _rounded(
            external_persistent_rank_consistent_cluster_gated.metric.real_vs_decoy_coupling_nucleus_enrichment_ratio
            / max_adversarial_persistent_nucleus_enrichment
        )
        if persistent_selected_event_count > 0
        and max_adversarial_persistent_nucleus_enrichment
        else None
    )
    persistent_beats_adversarial_calibrated_controls = (
        bool(adversarial_persistent_rank_consistent_controls)
        and persistent_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or external_persistent_rank_consistent_cluster_gated.metric.false_nucleus_rate
            < run.metric.false_nucleus_rate
            for run in adversarial_persistent_rank_consistent_controls
        )
        and all(
            run.metric.selected_event_count == 0
            or external_persistent_rank_consistent_cluster_gated.metric.contact_cluster_precision
            > run.metric.contact_cluster_precision
            for run in adversarial_persistent_rank_consistent_controls
        )
        and all(
            run.metric.selected_event_count == 0
            or external_persistent_rank_consistent_cluster_gated.metric.long_range_contact_recall
            > run.metric.long_range_contact_recall
            for run in adversarial_persistent_rank_consistent_controls
        )
    )
    persistent_enrichment_ratio_meets = (
        persistent_vs_control_enrichment_ratio is not None
        and persistent_vs_control_enrichment_ratio > 1.25
    )
    persistent_adversarial_enrichment_ratio_meets = (
        persistent_vs_adversarial_enrichment_ratio is not None
        and persistent_vs_adversarial_enrichment_ratio > 1.0
    )
    persistent_nucleus_score_enrichment_ratio_meets = (
        persistent_vs_control_nucleus_score_enrichment_ratio is not None
        and persistent_vs_control_nucleus_score_enrichment_ratio > 1.25
    )
    persistent_adversarial_nucleus_score_enrichment_ratio_meets = (
        persistent_vs_adversarial_nucleus_score_enrichment_ratio is not None
        and persistent_vs_adversarial_nucleus_score_enrichment_ratio > 1.0
    )
    persistent_recovered_events = tuple(
        event
        for event in external_persistent_rank_consistent_cluster_gated.selected_events
        if event.event_id
        not in {
            selected.event_id
            for selected in external_rank_consistent_cluster_gated.selected_events
        }
    )
    claim_mode_failures = coupling_claim_mode_validation_failures(
        import_result.dataset
    )
    enrichment_ratio_meets = (
        external_real_vs_control_enrichment_ratio is not None
        and external_real_vs_control_enrichment_ratio > 1.25
    )
    rank_consistent_enrichment_ratio_meets = (
        rank_consistent_vs_control_enrichment_ratio is not None
        and rank_consistent_vs_control_enrichment_ratio > 1.25
    )
    acceptance_criteria_met = (
        available_rows >= 4
        and external_real_beats_physical
        and external_real_beats_matched_controls
        and external_real_meets_oracle_recall_floor
        and enrichment_ratio_meets
        and not claim_mode_failures
    )
    rank_consistent_acceptance_criteria_met = (
        available_rows >= 4
        and rank_consistent_beats_physical
        and rank_consistent_beats_matched_controls
        and rank_consistent_meets_oracle_recall_floor
        and rank_consistent_enrichment_ratio_meets
        and not claim_mode_failures
    )
    persistent_acceptance_criteria_met = (
        available_rows >= 4
        and persistent_beats_physical
        and persistent_beats_matched_controls
        and persistent_meets_oracle_recall_floor
        and persistent_enrichment_ratio_meets
        and not claim_mode_failures
    )
    persistent_selector_score_acceptance_criteria_met = (
        available_rows >= 4
        and persistent_beats_physical
        and persistent_beats_matched_controls
        and persistent_meets_oracle_recall_floor
        and persistent_nucleus_score_enrichment_ratio_meets
        and not claim_mode_failures
    )
    hard_adversarial_calibrated_probe_passed = (
        rank_consistent_acceptance_criteria_met
        and rank_consistent_beats_adversarial_calibrated_controls
        and adversarial_enrichment_ratio_meets
    )
    persistent_hard_adversarial_calibrated_probe_passed = (
        persistent_acceptance_criteria_met
        and persistent_beats_adversarial_calibrated_controls
        and persistent_adversarial_enrichment_ratio_meets
    )
    persistent_hard_selector_score_probe_passed = (
        persistent_selector_score_acceptance_criteria_met
        and persistent_beats_adversarial_calibrated_controls
        and persistent_adversarial_nucleus_score_enrichment_ratio_meets
    )
    external_probe_passed = (
        acceptance_criteria_met
        or rank_consistent_acceptance_criteria_met
        or persistent_acceptance_criteria_met
        or persistent_selector_score_acceptance_criteria_met
    )
    result = classify_external_probe_result(
        available_rows=available_rows,
        external_constraint_count=len(import_result.dataset.constraints),
        external_real_beats_physical=(
            external_real_beats_physical
            or rank_consistent_beats_physical
            or persistent_beats_physical
        ),
        external_real_beats_matched_controls=(
            external_real_beats_matched_controls
            or rank_consistent_beats_matched_controls
            or persistent_beats_matched_controls
        ),
    )
    reason = (
        "persistent external couplings beat controls under the selector-score decoy metric"
        if persistent_selector_score_acceptance_criteria_met
        else
        "persistent external couplings recover a calibrated trace event while beating controls"
        if persistent_acceptance_criteria_met
        else
        "provenance-calibrated external couplings beat physical and matched controls"
        if rank_consistent_acceptance_criteria_met
        else _external_probe_reason(
            result=result,
            selected_event_count=external_selected_event_count,
            external_constraint_count=len(import_result.dataset.constraints),
        )
    )
    external_metric_defined = external_selected_event_count > 0
    frontier_native_contact_count = sum(
        int(row["native_contact_count_after_scoring"]) for row in frontier_rows
    )
    frontier_native_long_range_contact_count = sum(
        int(row["native_long_range_contact_count_after_scoring"])
        for row in frontier_rows
    )
    score_margin_expansion_rows = [
        row
        for row in recall_frontier_rows
        if bool(row["score_margin_expansion_gate_passed"])
    ]
    score_margin_expansion_native_contact_count = sum(
        int(row["native_contact_count_after_scoring"])
        for row in score_margin_expansion_rows
    )
    score_margin_expansion_native_long_range_contact_count = sum(
        int(row["native_long_range_contact_count_after_scoring"])
        for row in score_margin_expansion_rows
    )
    score_margin_expansion_false_candidate_count = sum(
        1
        for row in score_margin_expansion_rows
        if int(row["native_contact_count_after_scoring"]) == 0
    )
    score_margin_expansion_row_count = len(
        {str(row["row_id"]) for row in score_margin_expansion_rows}
    )
    score_margin_expansion_native_long_range_row_count = len(
        {
            str(row["row_id"])
            for row in score_margin_expansion_rows
            if int(row["native_long_range_contact_count_after_scoring"]) > 0
        }
    )
    max_matched_expansion_candidate_count = _max_summary_value(
        matched_control_recall_frontier_summaries,
        "candidate_count",
    )
    max_matched_expansion_row_count = _max_summary_value(
        matched_control_recall_frontier_summaries,
        "row_count",
    )
    max_matched_expansion_native_long_range_row_count = _max_summary_value(
        matched_control_recall_frontier_summaries,
        "native_long_range_row_count",
    )
    max_matched_expansion_native_long_range_contact_count = _max_summary_value(
        matched_control_recall_frontier_summaries,
        "native_long_range_contact_count",
    )
    max_matched_expansion_false_candidate_count = _max_summary_value(
        matched_control_recall_frontier_summaries,
        "false_candidate_count",
    )
    max_adversarial_expansion_candidate_count = _max_summary_value(
        adversarial_recall_frontier_summaries,
        "candidate_count",
    )
    max_adversarial_expansion_row_count = _max_summary_value(
        adversarial_recall_frontier_summaries,
        "row_count",
    )
    max_adversarial_expansion_native_long_range_row_count = _max_summary_value(
        adversarial_recall_frontier_summaries,
        "native_long_range_row_count",
    )
    max_adversarial_expansion_native_long_range_contact_count = _max_summary_value(
        adversarial_recall_frontier_summaries,
        "native_long_range_contact_count",
    )
    max_adversarial_expansion_false_candidate_count = _max_summary_value(
        adversarial_recall_frontier_summaries,
        "false_candidate_count",
    )
    persistent_selected_ids = {
        event.event_id
        for event in external_persistent_rank_consistent_cluster_gated.selected_events
    }
    score_margin_expanded_metric = external_score_margin_expanded.metric
    score_margin_expanded_selected_event_count = (
        score_margin_expanded_metric.selected_event_count
    )
    score_margin_expanded_added_events = tuple(
        event
        for event in external_score_margin_expanded.selected_events
        if event.event_id not in persistent_selected_ids
    )
    score_margin_expanded_max_matched_control_precision = _max_selected_metric(
        score_margin_expanded_controls,
        "contact_cluster_precision",
    )
    score_margin_expanded_max_matched_control_long_range_recall = (
        _max_selected_metric(
            score_margin_expanded_controls,
            "long_range_contact_recall",
        )
    )
    score_margin_expanded_max_adversarial_precision = _max_selected_metric(
        adversarial_score_margin_expanded_controls,
        "contact_cluster_precision",
    )
    score_margin_expanded_max_adversarial_long_range_recall = _max_selected_metric(
        adversarial_score_margin_expanded_controls,
        "long_range_contact_recall",
    )
    score_margin_expanded_noninferior_false_rate_vs_matched_controls = (
        bool(score_margin_expanded_controls)
        and score_margin_expanded_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or score_margin_expanded_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in score_margin_expanded_controls
        )
    )
    score_margin_expanded_noninferior_false_rate_vs_adversarial_controls = (
        bool(adversarial_score_margin_expanded_controls)
        and score_margin_expanded_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or score_margin_expanded_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in adversarial_score_margin_expanded_controls
        )
    )
    score_margin_expanded_beats_matched_controls = (
        score_margin_expanded_noninferior_false_rate_vs_matched_controls
        and score_margin_expanded_metric.contact_cluster_precision
        > score_margin_expanded_max_matched_control_precision
        and score_margin_expanded_metric.long_range_contact_recall
        > score_margin_expanded_max_matched_control_long_range_recall
    )
    score_margin_expanded_beats_adversarial_calibrated_controls = (
        score_margin_expanded_noninferior_false_rate_vs_adversarial_controls
        and score_margin_expanded_metric.contact_cluster_precision
        > score_margin_expanded_max_adversarial_precision
        and score_margin_expanded_metric.long_range_contact_recall
        > score_margin_expanded_max_adversarial_long_range_recall
    )
    score_margin_expanded_ids = {
        event.event_id for event in external_score_margin_expanded.selected_events
    }
    boundary_continuity_expanded_metric = (
        external_boundary_continuity_expanded.metric
    )
    boundary_continuity_expanded_selected_event_count = (
        boundary_continuity_expanded_metric.selected_event_count
    )
    boundary_continuity_expanded_added_events = tuple(
        event
        for event in external_boundary_continuity_expanded.selected_events
        if event.event_id not in score_margin_expanded_ids
    )
    boundary_continuity_expanded_max_matched_control_precision = (
        _max_selected_metric(
            boundary_continuity_expanded_controls,
            "contact_cluster_precision",
        )
    )
    boundary_continuity_expanded_max_matched_control_long_range_recall = (
        _max_selected_metric(
            boundary_continuity_expanded_controls,
            "long_range_contact_recall",
        )
    )
    boundary_continuity_expanded_max_adversarial_precision = (
        _max_selected_metric(
            adversarial_boundary_continuity_expanded_controls,
            "contact_cluster_precision",
        )
    )
    boundary_continuity_expanded_max_adversarial_long_range_recall = (
        _max_selected_metric(
            adversarial_boundary_continuity_expanded_controls,
            "long_range_contact_recall",
        )
    )
    boundary_continuity_expanded_noninferior_false_rate_vs_matched_controls = (
        bool(boundary_continuity_expanded_controls)
        and boundary_continuity_expanded_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or boundary_continuity_expanded_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in boundary_continuity_expanded_controls
        )
    )
    boundary_continuity_expanded_noninferior_false_rate_vs_adversarial_controls = (
        bool(adversarial_boundary_continuity_expanded_controls)
        and boundary_continuity_expanded_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or boundary_continuity_expanded_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in adversarial_boundary_continuity_expanded_controls
        )
    )
    boundary_continuity_expanded_beats_matched_controls = (
        boundary_continuity_expanded_noninferior_false_rate_vs_matched_controls
        and boundary_continuity_expanded_metric.contact_cluster_precision
        > boundary_continuity_expanded_max_matched_control_precision
        and boundary_continuity_expanded_metric.long_range_contact_recall
        > boundary_continuity_expanded_max_matched_control_long_range_recall
    )
    boundary_continuity_expanded_beats_adversarial_calibrated_controls = (
        boundary_continuity_expanded_noninferior_false_rate_vs_adversarial_controls
        and boundary_continuity_expanded_metric.contact_cluster_precision
        > boundary_continuity_expanded_max_adversarial_precision
        and boundary_continuity_expanded_metric.long_range_contact_recall
        > boundary_continuity_expanded_max_adversarial_long_range_recall
    )
    boundary_continuity_expanded_ids = {
        event.event_id
        for event in external_boundary_continuity_expanded.selected_events
    }
    edge_continuity_expanded_metric = external_edge_continuity_expanded.metric
    edge_continuity_expanded_selected_event_count = (
        edge_continuity_expanded_metric.selected_event_count
    )
    edge_continuity_expanded_added_events = tuple(
        event
        for event in external_edge_continuity_expanded.selected_events
        if event.event_id not in boundary_continuity_expanded_ids
    )
    edge_continuity_expanded_max_matched_control_precision = _max_selected_metric(
        edge_continuity_expanded_controls,
        "contact_cluster_precision",
    )
    edge_continuity_expanded_max_matched_control_long_range_recall = (
        _max_selected_metric(
            edge_continuity_expanded_controls,
            "long_range_contact_recall",
        )
    )
    edge_continuity_expanded_max_adversarial_precision = _max_selected_metric(
        adversarial_edge_continuity_expanded_controls,
        "contact_cluster_precision",
    )
    edge_continuity_expanded_max_adversarial_long_range_recall = _max_selected_metric(
        adversarial_edge_continuity_expanded_controls,
        "long_range_contact_recall",
    )
    edge_continuity_expanded_noninferior_false_rate_vs_matched_controls = (
        bool(edge_continuity_expanded_controls)
        and edge_continuity_expanded_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or edge_continuity_expanded_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in edge_continuity_expanded_controls
        )
    )
    edge_continuity_expanded_noninferior_false_rate_vs_adversarial_controls = (
        bool(adversarial_edge_continuity_expanded_controls)
        and edge_continuity_expanded_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or edge_continuity_expanded_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in adversarial_edge_continuity_expanded_controls
        )
    )
    edge_continuity_expanded_beats_matched_controls = (
        edge_continuity_expanded_noninferior_false_rate_vs_matched_controls
        and edge_continuity_expanded_metric.contact_cluster_precision
        > edge_continuity_expanded_max_matched_control_precision
        and edge_continuity_expanded_metric.long_range_contact_recall
        > edge_continuity_expanded_max_matched_control_long_range_recall
    )
    edge_continuity_expanded_beats_adversarial_calibrated_controls = (
        edge_continuity_expanded_noninferior_false_rate_vs_adversarial_controls
        and edge_continuity_expanded_metric.contact_cluster_precision
        > edge_continuity_expanded_max_adversarial_precision
        and edge_continuity_expanded_metric.long_range_contact_recall
        > edge_continuity_expanded_max_adversarial_long_range_recall
    )
    edge_continuity_expanded_ids = {
        event.event_id for event in external_edge_continuity_expanded.selected_events
    }
    pressure_release_expanded_metric = external_pressure_release_expanded.metric
    pressure_release_expanded_selected_event_count = (
        pressure_release_expanded_metric.selected_event_count
    )
    pressure_release_expanded_added_events = tuple(
        event
        for event in external_pressure_release_expanded.selected_events
        if event.event_id not in edge_continuity_expanded_ids
    )
    pressure_release_expanded_max_matched_control_precision = _max_selected_metric(
        pressure_release_expanded_controls,
        "contact_cluster_precision",
    )
    pressure_release_expanded_max_matched_control_long_range_recall = (
        _max_selected_metric(
            pressure_release_expanded_controls,
            "long_range_contact_recall",
        )
    )
    pressure_release_expanded_max_adversarial_precision = _max_selected_metric(
        adversarial_pressure_release_expanded_controls,
        "contact_cluster_precision",
    )
    pressure_release_expanded_max_adversarial_long_range_recall = (
        _max_selected_metric(
            adversarial_pressure_release_expanded_controls,
            "long_range_contact_recall",
        )
    )
    pressure_release_expanded_noninferior_false_rate_vs_matched_controls = (
        bool(pressure_release_expanded_controls)
        and pressure_release_expanded_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or pressure_release_expanded_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in pressure_release_expanded_controls
        )
    )
    pressure_release_expanded_noninferior_false_rate_vs_adversarial_controls = (
        bool(adversarial_pressure_release_expanded_controls)
        and pressure_release_expanded_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or pressure_release_expanded_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in adversarial_pressure_release_expanded_controls
        )
    )
    pressure_release_expanded_beats_matched_controls = (
        pressure_release_expanded_noninferior_false_rate_vs_matched_controls
        and pressure_release_expanded_metric.contact_cluster_precision
        > pressure_release_expanded_max_matched_control_precision
        and pressure_release_expanded_metric.long_range_contact_recall
        > pressure_release_expanded_max_matched_control_long_range_recall
    )
    pressure_release_expanded_beats_adversarial_calibrated_controls = (
        pressure_release_expanded_noninferior_false_rate_vs_adversarial_controls
        and pressure_release_expanded_metric.contact_cluster_precision
        > pressure_release_expanded_max_adversarial_precision
        and pressure_release_expanded_metric.long_range_contact_recall
        > pressure_release_expanded_max_adversarial_long_range_recall
    )
    pressure_release_expanded_ids = {
        event.event_id
        for event in external_pressure_release_expanded.selected_events
    }
    registry_extension_expanded_metric = external_registry_extension_expanded.metric
    registry_extension_expanded_selected_event_count = (
        registry_extension_expanded_metric.selected_event_count
    )
    registry_extension_expanded_added_events = tuple(
        event
        for event in external_registry_extension_expanded.selected_events
        if event.event_id not in pressure_release_expanded_ids
    )
    registry_extension_expanded_max_matched_control_precision = _max_selected_metric(
        registry_extension_expanded_controls,
        "contact_cluster_precision",
    )
    registry_extension_expanded_max_matched_control_long_range_recall = (
        _max_selected_metric(
            registry_extension_expanded_controls,
            "long_range_contact_recall",
        )
    )
    registry_extension_expanded_max_adversarial_precision = _max_selected_metric(
        adversarial_registry_extension_expanded_controls,
        "contact_cluster_precision",
    )
    registry_extension_expanded_max_adversarial_long_range_recall = (
        _max_selected_metric(
            adversarial_registry_extension_expanded_controls,
            "long_range_contact_recall",
        )
    )
    registry_extension_expanded_noninferior_false_rate_vs_matched_controls = (
        bool(registry_extension_expanded_controls)
        and registry_extension_expanded_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or registry_extension_expanded_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in registry_extension_expanded_controls
        )
    )
    registry_extension_expanded_noninferior_false_rate_vs_adversarial_controls = (
        bool(adversarial_registry_extension_expanded_controls)
        and registry_extension_expanded_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or registry_extension_expanded_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in adversarial_registry_extension_expanded_controls
        )
    )
    registry_extension_expanded_beats_matched_controls = (
        registry_extension_expanded_noninferior_false_rate_vs_matched_controls
        and registry_extension_expanded_metric.contact_cluster_precision
        > registry_extension_expanded_max_matched_control_precision
        and registry_extension_expanded_metric.long_range_contact_recall
        > registry_extension_expanded_max_matched_control_long_range_recall
    )
    registry_extension_expanded_beats_adversarial_calibrated_controls = (
        registry_extension_expanded_noninferior_false_rate_vs_adversarial_controls
        and registry_extension_expanded_metric.contact_cluster_precision
        > registry_extension_expanded_max_adversarial_precision
        and registry_extension_expanded_metric.long_range_contact_recall
        > registry_extension_expanded_max_adversarial_long_range_recall
    )
    registry_extension_expanded_ids = {
        event.event_id
        for event in external_registry_extension_expanded.selected_events
    }
    terminal_bridge_expanded_metric = external_terminal_bridge_expanded.metric
    terminal_bridge_expanded_selected_event_count = (
        terminal_bridge_expanded_metric.selected_event_count
    )
    terminal_bridge_expanded_added_events = tuple(
        event
        for event in external_terminal_bridge_expanded.selected_events
        if event.event_id not in registry_extension_expanded_ids
    )
    terminal_bridge_expanded_ids = {
        event.event_id for event in external_terminal_bridge_expanded.selected_events
    }
    boundary_field_replacement_probe_metric = (
        external_boundary_field_replacement_probe.metric
    )
    boundary_field_replacement_probe_added_events = tuple(
        event
        for event in external_boundary_field_replacement_probe.selected_events
        if event.event_id not in terminal_bridge_expanded_ids
    )
    macro_scale_future_preserved_metric = (
        external_macro_scale_future_preserved.metric
    )
    macro_scale_future_preserved_selected_event_count = (
        macro_scale_future_preserved_metric.selected_event_count
    )
    macro_scale_future_preserved_max_matched_control_precision = _max_selected_metric(
        macro_scale_future_preserved_controls,
        "contact_cluster_precision",
    )
    macro_scale_future_preserved_max_matched_control_long_range_recall = (
        _max_selected_metric(
            macro_scale_future_preserved_controls,
            "long_range_contact_recall",
        )
    )
    macro_scale_future_preserved_max_adversarial_precision = _max_selected_metric(
        adversarial_macro_scale_future_preserved_controls,
        "contact_cluster_precision",
    )
    macro_scale_future_preserved_max_adversarial_long_range_recall = (
        _max_selected_metric(
            adversarial_macro_scale_future_preserved_controls,
            "long_range_contact_recall",
        )
    )
    macro_scale_future_preserved_noninferior_false_rate_vs_matched_controls = (
        bool(macro_scale_future_preserved_controls)
        and macro_scale_future_preserved_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or macro_scale_future_preserved_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in macro_scale_future_preserved_controls
        )
    )
    macro_scale_future_preserved_noninferior_false_rate_vs_adversarial_controls = (
        bool(adversarial_macro_scale_future_preserved_controls)
        and macro_scale_future_preserved_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or macro_scale_future_preserved_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in adversarial_macro_scale_future_preserved_controls
        )
    )
    macro_scale_future_preserved_beats_matched_controls = (
        macro_scale_future_preserved_noninferior_false_rate_vs_matched_controls
        and macro_scale_future_preserved_metric.contact_cluster_precision
        > macro_scale_future_preserved_max_matched_control_precision
        and macro_scale_future_preserved_metric.long_range_contact_recall
        > macro_scale_future_preserved_max_matched_control_long_range_recall
    )
    macro_scale_future_preserved_beats_adversarial_calibrated_controls = (
        macro_scale_future_preserved_noninferior_false_rate_vs_adversarial_controls
        and macro_scale_future_preserved_metric.contact_cluster_precision
        > macro_scale_future_preserved_max_adversarial_precision
        and macro_scale_future_preserved_metric.long_range_contact_recall
        > macro_scale_future_preserved_max_adversarial_long_range_recall
    )
    multiscale_future_preserved_metric = (
        external_multiscale_future_preserved.metric
    )
    multiscale_future_preserved_selected_event_count = (
        multiscale_future_preserved_metric.selected_event_count
    )
    multiscale_future_preserved_max_matched_control_precision = (
        _max_selected_metric(
            multiscale_future_preserved_controls,
            "contact_cluster_precision",
        )
    )
    multiscale_future_preserved_max_matched_control_long_range_recall = (
        _max_selected_metric(
            multiscale_future_preserved_controls,
            "long_range_contact_recall",
        )
    )
    multiscale_future_preserved_max_adversarial_precision = _max_selected_metric(
        adversarial_multiscale_future_preserved_controls,
        "contact_cluster_precision",
    )
    multiscale_future_preserved_max_adversarial_long_range_recall = (
        _max_selected_metric(
            adversarial_multiscale_future_preserved_controls,
            "long_range_contact_recall",
        )
    )
    multiscale_future_preserved_noninferior_false_rate_vs_matched_controls = (
        bool(multiscale_future_preserved_controls)
        and multiscale_future_preserved_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or multiscale_future_preserved_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in multiscale_future_preserved_controls
        )
    )
    multiscale_future_preserved_noninferior_false_rate_vs_adversarial_controls = (
        bool(adversarial_multiscale_future_preserved_controls)
        and multiscale_future_preserved_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or multiscale_future_preserved_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in adversarial_multiscale_future_preserved_controls
        )
    )
    multiscale_future_preserved_beats_matched_controls = (
        multiscale_future_preserved_noninferior_false_rate_vs_matched_controls
        and multiscale_future_preserved_metric.contact_cluster_precision
        > multiscale_future_preserved_max_matched_control_precision
        and multiscale_future_preserved_metric.long_range_contact_recall
        > multiscale_future_preserved_max_matched_control_long_range_recall
    )
    multiscale_future_preserved_beats_adversarial_calibrated_controls = (
        multiscale_future_preserved_noninferior_false_rate_vs_adversarial_controls
        and multiscale_future_preserved_metric.contact_cluster_precision
        > multiscale_future_preserved_max_adversarial_precision
        and multiscale_future_preserved_metric.long_range_contact_recall
        > multiscale_future_preserved_max_adversarial_long_range_recall
    )
    replacement_frontier_native_long_range_delta_sum = sum(
        int(row["replacement_native_long_range_delta_after_scoring"])
        for row in replacement_frontier_rows
    )
    replacement_frontier_probe_selected_rows = [
        row for row in replacement_frontier_rows if bool(row["replacement_probe_selected"])
    ]
    replacement_frontier_probe_selected_native_long_range_delta_sum = sum(
        int(row["replacement_native_long_range_delta_after_scoring"])
        for row in replacement_frontier_probe_selected_rows
    )
    replacement_frontier_external_count_delta_positive_rows = [
        row
        for row in replacement_frontier_rows
        if int(row["replacement_external_constraint_coverage_delta"]) > 0
    ]
    replacement_frontier_external_confidence_delta_positive_rows = [
        row
        for row in replacement_frontier_rows
        if float(row["replacement_external_constraint_confidence_delta"]) > 0.0
    ]
    replacement_frontier_external_count_delta_sum = sum(
        int(row["replacement_external_constraint_coverage_delta"])
        for row in replacement_frontier_rows
    )
    replacement_frontier_external_confidence_delta_sum = _rounded(
        sum(
            float(row["replacement_external_constraint_confidence_delta"])
            for row in replacement_frontier_rows
        )
    )
    replacement_frontier_native_long_range_delta_sum_with_external_count_gain = sum(
        int(row["replacement_native_long_range_delta_after_scoring"])
        for row in replacement_frontier_external_count_delta_positive_rows
    )
    replacement_frontier_native_long_range_delta_sum_with_external_confidence_gain = sum(
        int(row["replacement_native_long_range_delta_after_scoring"])
        for row in replacement_frontier_external_confidence_delta_positive_rows
    )
    terminal_bridge_expanded_max_matched_control_precision = _max_selected_metric(
        terminal_bridge_expanded_controls,
        "contact_cluster_precision",
    )
    terminal_bridge_expanded_max_matched_control_long_range_recall = (
        _max_selected_metric(
            terminal_bridge_expanded_controls,
            "long_range_contact_recall",
        )
    )
    terminal_bridge_expanded_max_adversarial_precision = _max_selected_metric(
        adversarial_terminal_bridge_expanded_controls,
        "contact_cluster_precision",
    )
    terminal_bridge_expanded_max_adversarial_long_range_recall = (
        _max_selected_metric(
            adversarial_terminal_bridge_expanded_controls,
            "long_range_contact_recall",
        )
    )
    terminal_bridge_expanded_noninferior_false_rate_vs_matched_controls = (
        bool(terminal_bridge_expanded_controls)
        and terminal_bridge_expanded_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or terminal_bridge_expanded_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in terminal_bridge_expanded_controls
        )
    )
    terminal_bridge_expanded_noninferior_false_rate_vs_adversarial_controls = (
        bool(adversarial_terminal_bridge_expanded_controls)
        and terminal_bridge_expanded_selected_event_count > 0
        and all(
            run.metric.selected_event_count == 0
            or terminal_bridge_expanded_metric.false_nucleus_rate
            <= run.metric.false_nucleus_rate
            for run in adversarial_terminal_bridge_expanded_controls
        )
    )
    terminal_bridge_expanded_beats_matched_controls = (
        terminal_bridge_expanded_noninferior_false_rate_vs_matched_controls
        and terminal_bridge_expanded_metric.contact_cluster_precision
        > terminal_bridge_expanded_max_matched_control_precision
        and terminal_bridge_expanded_metric.long_range_contact_recall
        > terminal_bridge_expanded_max_matched_control_long_range_recall
    )
    terminal_bridge_expanded_beats_adversarial_calibrated_controls = (
        terminal_bridge_expanded_noninferior_false_rate_vs_adversarial_controls
        and terminal_bridge_expanded_metric.contact_cluster_precision
        > terminal_bridge_expanded_max_adversarial_precision
        and terminal_bridge_expanded_metric.long_range_contact_recall
        > terminal_bridge_expanded_max_adversarial_long_range_recall
    )
    return {
        "report_kind": EXTERNAL_COUPLING_TRACE_LOOP_REPORT_KIND,
        "batch_id": EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID,
        "source_benchmark_file": str(source_benchmark_file),
        "external_coupling_file": str(external_coupling_file),
        "oracle_positive_control_file": str(oracle_coupling_file),
        "benchmark_size": len(rows),
        "result": result,
        "reason": reason,
        "external_probe_passed": external_probe_passed,
        "external_real_probe_passed": acceptance_criteria_met,
        "external_rank_consistent_cluster_gated_probe_passed": (
            rank_consistent_acceptance_criteria_met
        ),
        "hard_adversarial_calibrated_probe_passed": (
            hard_adversarial_calibrated_probe_passed
        ),
        "external_couplings_available_rows": available_rows,
        "external_rows_rejected_low_depth": rejected_low_depth,
        "external_rows_rejected_low_coverage": rejected_low_coverage,
        "external_rows_rejected_mapping": rejected_mapping,
        "external_rows_rejected_coordinate_taint": rejected_coordinate_taint,
        "usable_external_rows": available_rows,
        "external_constraint_count": len(import_result.dataset.constraints),
        "external_real_false_nucleus_rate": (
            external_real.metric.false_nucleus_rate if external_metric_defined else None
        ),
        "physical_rerank_false_nucleus_rate": (
            physical_baseline.metric.false_nucleus_rate
        ),
        "external_real_cluster_precision": (
            external_real.metric.contact_cluster_precision
            if external_metric_defined
            else None
        ),
        "physical_rerank_cluster_precision": (
            physical_baseline.metric.contact_cluster_precision
        ),
        "external_real_long_range_recall": (
            external_real.metric.long_range_contact_recall
            if external_metric_defined
            else None
        ),
        "oracle_trace_loop_long_range_recall": (
            oracle_positive_control.metric.long_range_contact_recall
        ),
        "oracle_trace_loop_long_range_recall_floor_50pct": oracle_recall_floor,
        "external_real_vs_control_enrichment_ratio": (
            external_real_vs_control_enrichment_ratio
        ),
        "external_real_mean_selected_coupling_selectivity_score": (
            external_real.metric.mean_selected_coupling_selectivity_score
        ),
        "max_matched_control_mean_selected_coupling_selectivity_score": (
            _max_metric(
                matched_controls,
                "mean_selected_coupling_selectivity_score",
            )
        ),
        "external_real_mean_decoy_coupling_selectivity_score": (
            external_real.metric.mean_decoy_coupling_selectivity_score
        ),
        "max_matched_control_mean_decoy_coupling_selectivity_score": (
            _max_metric(
                matched_controls,
                "mean_decoy_coupling_selectivity_score",
            )
        ),
        "external_real_mean_coupling_decoy_selectivity_margin": (
            external_real.metric.mean_coupling_decoy_selectivity_margin
        ),
        "max_matched_control_mean_coupling_decoy_selectivity_margin": (
            _max_metric(
                matched_controls,
                "mean_coupling_decoy_selectivity_margin",
            )
        ),
        "mean_matched_control_false_nucleus_rate": _mean(control_false_rates),
        "mean_matched_control_cluster_precision": _mean(control_precisions),
        "external_real_beats_physical": external_real_beats_physical,
        "external_real_beats_matched_controls": (
            external_real_beats_matched_controls
        ),
        "external_real_meets_oracle_recall_floor": (
            external_real_meets_oracle_recall_floor
        ),
        "external_margin_gated_selected_event_count": margin_selected_event_count,
        "external_margin_gated_false_nucleus_rate": (
            external_margin_gated.metric.false_nucleus_rate
            if margin_selected_event_count
            else None
        ),
        "external_margin_gated_cluster_precision": (
            external_margin_gated.metric.contact_cluster_precision
            if margin_selected_event_count
            else None
        ),
        "external_margin_gated_long_range_recall": (
            external_margin_gated.metric.long_range_contact_recall
            if margin_selected_event_count
            else None
        ),
        "external_margin_gated_vs_control_enrichment_ratio": (
            margin_vs_control_enrichment_ratio
        ),
        "external_margin_gated_beats_physical": margin_beats_physical,
        "external_margin_gated_beats_matched_controls": (
            margin_beats_matched_controls
        ),
        "external_margin_gated_meets_oracle_recall_floor": (
            margin_meets_oracle_recall_floor
        ),
        "external_margin_gated_claim_allowed": False,
        "external_top_rank_gated_selected_event_count": (
            top_rank_selected_event_count
        ),
        "external_top_rank_gated_false_nucleus_rate": (
            external_top_rank_gated.metric.false_nucleus_rate
            if top_rank_selected_event_count
            else None
        ),
        "external_top_rank_gated_cluster_precision": (
            external_top_rank_gated.metric.contact_cluster_precision
            if top_rank_selected_event_count
            else None
        ),
        "external_top_rank_gated_long_range_recall": (
            external_top_rank_gated.metric.long_range_contact_recall
            if top_rank_selected_event_count
            else None
        ),
        "external_top_rank_gated_vs_control_enrichment_ratio": (
            top_rank_vs_control_enrichment_ratio
        ),
        "external_top_rank_gated_beats_physical": top_rank_beats_physical,
        "external_top_rank_gated_beats_matched_controls": (
            top_rank_beats_matched_controls
        ),
        "external_top_rank_gated_meets_oracle_recall_floor": (
            top_rank_meets_oracle_recall_floor
        ),
        "external_top_rank_gated_claim_allowed": False,
        "external_core_expanded_selected_event_count": (
            core_expanded_selected_event_count
        ),
        "external_core_expanded_false_nucleus_rate": (
            external_core_expanded.metric.false_nucleus_rate
            if core_expanded_selected_event_count
            else None
        ),
        "external_core_expanded_cluster_precision": (
            external_core_expanded.metric.contact_cluster_precision
            if core_expanded_selected_event_count
            else None
        ),
        "external_core_expanded_long_range_recall": (
            external_core_expanded.metric.long_range_contact_recall
            if core_expanded_selected_event_count
            else None
        ),
        "external_core_expanded_vs_control_enrichment_ratio": (
            core_expanded_vs_control_enrichment_ratio
        ),
        "external_core_expanded_beats_physical": core_expanded_beats_physical,
        "external_core_expanded_beats_matched_controls": (
            core_expanded_beats_matched_controls
        ),
        "external_core_expanded_meets_oracle_recall_floor": (
            core_expanded_meets_oracle_recall_floor
        ),
        "external_core_expanded_claim_allowed": False,
        "external_cluster_gated_core_expanded_selected_event_count": (
            cluster_gated_selected_event_count
        ),
        "external_cluster_gated_core_expanded_false_nucleus_rate": (
            external_cluster_gated_core_expanded.metric.false_nucleus_rate
            if cluster_gated_selected_event_count
            else None
        ),
        "external_cluster_gated_core_expanded_cluster_precision": (
            external_cluster_gated_core_expanded.metric.contact_cluster_precision
            if cluster_gated_selected_event_count
            else None
        ),
        "external_cluster_gated_core_expanded_long_range_recall": (
            external_cluster_gated_core_expanded.metric.long_range_contact_recall
            if cluster_gated_selected_event_count
            else None
        ),
        "external_cluster_gated_core_expanded_vs_control_enrichment_ratio": (
            cluster_gated_vs_control_enrichment_ratio
        ),
        "external_cluster_gated_core_expanded_mean_selected_coupling_selectivity_score": (
            cluster_metric.mean_selected_coupling_selectivity_score
        ),
        "external_cluster_gated_core_expanded_max_control_mean_selected_coupling_selectivity_score": (
            _max_metric(
                cluster_gated_core_expanded_controls,
                "mean_selected_coupling_selectivity_score",
            )
        ),
        "external_cluster_gated_core_expanded_mean_decoy_coupling_selectivity_score": (
            cluster_metric.mean_decoy_coupling_selectivity_score
        ),
        "external_cluster_gated_core_expanded_max_control_mean_decoy_coupling_selectivity_score": (
            _max_metric(
                cluster_gated_core_expanded_controls,
                "mean_decoy_coupling_selectivity_score",
            )
        ),
        "external_cluster_gated_core_expanded_mean_coupling_decoy_selectivity_margin": (
            cluster_metric.mean_coupling_decoy_selectivity_margin
        ),
        "external_cluster_gated_core_expanded_max_control_mean_coupling_decoy_selectivity_margin": (
            _max_metric(
                cluster_gated_core_expanded_controls,
                "mean_coupling_decoy_selectivity_margin",
            )
        ),
        "external_cluster_gated_core_expanded_beats_physical": (
            cluster_gated_beats_physical
        ),
        "external_cluster_gated_core_expanded_beats_matched_controls": (
            cluster_gated_beats_matched_controls
        ),
        "external_cluster_gated_core_expanded_meets_oracle_recall_floor": (
            cluster_gated_meets_oracle_recall_floor
        ),
        "external_cluster_gated_core_expanded_claim_allowed": False,
        "external_rank_consistent_cluster_gated_selected_event_count": (
            rank_consistent_selected_event_count
        ),
        "external_rank_consistent_cluster_gated_false_nucleus_rate": (
            external_rank_consistent_cluster_gated.metric.false_nucleus_rate
            if rank_consistent_selected_event_count
            else None
        ),
        "external_rank_consistent_cluster_gated_cluster_precision": (
            external_rank_consistent_cluster_gated.metric.contact_cluster_precision
            if rank_consistent_selected_event_count
            else None
        ),
        "external_rank_consistent_cluster_gated_long_range_recall": (
            external_rank_consistent_cluster_gated.metric.long_range_contact_recall
            if rank_consistent_selected_event_count
            else None
        ),
        "external_rank_consistent_cluster_gated_vs_control_enrichment_ratio": (
            rank_consistent_vs_control_enrichment_ratio
        ),
        "external_rank_consistent_cluster_gated_vs_adversarial_calibrated_enrichment_ratio": (
            rank_consistent_vs_adversarial_enrichment_ratio
        ),
        "adversarial_calibrated_enrichment_min_selected_events": (
            ADVERSARIAL_ENRICHMENT_MIN_SELECTED_EVENTS
        ),
        "external_rank_consistent_cluster_gated_mean_selected_coupling_selectivity_score": (
            external_rank_consistent_cluster_gated.metric.mean_selected_coupling_selectivity_score
        ),
        "external_rank_consistent_cluster_gated_max_control_mean_selected_coupling_selectivity_score": (
            _max_metric(
                rank_consistent_cluster_gated_controls,
                "mean_selected_coupling_selectivity_score",
            )
        ),
        "external_rank_consistent_cluster_gated_mean_coupling_decoy_selectivity_margin": (
            external_rank_consistent_cluster_gated.metric.mean_coupling_decoy_selectivity_margin
        ),
        "external_rank_consistent_cluster_gated_max_control_mean_coupling_decoy_selectivity_margin": (
            _max_metric(
                rank_consistent_cluster_gated_controls,
                "mean_coupling_decoy_selectivity_margin",
            )
        ),
        "external_rank_consistent_cluster_gated_beats_physical": (
            rank_consistent_beats_physical
        ),
        "external_rank_consistent_cluster_gated_beats_matched_controls": (
            rank_consistent_beats_matched_controls
        ),
        "external_rank_consistent_cluster_gated_beats_adversarial_calibrated_controls": (
            rank_consistent_beats_adversarial_calibrated_controls
        ),
        "external_rank_consistent_cluster_gated_meets_oracle_recall_floor": (
            rank_consistent_meets_oracle_recall_floor
        ),
        "external_rank_consistent_cluster_gated_claim_allowed": False,
        "external_persistent_rank_consistent_cluster_gated_selected_event_count": (
            persistent_selected_event_count
        ),
        "external_persistent_rank_consistent_cluster_gated_false_nucleus_rate": (
            external_persistent_rank_consistent_cluster_gated.metric.false_nucleus_rate
            if persistent_selected_event_count
            else None
        ),
        "external_persistent_rank_consistent_cluster_gated_cluster_precision": (
            external_persistent_rank_consistent_cluster_gated.metric.contact_cluster_precision
            if persistent_selected_event_count
            else None
        ),
        "external_persistent_rank_consistent_cluster_gated_long_range_recall": (
            external_persistent_rank_consistent_cluster_gated.metric.long_range_contact_recall
            if persistent_selected_event_count
            else None
        ),
        "external_persistent_rank_consistent_cluster_gated_vs_control_enrichment_ratio": (
            persistent_vs_control_enrichment_ratio
        ),
        "external_persistent_rank_consistent_cluster_gated_vs_adversarial_calibrated_enrichment_ratio": (
            persistent_vs_adversarial_enrichment_ratio
        ),
        "external_persistent_rank_consistent_cluster_gated_real_vs_decoy_coupling_nucleus_enrichment_ratio": (
            external_persistent_rank_consistent_cluster_gated.metric.real_vs_decoy_coupling_nucleus_enrichment_ratio
        ),
        "external_persistent_rank_consistent_cluster_gated_vs_control_nucleus_score_enrichment_ratio": (
            persistent_vs_control_nucleus_score_enrichment_ratio
        ),
        "external_persistent_rank_consistent_cluster_gated_vs_adversarial_calibrated_nucleus_score_enrichment_ratio": (
            persistent_vs_adversarial_nucleus_score_enrichment_ratio
        ),
        "external_persistent_rank_consistent_cluster_gated_mean_selected_coupling_selectivity_score": (
            external_persistent_rank_consistent_cluster_gated.metric.mean_selected_coupling_selectivity_score
        ),
        "external_persistent_rank_consistent_cluster_gated_max_control_mean_selected_coupling_selectivity_score": (
            _max_metric(
                persistent_rank_consistent_controls,
                "mean_selected_coupling_selectivity_score",
            )
        ),
        "external_persistent_rank_consistent_cluster_gated_mean_coupling_decoy_selectivity_margin": (
            external_persistent_rank_consistent_cluster_gated.metric.mean_coupling_decoy_selectivity_margin
        ),
        "external_persistent_rank_consistent_cluster_gated_max_control_mean_coupling_decoy_selectivity_margin": (
            _max_metric(
                persistent_rank_consistent_controls,
                "mean_coupling_decoy_selectivity_margin",
            )
        ),
        "external_persistent_rank_consistent_cluster_gated_mean_coupling_nucleus_score": (
            external_persistent_rank_consistent_cluster_gated.metric.mean_coupling_nucleus_score
        ),
        "external_persistent_rank_consistent_cluster_gated_mean_decoy_coupling_nucleus_score": (
            external_persistent_rank_consistent_cluster_gated.metric.mean_decoy_coupling_nucleus_score
        ),
        "external_persistent_rank_consistent_cluster_gated_mean_coupling_nucleus_decoy_margin": (
            external_persistent_rank_consistent_cluster_gated.metric.mean_coupling_nucleus_decoy_margin
        ),
        "external_persistent_rank_consistent_cluster_gated_max_control_coupling_nucleus_enrichment_ratio": (
            max_persistent_control_nucleus_enrichment
        ),
        "external_persistent_rank_consistent_cluster_gated_max_adversarial_calibrated_coupling_nucleus_enrichment_ratio": (
            max_adversarial_persistent_nucleus_enrichment
        ),
        "external_persistent_rank_consistent_cluster_gated_beats_physical": (
            persistent_beats_physical
        ),
        "external_persistent_rank_consistent_cluster_gated_beats_matched_controls": (
            persistent_beats_matched_controls
        ),
        "external_persistent_rank_consistent_cluster_gated_beats_adversarial_calibrated_controls": (
            persistent_beats_adversarial_calibrated_controls
        ),
        "external_persistent_rank_consistent_cluster_gated_meets_oracle_recall_floor": (
            persistent_meets_oracle_recall_floor
        ),
        "external_persistent_rank_consistent_cluster_gated_probe_passed": (
            persistent_acceptance_criteria_met
        ),
        "external_persistent_rank_consistent_cluster_gated_hard_adversarial_calibrated_probe_passed": (
            persistent_hard_adversarial_calibrated_probe_passed
        ),
        "external_persistent_rank_consistent_cluster_gated_selector_score_probe_passed": (
            persistent_selector_score_acceptance_criteria_met
        ),
        "external_persistent_rank_consistent_cluster_gated_hard_selector_score_probe_passed": (
            persistent_hard_selector_score_probe_passed
        ),
        "external_persistent_rank_consistent_cluster_gated_claim_allowed": False,
        "external_score_margin_expanded_selected_event_count": (
            score_margin_expanded_selected_event_count
        ),
        "external_score_margin_expanded_added_event_count": (
            len(score_margin_expanded_added_events)
        ),
        "external_score_margin_expanded_added_native_contact_count": (
            sum(
                event.native_contact_count_after_scoring
                for event in score_margin_expanded_added_events
            )
        ),
        "external_score_margin_expanded_added_native_long_range_contact_count": (
            sum(
                event.native_long_range_contact_count_after_scoring
                for event in score_margin_expanded_added_events
            )
        ),
        "external_score_margin_expanded_added_false_event_count": (
            sum(
                1
                for event in score_margin_expanded_added_events
                if event.native_contact_count_after_scoring == 0
            )
        ),
        "external_score_margin_expanded_false_nucleus_rate": (
            score_margin_expanded_metric.false_nucleus_rate
            if score_margin_expanded_selected_event_count
            else None
        ),
        "external_score_margin_expanded_cluster_precision": (
            score_margin_expanded_metric.contact_cluster_precision
            if score_margin_expanded_selected_event_count
            else None
        ),
        "external_score_margin_expanded_long_range_recall": (
            score_margin_expanded_metric.long_range_contact_recall
            if score_margin_expanded_selected_event_count
            else None
        ),
        "external_score_margin_expanded_long_range_recall_delta_vs_persistent": (
            _rounded(
                score_margin_expanded_metric.long_range_contact_recall
                - external_persistent_rank_consistent_cluster_gated.metric.long_range_contact_recall
            )
        ),
        "external_score_margin_expanded_max_matched_control_cluster_precision": (
            score_margin_expanded_max_matched_control_precision
        ),
        "external_score_margin_expanded_max_matched_control_long_range_recall": (
            score_margin_expanded_max_matched_control_long_range_recall
        ),
        "external_score_margin_expanded_cluster_precision_margin_vs_matched_controls": (
            _rounded(
                score_margin_expanded_metric.contact_cluster_precision
                - score_margin_expanded_max_matched_control_precision
            )
        ),
        "external_score_margin_expanded_long_range_recall_margin_vs_matched_controls": (
            _rounded(
                score_margin_expanded_metric.long_range_contact_recall
                - score_margin_expanded_max_matched_control_long_range_recall
            )
        ),
        "external_score_margin_expanded_max_adversarial_cluster_precision": (
            score_margin_expanded_max_adversarial_precision
        ),
        "external_score_margin_expanded_max_adversarial_long_range_recall": (
            score_margin_expanded_max_adversarial_long_range_recall
        ),
        "external_score_margin_expanded_cluster_precision_margin_vs_adversarial_controls": (
            _rounded(
                score_margin_expanded_metric.contact_cluster_precision
                - score_margin_expanded_max_adversarial_precision
            )
        ),
        "external_score_margin_expanded_long_range_recall_margin_vs_adversarial_controls": (
            _rounded(
                score_margin_expanded_metric.long_range_contact_recall
                - score_margin_expanded_max_adversarial_long_range_recall
            )
        ),
        "external_score_margin_expanded_noninferior_false_rate_vs_matched_controls": (
            score_margin_expanded_noninferior_false_rate_vs_matched_controls
        ),
        "external_score_margin_expanded_noninferior_false_rate_vs_adversarial_controls": (
            score_margin_expanded_noninferior_false_rate_vs_adversarial_controls
        ),
        "external_score_margin_expanded_beats_matched_controls": (
            score_margin_expanded_beats_matched_controls
        ),
        "external_score_margin_expanded_beats_adversarial_calibrated_controls": (
            score_margin_expanded_beats_adversarial_calibrated_controls
        ),
        "external_score_margin_expanded_claim_allowed": False,
        "external_boundary_continuity_expanded_selected_event_count": (
            boundary_continuity_expanded_selected_event_count
        ),
        "external_boundary_continuity_expanded_added_event_count": (
            len(boundary_continuity_expanded_added_events)
        ),
        "external_boundary_continuity_expanded_added_native_contact_count": (
            sum(
                event.native_contact_count_after_scoring
                for event in boundary_continuity_expanded_added_events
            )
        ),
        "external_boundary_continuity_expanded_added_native_long_range_contact_count": (
            sum(
                event.native_long_range_contact_count_after_scoring
                for event in boundary_continuity_expanded_added_events
            )
        ),
        "external_boundary_continuity_expanded_added_false_event_count": (
            sum(
                1
                for event in boundary_continuity_expanded_added_events
                if event.native_contact_count_after_scoring == 0
            )
        ),
        "external_boundary_continuity_expanded_false_nucleus_rate": (
            boundary_continuity_expanded_metric.false_nucleus_rate
            if boundary_continuity_expanded_selected_event_count
            else None
        ),
        "external_boundary_continuity_expanded_cluster_precision": (
            boundary_continuity_expanded_metric.contact_cluster_precision
            if boundary_continuity_expanded_selected_event_count
            else None
        ),
        "external_boundary_continuity_expanded_long_range_recall": (
            boundary_continuity_expanded_metric.long_range_contact_recall
            if boundary_continuity_expanded_selected_event_count
            else None
        ),
        "external_boundary_continuity_expanded_long_range_recall_delta_vs_score_margin": (
            _rounded(
                boundary_continuity_expanded_metric.long_range_contact_recall
                - score_margin_expanded_metric.long_range_contact_recall
            )
        ),
        "external_boundary_continuity_expanded_max_matched_control_cluster_precision": (
            boundary_continuity_expanded_max_matched_control_precision
        ),
        "external_boundary_continuity_expanded_max_matched_control_long_range_recall": (
            boundary_continuity_expanded_max_matched_control_long_range_recall
        ),
        "external_boundary_continuity_expanded_cluster_precision_margin_vs_matched_controls": (
            _rounded(
                boundary_continuity_expanded_metric.contact_cluster_precision
                - boundary_continuity_expanded_max_matched_control_precision
            )
        ),
        "external_boundary_continuity_expanded_long_range_recall_margin_vs_matched_controls": (
            _rounded(
                boundary_continuity_expanded_metric.long_range_contact_recall
                - boundary_continuity_expanded_max_matched_control_long_range_recall
            )
        ),
        "external_boundary_continuity_expanded_max_adversarial_cluster_precision": (
            boundary_continuity_expanded_max_adversarial_precision
        ),
        "external_boundary_continuity_expanded_max_adversarial_long_range_recall": (
            boundary_continuity_expanded_max_adversarial_long_range_recall
        ),
        "external_boundary_continuity_expanded_cluster_precision_margin_vs_adversarial_controls": (
            _rounded(
                boundary_continuity_expanded_metric.contact_cluster_precision
                - boundary_continuity_expanded_max_adversarial_precision
            )
        ),
        "external_boundary_continuity_expanded_long_range_recall_margin_vs_adversarial_controls": (
            _rounded(
                boundary_continuity_expanded_metric.long_range_contact_recall
                - boundary_continuity_expanded_max_adversarial_long_range_recall
            )
        ),
        "external_boundary_continuity_expanded_noninferior_false_rate_vs_matched_controls": (
            boundary_continuity_expanded_noninferior_false_rate_vs_matched_controls
        ),
        "external_boundary_continuity_expanded_noninferior_false_rate_vs_adversarial_controls": (
            boundary_continuity_expanded_noninferior_false_rate_vs_adversarial_controls
        ),
        "external_boundary_continuity_expanded_beats_matched_controls": (
            boundary_continuity_expanded_beats_matched_controls
        ),
        "external_boundary_continuity_expanded_beats_adversarial_calibrated_controls": (
            boundary_continuity_expanded_beats_adversarial_calibrated_controls
        ),
        "external_boundary_continuity_expanded_claim_allowed": False,
        "external_edge_continuity_expanded_selected_event_count": (
            edge_continuity_expanded_selected_event_count
        ),
        "external_edge_continuity_expanded_added_event_count": (
            len(edge_continuity_expanded_added_events)
        ),
        "external_edge_continuity_expanded_added_native_contact_count": (
            sum(
                event.native_contact_count_after_scoring
                for event in edge_continuity_expanded_added_events
            )
        ),
        "external_edge_continuity_expanded_added_native_long_range_contact_count": (
            sum(
                event.native_long_range_contact_count_after_scoring
                for event in edge_continuity_expanded_added_events
            )
        ),
        "external_edge_continuity_expanded_added_false_event_count": (
            sum(
                1
                for event in edge_continuity_expanded_added_events
                if event.native_contact_count_after_scoring == 0
            )
        ),
        "external_edge_continuity_expanded_false_nucleus_rate": (
            edge_continuity_expanded_metric.false_nucleus_rate
            if edge_continuity_expanded_selected_event_count
            else None
        ),
        "external_edge_continuity_expanded_cluster_precision": (
            edge_continuity_expanded_metric.contact_cluster_precision
            if edge_continuity_expanded_selected_event_count
            else None
        ),
        "external_edge_continuity_expanded_long_range_recall": (
            edge_continuity_expanded_metric.long_range_contact_recall
            if edge_continuity_expanded_selected_event_count
            else None
        ),
        "external_edge_continuity_expanded_long_range_recall_delta_vs_boundary_continuity": (
            _rounded(
                edge_continuity_expanded_metric.long_range_contact_recall
                - boundary_continuity_expanded_metric.long_range_contact_recall
            )
        ),
        "external_edge_continuity_expanded_max_matched_control_cluster_precision": (
            edge_continuity_expanded_max_matched_control_precision
        ),
        "external_edge_continuity_expanded_max_matched_control_long_range_recall": (
            edge_continuity_expanded_max_matched_control_long_range_recall
        ),
        "external_edge_continuity_expanded_cluster_precision_margin_vs_matched_controls": (
            _rounded(
                edge_continuity_expanded_metric.contact_cluster_precision
                - edge_continuity_expanded_max_matched_control_precision
            )
        ),
        "external_edge_continuity_expanded_long_range_recall_margin_vs_matched_controls": (
            _rounded(
                edge_continuity_expanded_metric.long_range_contact_recall
                - edge_continuity_expanded_max_matched_control_long_range_recall
            )
        ),
        "external_edge_continuity_expanded_max_adversarial_cluster_precision": (
            edge_continuity_expanded_max_adversarial_precision
        ),
        "external_edge_continuity_expanded_max_adversarial_long_range_recall": (
            edge_continuity_expanded_max_adversarial_long_range_recall
        ),
        "external_edge_continuity_expanded_cluster_precision_margin_vs_adversarial_controls": (
            _rounded(
                edge_continuity_expanded_metric.contact_cluster_precision
                - edge_continuity_expanded_max_adversarial_precision
            )
        ),
        "external_edge_continuity_expanded_long_range_recall_margin_vs_adversarial_controls": (
            _rounded(
                edge_continuity_expanded_metric.long_range_contact_recall
                - edge_continuity_expanded_max_adversarial_long_range_recall
            )
        ),
        "external_edge_continuity_expanded_noninferior_false_rate_vs_matched_controls": (
            edge_continuity_expanded_noninferior_false_rate_vs_matched_controls
        ),
        "external_edge_continuity_expanded_noninferior_false_rate_vs_adversarial_controls": (
            edge_continuity_expanded_noninferior_false_rate_vs_adversarial_controls
        ),
        "external_edge_continuity_expanded_beats_matched_controls": (
            edge_continuity_expanded_beats_matched_controls
        ),
        "external_edge_continuity_expanded_beats_adversarial_calibrated_controls": (
            edge_continuity_expanded_beats_adversarial_calibrated_controls
        ),
        "external_edge_continuity_expanded_claim_allowed": False,
        "external_pressure_release_expanded_selected_event_count": (
            pressure_release_expanded_selected_event_count
        ),
        "external_pressure_release_expanded_added_event_count": (
            len(pressure_release_expanded_added_events)
        ),
        "external_pressure_release_expanded_added_native_contact_count": (
            sum(
                event.native_contact_count_after_scoring
                for event in pressure_release_expanded_added_events
            )
        ),
        "external_pressure_release_expanded_added_native_long_range_contact_count": (
            sum(
                event.native_long_range_contact_count_after_scoring
                for event in pressure_release_expanded_added_events
            )
        ),
        "external_pressure_release_expanded_added_false_event_count": (
            sum(
                1
                for event in pressure_release_expanded_added_events
                if event.native_contact_count_after_scoring == 0
            )
        ),
        "external_pressure_release_expanded_false_nucleus_rate": (
            pressure_release_expanded_metric.false_nucleus_rate
            if pressure_release_expanded_selected_event_count
            else None
        ),
        "external_pressure_release_expanded_cluster_precision": (
            pressure_release_expanded_metric.contact_cluster_precision
            if pressure_release_expanded_selected_event_count
            else None
        ),
        "external_pressure_release_expanded_long_range_recall": (
            pressure_release_expanded_metric.long_range_contact_recall
            if pressure_release_expanded_selected_event_count
            else None
        ),
        "external_pressure_release_expanded_long_range_recall_delta_vs_edge_continuity": (
            _rounded(
                pressure_release_expanded_metric.long_range_contact_recall
                - edge_continuity_expanded_metric.long_range_contact_recall
            )
        ),
        "external_pressure_release_expanded_max_matched_control_cluster_precision": (
            pressure_release_expanded_max_matched_control_precision
        ),
        "external_pressure_release_expanded_max_matched_control_long_range_recall": (
            pressure_release_expanded_max_matched_control_long_range_recall
        ),
        "external_pressure_release_expanded_cluster_precision_margin_vs_matched_controls": (
            _rounded(
                pressure_release_expanded_metric.contact_cluster_precision
                - pressure_release_expanded_max_matched_control_precision
            )
        ),
        "external_pressure_release_expanded_long_range_recall_margin_vs_matched_controls": (
            _rounded(
                pressure_release_expanded_metric.long_range_contact_recall
                - pressure_release_expanded_max_matched_control_long_range_recall
            )
        ),
        "external_pressure_release_expanded_max_adversarial_cluster_precision": (
            pressure_release_expanded_max_adversarial_precision
        ),
        "external_pressure_release_expanded_max_adversarial_long_range_recall": (
            pressure_release_expanded_max_adversarial_long_range_recall
        ),
        "external_pressure_release_expanded_cluster_precision_margin_vs_adversarial_controls": (
            _rounded(
                pressure_release_expanded_metric.contact_cluster_precision
                - pressure_release_expanded_max_adversarial_precision
            )
        ),
        "external_pressure_release_expanded_long_range_recall_margin_vs_adversarial_controls": (
            _rounded(
                pressure_release_expanded_metric.long_range_contact_recall
                - pressure_release_expanded_max_adversarial_long_range_recall
            )
        ),
        "external_pressure_release_expanded_noninferior_false_rate_vs_matched_controls": (
            pressure_release_expanded_noninferior_false_rate_vs_matched_controls
        ),
        "external_pressure_release_expanded_noninferior_false_rate_vs_adversarial_controls": (
            pressure_release_expanded_noninferior_false_rate_vs_adversarial_controls
        ),
        "external_pressure_release_expanded_beats_matched_controls": (
            pressure_release_expanded_beats_matched_controls
        ),
        "external_pressure_release_expanded_beats_adversarial_calibrated_controls": (
            pressure_release_expanded_beats_adversarial_calibrated_controls
        ),
        "external_pressure_release_expanded_claim_allowed": False,
        "external_registry_extension_expanded_selected_event_count": (
            registry_extension_expanded_selected_event_count
        ),
        "external_registry_extension_expanded_added_event_count": (
            len(registry_extension_expanded_added_events)
        ),
        "external_registry_extension_expanded_added_native_contact_count": (
            sum(
                event.native_contact_count_after_scoring
                for event in registry_extension_expanded_added_events
            )
        ),
        "external_registry_extension_expanded_added_native_long_range_contact_count": (
            sum(
                event.native_long_range_contact_count_after_scoring
                for event in registry_extension_expanded_added_events
            )
        ),
        "external_registry_extension_expanded_added_false_event_count": (
            sum(
                1
                for event in registry_extension_expanded_added_events
                if event.native_contact_count_after_scoring == 0
            )
        ),
        "external_registry_extension_expanded_false_nucleus_rate": (
            registry_extension_expanded_metric.false_nucleus_rate
            if registry_extension_expanded_selected_event_count
            else None
        ),
        "external_registry_extension_expanded_cluster_precision": (
            registry_extension_expanded_metric.contact_cluster_precision
            if registry_extension_expanded_selected_event_count
            else None
        ),
        "external_registry_extension_expanded_long_range_recall": (
            registry_extension_expanded_metric.long_range_contact_recall
            if registry_extension_expanded_selected_event_count
            else None
        ),
        "external_registry_extension_expanded_long_range_recall_delta_vs_pressure_release": (
            _rounded(
                registry_extension_expanded_metric.long_range_contact_recall
                - pressure_release_expanded_metric.long_range_contact_recall
            )
        ),
        "external_registry_extension_expanded_max_matched_control_cluster_precision": (
            registry_extension_expanded_max_matched_control_precision
        ),
        "external_registry_extension_expanded_max_matched_control_long_range_recall": (
            registry_extension_expanded_max_matched_control_long_range_recall
        ),
        "external_registry_extension_expanded_cluster_precision_margin_vs_matched_controls": (
            _rounded(
                registry_extension_expanded_metric.contact_cluster_precision
                - registry_extension_expanded_max_matched_control_precision
            )
        ),
        "external_registry_extension_expanded_long_range_recall_margin_vs_matched_controls": (
            _rounded(
                registry_extension_expanded_metric.long_range_contact_recall
                - registry_extension_expanded_max_matched_control_long_range_recall
            )
        ),
        "external_registry_extension_expanded_max_adversarial_cluster_precision": (
            registry_extension_expanded_max_adversarial_precision
        ),
        "external_registry_extension_expanded_max_adversarial_long_range_recall": (
            registry_extension_expanded_max_adversarial_long_range_recall
        ),
        "external_registry_extension_expanded_cluster_precision_margin_vs_adversarial_controls": (
            _rounded(
                registry_extension_expanded_metric.contact_cluster_precision
                - registry_extension_expanded_max_adversarial_precision
            )
        ),
        "external_registry_extension_expanded_long_range_recall_margin_vs_adversarial_controls": (
            _rounded(
                registry_extension_expanded_metric.long_range_contact_recall
                - registry_extension_expanded_max_adversarial_long_range_recall
            )
        ),
        "external_registry_extension_expanded_noninferior_false_rate_vs_matched_controls": (
            registry_extension_expanded_noninferior_false_rate_vs_matched_controls
        ),
        "external_registry_extension_expanded_noninferior_false_rate_vs_adversarial_controls": (
            registry_extension_expanded_noninferior_false_rate_vs_adversarial_controls
        ),
        "external_registry_extension_expanded_beats_matched_controls": (
            registry_extension_expanded_beats_matched_controls
        ),
        "external_registry_extension_expanded_beats_adversarial_calibrated_controls": (
            registry_extension_expanded_beats_adversarial_calibrated_controls
        ),
        "external_registry_extension_expanded_claim_allowed": False,
        "external_terminal_bridge_expanded_selected_event_count": (
            terminal_bridge_expanded_selected_event_count
        ),
        "external_terminal_bridge_expanded_added_event_count": (
            len(terminal_bridge_expanded_added_events)
        ),
        "external_terminal_bridge_expanded_added_native_contact_count": (
            sum(
                event.native_contact_count_after_scoring
                for event in terminal_bridge_expanded_added_events
            )
        ),
        "external_terminal_bridge_expanded_added_native_long_range_contact_count": (
            sum(
                event.native_long_range_contact_count_after_scoring
                for event in terminal_bridge_expanded_added_events
            )
        ),
        "external_terminal_bridge_expanded_added_false_event_count": (
            sum(
                1
                for event in terminal_bridge_expanded_added_events
                if event.native_contact_count_after_scoring == 0
            )
        ),
        "external_terminal_bridge_expanded_false_nucleus_rate": (
            terminal_bridge_expanded_metric.false_nucleus_rate
            if terminal_bridge_expanded_selected_event_count
            else None
        ),
        "external_terminal_bridge_expanded_cluster_precision": (
            terminal_bridge_expanded_metric.contact_cluster_precision
            if terminal_bridge_expanded_selected_event_count
            else None
        ),
        "external_terminal_bridge_expanded_long_range_recall": (
            terminal_bridge_expanded_metric.long_range_contact_recall
            if terminal_bridge_expanded_selected_event_count
            else None
        ),
        "external_terminal_bridge_expanded_long_range_recall_delta_vs_registry_extension": (
            _rounded(
                terminal_bridge_expanded_metric.long_range_contact_recall
                - registry_extension_expanded_metric.long_range_contact_recall
            )
        ),
        "external_terminal_bridge_expanded_max_matched_control_cluster_precision": (
            terminal_bridge_expanded_max_matched_control_precision
        ),
        "external_terminal_bridge_expanded_max_matched_control_long_range_recall": (
            terminal_bridge_expanded_max_matched_control_long_range_recall
        ),
        "external_terminal_bridge_expanded_cluster_precision_margin_vs_matched_controls": (
            _rounded(
                terminal_bridge_expanded_metric.contact_cluster_precision
                - terminal_bridge_expanded_max_matched_control_precision
            )
        ),
        "external_terminal_bridge_expanded_long_range_recall_margin_vs_matched_controls": (
            _rounded(
                terminal_bridge_expanded_metric.long_range_contact_recall
                - terminal_bridge_expanded_max_matched_control_long_range_recall
            )
        ),
        "external_terminal_bridge_expanded_max_adversarial_cluster_precision": (
            terminal_bridge_expanded_max_adversarial_precision
        ),
        "external_terminal_bridge_expanded_max_adversarial_long_range_recall": (
            terminal_bridge_expanded_max_adversarial_long_range_recall
        ),
        "external_terminal_bridge_expanded_cluster_precision_margin_vs_adversarial_controls": (
            _rounded(
                terminal_bridge_expanded_metric.contact_cluster_precision
                - terminal_bridge_expanded_max_adversarial_precision
            )
        ),
        "external_terminal_bridge_expanded_long_range_recall_margin_vs_adversarial_controls": (
            _rounded(
                terminal_bridge_expanded_metric.long_range_contact_recall
                - terminal_bridge_expanded_max_adversarial_long_range_recall
            )
        ),
        "external_terminal_bridge_expanded_noninferior_false_rate_vs_matched_controls": (
            terminal_bridge_expanded_noninferior_false_rate_vs_matched_controls
        ),
        "external_terminal_bridge_expanded_noninferior_false_rate_vs_adversarial_controls": (
            terminal_bridge_expanded_noninferior_false_rate_vs_adversarial_controls
        ),
        "external_terminal_bridge_expanded_beats_matched_controls": (
            terminal_bridge_expanded_beats_matched_controls
        ),
        "external_terminal_bridge_expanded_beats_adversarial_calibrated_controls": (
            terminal_bridge_expanded_beats_adversarial_calibrated_controls
        ),
        "external_terminal_bridge_expanded_claim_allowed": False,
        "external_boundary_field_replacement_probe_selected_event_count": (
            boundary_field_replacement_probe_metric.selected_event_count
        ),
        "external_boundary_field_replacement_probe_added_event_count": (
            len(boundary_field_replacement_probe_added_events)
        ),
        "external_boundary_field_replacement_probe_added_native_contact_count": (
            sum(
                event.native_contact_count_after_scoring
                for event in boundary_field_replacement_probe_added_events
            )
        ),
        "external_boundary_field_replacement_probe_added_native_long_range_contact_count": (
            sum(
                event.native_long_range_contact_count_after_scoring
                for event in boundary_field_replacement_probe_added_events
            )
        ),
        "external_boundary_field_replacement_probe_added_false_event_count": (
            sum(
                1
                for event in boundary_field_replacement_probe_added_events
                if event.native_contact_count_after_scoring == 0
            )
        ),
        "external_boundary_field_replacement_probe_false_nucleus_rate": (
            boundary_field_replacement_probe_metric.false_nucleus_rate
            if boundary_field_replacement_probe_metric.selected_event_count
            else 0.0
        ),
        "external_boundary_field_replacement_probe_cluster_precision": (
            boundary_field_replacement_probe_metric.contact_cluster_precision
            if boundary_field_replacement_probe_metric.selected_event_count
            else 0.0
        ),
        "external_boundary_field_replacement_probe_long_range_recall": (
            boundary_field_replacement_probe_metric.long_range_contact_recall
            if boundary_field_replacement_probe_metric.selected_event_count
            else 0.0
        ),
        "external_boundary_field_replacement_probe_long_range_recall_delta_vs_terminal_bridge": (
            _rounded(
                boundary_field_replacement_probe_metric.long_range_contact_recall
                - terminal_bridge_expanded_metric.long_range_contact_recall
            )
        ),
        "external_boundary_field_replacement_probe_claim_allowed": False,
        "external_macro_scale_future_preserved_segment_length": (
            MACRO_SCALE_SEGMENT_LENGTH
        ),
        "external_macro_scale_future_preserved_segment_stride": (
            MACRO_SCALE_SEGMENT_STRIDE
        ),
        "external_macro_scale_future_preserved_selected_event_count": (
            macro_scale_future_preserved_metric.selected_event_count
        ),
        "external_macro_scale_future_preserved_false_nucleus_rate": (
            macro_scale_future_preserved_metric.false_nucleus_rate
            if macro_scale_future_preserved_metric.selected_event_count
            else 0.0
        ),
        "external_macro_scale_future_preserved_cluster_precision": (
            macro_scale_future_preserved_metric.contact_cluster_precision
            if macro_scale_future_preserved_metric.selected_event_count
            else 0.0
        ),
        "external_macro_scale_future_preserved_long_range_recall": (
            macro_scale_future_preserved_metric.long_range_contact_recall
            if macro_scale_future_preserved_metric.selected_event_count
            else 0.0
        ),
        "external_macro_scale_future_preserved_long_range_recall_delta_vs_boundary_replacement": (
            _rounded(
                macro_scale_future_preserved_metric.long_range_contact_recall
                - boundary_field_replacement_probe_metric.long_range_contact_recall
            )
        ),
        "external_macro_scale_future_preserved_max_matched_control_cluster_precision": (
            macro_scale_future_preserved_max_matched_control_precision
        ),
        "external_macro_scale_future_preserved_max_matched_control_long_range_recall": (
            macro_scale_future_preserved_max_matched_control_long_range_recall
        ),
        "external_macro_scale_future_preserved_cluster_precision_margin_vs_matched_controls": (
            _rounded(
                macro_scale_future_preserved_metric.contact_cluster_precision
                - macro_scale_future_preserved_max_matched_control_precision
            )
        ),
        "external_macro_scale_future_preserved_long_range_recall_margin_vs_matched_controls": (
            _rounded(
                macro_scale_future_preserved_metric.long_range_contact_recall
                - macro_scale_future_preserved_max_matched_control_long_range_recall
            )
        ),
        "external_macro_scale_future_preserved_max_adversarial_cluster_precision": (
            macro_scale_future_preserved_max_adversarial_precision
        ),
        "external_macro_scale_future_preserved_max_adversarial_long_range_recall": (
            macro_scale_future_preserved_max_adversarial_long_range_recall
        ),
        "external_macro_scale_future_preserved_cluster_precision_margin_vs_adversarial_controls": (
            _rounded(
                macro_scale_future_preserved_metric.contact_cluster_precision
                - macro_scale_future_preserved_max_adversarial_precision
            )
        ),
        "external_macro_scale_future_preserved_long_range_recall_margin_vs_adversarial_controls": (
            _rounded(
                macro_scale_future_preserved_metric.long_range_contact_recall
                - macro_scale_future_preserved_max_adversarial_long_range_recall
            )
        ),
        "external_macro_scale_future_preserved_noninferior_false_rate_vs_matched_controls": (
            macro_scale_future_preserved_noninferior_false_rate_vs_matched_controls
        ),
        "external_macro_scale_future_preserved_noninferior_false_rate_vs_adversarial_controls": (
            macro_scale_future_preserved_noninferior_false_rate_vs_adversarial_controls
        ),
        "external_macro_scale_future_preserved_beats_matched_controls": (
            macro_scale_future_preserved_beats_matched_controls
        ),
        "external_macro_scale_future_preserved_beats_adversarial_calibrated_controls": (
            macro_scale_future_preserved_beats_adversarial_calibrated_controls
        ),
        "external_macro_scale_future_preserved_claim_allowed": False,
        "external_multiscale_future_preserved_segment_lengths": tuple(
            segment_length
            for segment_length, _ in MULTISCALE_FUTURE_PRESERVED_CONFIGS
        ),
        "external_multiscale_future_preserved_future_preservation_mins": tuple(
            future_preservation_min
            for _, future_preservation_min in MULTISCALE_FUTURE_PRESERVED_CONFIGS
        ),
        "external_multiscale_future_preserved_segment_stride": (
            MULTISCALE_FUTURE_PRESERVED_SEGMENT_STRIDE
        ),
        "external_multiscale_future_preserved_max_events_per_row": (
            MULTISCALE_FUTURE_PRESERVED_MAX_EVENTS_PER_ROW
        ),
        "external_multiscale_future_preserved_selected_event_count": (
            multiscale_future_preserved_metric.selected_event_count
        ),
        "external_multiscale_future_preserved_false_nucleus_rate": (
            multiscale_future_preserved_metric.false_nucleus_rate
            if multiscale_future_preserved_selected_event_count
            else 0.0
        ),
        "external_multiscale_future_preserved_cluster_precision": (
            multiscale_future_preserved_metric.contact_cluster_precision
            if multiscale_future_preserved_selected_event_count
            else 0.0
        ),
        "external_multiscale_future_preserved_long_range_recall": (
            multiscale_future_preserved_metric.long_range_contact_recall
            if multiscale_future_preserved_selected_event_count
            else 0.0
        ),
        "external_multiscale_future_preserved_long_range_recall_delta_vs_macro": (
            _rounded(
                multiscale_future_preserved_metric.long_range_contact_recall
                - macro_scale_future_preserved_metric.long_range_contact_recall
            )
        ),
        "external_multiscale_future_preserved_max_matched_control_cluster_precision": (
            multiscale_future_preserved_max_matched_control_precision
        ),
        "external_multiscale_future_preserved_max_matched_control_long_range_recall": (
            multiscale_future_preserved_max_matched_control_long_range_recall
        ),
        "external_multiscale_future_preserved_cluster_precision_margin_vs_matched_controls": (
            _rounded(
                multiscale_future_preserved_metric.contact_cluster_precision
                - multiscale_future_preserved_max_matched_control_precision
            )
        ),
        "external_multiscale_future_preserved_long_range_recall_margin_vs_matched_controls": (
            _rounded(
                multiscale_future_preserved_metric.long_range_contact_recall
                - multiscale_future_preserved_max_matched_control_long_range_recall
            )
        ),
        "external_multiscale_future_preserved_max_adversarial_cluster_precision": (
            multiscale_future_preserved_max_adversarial_precision
        ),
        "external_multiscale_future_preserved_max_adversarial_long_range_recall": (
            multiscale_future_preserved_max_adversarial_long_range_recall
        ),
        "external_multiscale_future_preserved_cluster_precision_margin_vs_adversarial_controls": (
            _rounded(
                multiscale_future_preserved_metric.contact_cluster_precision
                - multiscale_future_preserved_max_adversarial_precision
            )
        ),
        "external_multiscale_future_preserved_long_range_recall_margin_vs_adversarial_controls": (
            _rounded(
                multiscale_future_preserved_metric.long_range_contact_recall
                - multiscale_future_preserved_max_adversarial_long_range_recall
            )
        ),
        "external_multiscale_future_preserved_noninferior_false_rate_vs_matched_controls": (
            multiscale_future_preserved_noninferior_false_rate_vs_matched_controls
        ),
        "external_multiscale_future_preserved_noninferior_false_rate_vs_adversarial_controls": (
            multiscale_future_preserved_noninferior_false_rate_vs_adversarial_controls
        ),
        "external_multiscale_future_preserved_beats_matched_controls": (
            multiscale_future_preserved_beats_matched_controls
        ),
        "external_multiscale_future_preserved_beats_adversarial_calibrated_controls": (
            multiscale_future_preserved_beats_adversarial_calibrated_controls
        ),
        "external_multiscale_future_preserved_claim_allowed": False,
        "external_terminal_bridge_replacement_frontier_kind": (
            "terminal_bridge_replacement_frontier_v0"
        ),
        "external_terminal_bridge_replacement_frontier_count": (
            len(replacement_frontier_rows)
        ),
        "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum": (
            replacement_frontier_native_long_range_delta_sum
        ),
        "external_terminal_bridge_replacement_frontier_probe_selected_count": (
            len(replacement_frontier_probe_selected_rows)
        ),
        "external_terminal_bridge_replacement_frontier_probe_selected_native_long_range_delta_sum": (
            replacement_frontier_probe_selected_native_long_range_delta_sum
        ),
        "external_terminal_bridge_replacement_frontier_external_count_delta_positive_count": (
            len(replacement_frontier_external_count_delta_positive_rows)
        ),
        "external_terminal_bridge_replacement_frontier_external_confidence_delta_positive_count": (
            len(replacement_frontier_external_confidence_delta_positive_rows)
        ),
        "external_terminal_bridge_replacement_frontier_external_count_delta_sum": (
            replacement_frontier_external_count_delta_sum
        ),
        "external_terminal_bridge_replacement_frontier_external_confidence_delta_sum": (
            replacement_frontier_external_confidence_delta_sum
        ),
        "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum_with_external_count_gain": (
            replacement_frontier_native_long_range_delta_sum_with_external_count_gain
        ),
        "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum_with_external_confidence_gain": (
            replacement_frontier_native_long_range_delta_sum_with_external_confidence_gain
        ),
        "external_terminal_bridge_replacement_frontier_claim_allowed": False,
        "external_persistent_rank_consistent_cluster_gated_recovery_diagnostic_kind": (
            "persistent_trace_recovery_after_rank_consistent_gate_v0"
        ),
        "external_persistent_rank_consistent_cluster_gated_recovered_event_count": (
            len(persistent_recovered_events)
        ),
        "external_persistent_rank_consistent_cluster_gated_recovered_native_contact_count": (
            sum(
                event.native_contact_count_after_scoring
                for event in persistent_recovered_events
            )
        ),
        "external_persistent_rank_consistent_cluster_gated_recovered_native_long_range_contact_count": (
            sum(
                event.native_long_range_contact_count_after_scoring
                for event in persistent_recovered_events
            )
        ),
        "external_persistent_rank_consistent_cluster_gated_recall_frontier_kind": (
            "persistent_score_margin_recall_frontier_v0"
        ),
        "external_persistent_rank_consistent_cluster_gated_recall_frontier_count": (
            len(recall_frontier_rows)
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count": (
            len(score_margin_expansion_rows)
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count": (
            score_margin_expansion_row_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_row_count": (
            score_margin_expansion_native_long_range_row_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_native_contact_count": (
            score_margin_expansion_native_contact_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_native_long_range_contact_count": (
            score_margin_expansion_native_long_range_contact_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_false_candidate_count": (
            score_margin_expansion_false_candidate_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_candidate_count": (
            max_matched_expansion_candidate_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_row_count": (
            max_matched_expansion_row_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_native_long_range_row_count": (
            max_matched_expansion_native_long_range_row_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_native_long_range_contact_count": (
            max_matched_expansion_native_long_range_contact_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_false_candidate_count": (
            max_matched_expansion_false_candidate_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_candidate_count": (
            max_adversarial_expansion_candidate_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_row_count": (
            max_adversarial_expansion_row_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_native_long_range_row_count": (
            max_adversarial_expansion_native_long_range_row_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_native_long_range_contact_count": (
            max_adversarial_expansion_native_long_range_contact_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_false_candidate_count": (
            max_adversarial_expansion_false_candidate_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count_margin_vs_matched_controls": (
            len(score_margin_expansion_rows) - max_matched_expansion_candidate_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count_margin_vs_matched_controls": (
            score_margin_expansion_row_count - max_matched_expansion_row_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_row_margin_vs_matched_controls": (
            score_margin_expansion_native_long_range_row_count
            - max_matched_expansion_native_long_range_row_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_matched_controls": (
            score_margin_expansion_native_long_range_contact_count
            - max_matched_expansion_native_long_range_contact_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count_margin_vs_adversarial_controls": (
            len(score_margin_expansion_rows)
            - max_adversarial_expansion_candidate_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count_margin_vs_adversarial_controls": (
            score_margin_expansion_row_count - max_adversarial_expansion_row_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_row_margin_vs_adversarial_controls": (
            score_margin_expansion_native_long_range_row_count
            - max_adversarial_expansion_native_long_range_row_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_adversarial_controls": (
            score_margin_expansion_native_long_range_contact_count
            - max_adversarial_expansion_native_long_range_contact_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_signal_seen": (
            score_margin_expansion_row_count >= 2
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_beats_matched_controls": (
            score_margin_expansion_row_count > max_matched_expansion_row_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_beats_adversarial_controls": (
            score_margin_expansion_row_count > max_adversarial_expansion_row_count
        ),
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_claim_allowed": False,
        "external_rank_consistent_cluster_gated_frontier_diagnostic_kind": (
            "rank_consistent_native_after_scoring_frontier_v0"
        ),
        "external_rank_consistent_cluster_gated_native_positive_frontier_count": (
            len(frontier_rows)
        ),
        "external_rank_consistent_cluster_gated_frontier_native_contact_count": (
            frontier_native_contact_count
        ),
        "external_rank_consistent_cluster_gated_frontier_native_long_range_contact_count": (
            frontier_native_long_range_contact_count
        ),
        "external_rank_consistent_cluster_gated_frontier_claim_allowed": False,
        "matched_negative_controls_present": bool(matched_controls),
        "external_claim_mode_validation_failures": claim_mode_failures,
        "coordinate_truth_used_to_build_constraints": (
            import_result.dataset.coordinate_truth_tainted
        ),
        "native_truth_used_before_coupling_selection": (
            import_result.dataset.native_truth_tainted
        ),
        "structure_model_used": import_result.dataset.structure_model_tainted,
        "oracle_constraint_control": import_result.dataset.oracle_constraint_control,
        "mechanism_discovery_claim_allowed": False,
        "global_folding_claim_allowed": False,
        "folding_problem_solved": False,
        "claim_allowed": False,
        "best_possible_v0_claim": (
            "external evolutionary couplings preserve part of the "
            "anti-fake-nucleus signal under matched controls"
        ),
        "boundary_statement": (
            "V0 tests external MSA/DCA-style coupling signal against the frozen "
            "8-row coordinate benchmark, matched negative controls, and the "
            "oracle coordinate positive control. It cannot solve folding or "
            "unlock mechanism claims on eight rows."
        ),
    }


def render_external_coupling_trace_loop_dashboard(
    report: Mapping[str, object],
) -> str:
    labels = (
        "result",
        "reason",
        "external_probe_passed",
        "external_couplings_available_rows",
        "external_constraint_count",
        "external_rows_rejected_low_depth",
        "external_rows_rejected_mapping",
        "external_real_false_nucleus_rate",
        "physical_rerank_false_nucleus_rate",
        "external_real_cluster_precision",
        "mean_matched_control_cluster_precision",
        "external_real_beats_physical",
        "external_real_beats_matched_controls",
        "external_margin_gated_false_nucleus_rate",
        "external_margin_gated_cluster_precision",
        "external_margin_gated_long_range_recall",
        "external_margin_gated_beats_matched_controls",
        "external_top_rank_gated_false_nucleus_rate",
        "external_top_rank_gated_cluster_precision",
        "external_top_rank_gated_long_range_recall",
        "external_top_rank_gated_vs_control_enrichment_ratio",
        "external_top_rank_gated_beats_matched_controls",
        "external_core_expanded_false_nucleus_rate",
        "external_core_expanded_cluster_precision",
        "external_core_expanded_long_range_recall",
        "external_core_expanded_vs_control_enrichment_ratio",
        "external_core_expanded_beats_matched_controls",
        "external_cluster_gated_core_expanded_false_nucleus_rate",
        "external_cluster_gated_core_expanded_cluster_precision",
        "external_cluster_gated_core_expanded_long_range_recall",
        "external_cluster_gated_core_expanded_vs_control_enrichment_ratio",
        "external_cluster_gated_core_expanded_mean_selected_coupling_selectivity_score",
        "external_cluster_gated_core_expanded_max_control_mean_selected_coupling_selectivity_score",
        "external_cluster_gated_core_expanded_mean_decoy_coupling_selectivity_score",
        "external_cluster_gated_core_expanded_max_control_mean_decoy_coupling_selectivity_score",
        "external_cluster_gated_core_expanded_mean_coupling_decoy_selectivity_margin",
        "external_cluster_gated_core_expanded_max_control_mean_coupling_decoy_selectivity_margin",
        "external_cluster_gated_core_expanded_beats_matched_controls",
        "external_rank_consistent_cluster_gated_false_nucleus_rate",
        "external_rank_consistent_cluster_gated_cluster_precision",
        "external_rank_consistent_cluster_gated_long_range_recall",
        "external_rank_consistent_cluster_gated_vs_control_enrichment_ratio",
        "external_rank_consistent_cluster_gated_vs_adversarial_calibrated_enrichment_ratio",
        "external_rank_consistent_cluster_gated_mean_selected_coupling_selectivity_score",
        "external_rank_consistent_cluster_gated_max_control_mean_selected_coupling_selectivity_score",
        "external_rank_consistent_cluster_gated_mean_coupling_decoy_selectivity_margin",
        "external_rank_consistent_cluster_gated_max_control_mean_coupling_decoy_selectivity_margin",
        "external_rank_consistent_cluster_gated_beats_matched_controls",
        "external_rank_consistent_cluster_gated_beats_adversarial_calibrated_controls",
        "external_rank_consistent_cluster_gated_probe_passed",
        "external_persistent_rank_consistent_cluster_gated_false_nucleus_rate",
        "external_persistent_rank_consistent_cluster_gated_cluster_precision",
        "external_persistent_rank_consistent_cluster_gated_long_range_recall",
        "external_persistent_rank_consistent_cluster_gated_vs_control_enrichment_ratio",
        "external_persistent_rank_consistent_cluster_gated_vs_adversarial_calibrated_enrichment_ratio",
        "external_persistent_rank_consistent_cluster_gated_real_vs_decoy_coupling_nucleus_enrichment_ratio",
        "external_persistent_rank_consistent_cluster_gated_vs_control_nucleus_score_enrichment_ratio",
        "external_persistent_rank_consistent_cluster_gated_vs_adversarial_calibrated_nucleus_score_enrichment_ratio",
        "external_persistent_rank_consistent_cluster_gated_beats_matched_controls",
        "external_persistent_rank_consistent_cluster_gated_beats_adversarial_calibrated_controls",
        "external_persistent_rank_consistent_cluster_gated_probe_passed",
        "external_persistent_rank_consistent_cluster_gated_selector_score_probe_passed",
        "external_persistent_rank_consistent_cluster_gated_hard_selector_score_probe_passed",
        "external_score_margin_expanded_selected_event_count",
        "external_score_margin_expanded_added_event_count",
        "external_score_margin_expanded_added_native_long_range_contact_count",
        "external_score_margin_expanded_added_false_event_count",
        "external_score_margin_expanded_false_nucleus_rate",
        "external_score_margin_expanded_cluster_precision",
        "external_score_margin_expanded_long_range_recall",
        "external_score_margin_expanded_long_range_recall_delta_vs_persistent",
        "external_score_margin_expanded_cluster_precision_margin_vs_matched_controls",
        "external_score_margin_expanded_long_range_recall_margin_vs_matched_controls",
        "external_score_margin_expanded_cluster_precision_margin_vs_adversarial_controls",
        "external_score_margin_expanded_long_range_recall_margin_vs_adversarial_controls",
        "external_score_margin_expanded_beats_matched_controls",
        "external_score_margin_expanded_beats_adversarial_calibrated_controls",
        "external_score_margin_expanded_claim_allowed",
        "external_boundary_continuity_expanded_selected_event_count",
        "external_boundary_continuity_expanded_added_event_count",
        "external_boundary_continuity_expanded_added_native_long_range_contact_count",
        "external_boundary_continuity_expanded_added_false_event_count",
        "external_boundary_continuity_expanded_false_nucleus_rate",
        "external_boundary_continuity_expanded_cluster_precision",
        "external_boundary_continuity_expanded_long_range_recall",
        "external_boundary_continuity_expanded_long_range_recall_delta_vs_score_margin",
        "external_boundary_continuity_expanded_cluster_precision_margin_vs_matched_controls",
        "external_boundary_continuity_expanded_long_range_recall_margin_vs_matched_controls",
        "external_boundary_continuity_expanded_cluster_precision_margin_vs_adversarial_controls",
        "external_boundary_continuity_expanded_long_range_recall_margin_vs_adversarial_controls",
        "external_boundary_continuity_expanded_beats_matched_controls",
        "external_boundary_continuity_expanded_beats_adversarial_calibrated_controls",
        "external_boundary_continuity_expanded_claim_allowed",
        "external_edge_continuity_expanded_selected_event_count",
        "external_edge_continuity_expanded_added_event_count",
        "external_edge_continuity_expanded_added_native_long_range_contact_count",
        "external_edge_continuity_expanded_added_false_event_count",
        "external_edge_continuity_expanded_false_nucleus_rate",
        "external_edge_continuity_expanded_cluster_precision",
        "external_edge_continuity_expanded_long_range_recall",
        "external_edge_continuity_expanded_long_range_recall_delta_vs_boundary_continuity",
        "external_edge_continuity_expanded_cluster_precision_margin_vs_matched_controls",
        "external_edge_continuity_expanded_long_range_recall_margin_vs_matched_controls",
        "external_edge_continuity_expanded_cluster_precision_margin_vs_adversarial_controls",
        "external_edge_continuity_expanded_long_range_recall_margin_vs_adversarial_controls",
        "external_edge_continuity_expanded_beats_matched_controls",
        "external_edge_continuity_expanded_beats_adversarial_calibrated_controls",
        "external_edge_continuity_expanded_claim_allowed",
        "external_pressure_release_expanded_selected_event_count",
        "external_pressure_release_expanded_added_event_count",
        "external_pressure_release_expanded_added_native_long_range_contact_count",
        "external_pressure_release_expanded_added_false_event_count",
        "external_pressure_release_expanded_false_nucleus_rate",
        "external_pressure_release_expanded_cluster_precision",
        "external_pressure_release_expanded_long_range_recall",
        "external_pressure_release_expanded_long_range_recall_delta_vs_edge_continuity",
        "external_pressure_release_expanded_cluster_precision_margin_vs_matched_controls",
        "external_pressure_release_expanded_long_range_recall_margin_vs_matched_controls",
        "external_pressure_release_expanded_cluster_precision_margin_vs_adversarial_controls",
        "external_pressure_release_expanded_long_range_recall_margin_vs_adversarial_controls",
        "external_pressure_release_expanded_beats_matched_controls",
        "external_pressure_release_expanded_beats_adversarial_calibrated_controls",
        "external_pressure_release_expanded_claim_allowed",
        "external_registry_extension_expanded_selected_event_count",
        "external_registry_extension_expanded_added_event_count",
        "external_registry_extension_expanded_added_native_long_range_contact_count",
        "external_registry_extension_expanded_added_false_event_count",
        "external_registry_extension_expanded_false_nucleus_rate",
        "external_registry_extension_expanded_cluster_precision",
        "external_registry_extension_expanded_long_range_recall",
        "external_registry_extension_expanded_long_range_recall_delta_vs_pressure_release",
        "external_registry_extension_expanded_cluster_precision_margin_vs_matched_controls",
        "external_registry_extension_expanded_long_range_recall_margin_vs_matched_controls",
        "external_registry_extension_expanded_cluster_precision_margin_vs_adversarial_controls",
        "external_registry_extension_expanded_long_range_recall_margin_vs_adversarial_controls",
        "external_registry_extension_expanded_beats_matched_controls",
        "external_registry_extension_expanded_beats_adversarial_calibrated_controls",
        "external_registry_extension_expanded_claim_allowed",
        "external_terminal_bridge_expanded_selected_event_count",
        "external_terminal_bridge_expanded_added_event_count",
        "external_terminal_bridge_expanded_added_native_long_range_contact_count",
        "external_terminal_bridge_expanded_added_false_event_count",
        "external_terminal_bridge_expanded_false_nucleus_rate",
        "external_terminal_bridge_expanded_cluster_precision",
        "external_terminal_bridge_expanded_long_range_recall",
        "external_terminal_bridge_expanded_long_range_recall_delta_vs_registry_extension",
        "external_terminal_bridge_expanded_cluster_precision_margin_vs_matched_controls",
        "external_terminal_bridge_expanded_long_range_recall_margin_vs_matched_controls",
        "external_terminal_bridge_expanded_cluster_precision_margin_vs_adversarial_controls",
        "external_terminal_bridge_expanded_long_range_recall_margin_vs_adversarial_controls",
        "external_terminal_bridge_expanded_beats_matched_controls",
        "external_terminal_bridge_expanded_beats_adversarial_calibrated_controls",
        "external_terminal_bridge_expanded_claim_allowed",
        "external_boundary_field_replacement_probe_added_event_count",
        "external_boundary_field_replacement_probe_false_nucleus_rate",
        "external_boundary_field_replacement_probe_cluster_precision",
        "external_boundary_field_replacement_probe_long_range_recall",
        "external_boundary_field_replacement_probe_long_range_recall_delta_vs_terminal_bridge",
        "external_boundary_field_replacement_probe_claim_allowed",
        "external_macro_scale_future_preserved_segment_length",
        "external_macro_scale_future_preserved_segment_stride",
        "external_macro_scale_future_preserved_selected_event_count",
        "external_macro_scale_future_preserved_false_nucleus_rate",
        "external_macro_scale_future_preserved_cluster_precision",
        "external_macro_scale_future_preserved_long_range_recall",
        "external_macro_scale_future_preserved_long_range_recall_delta_vs_boundary_replacement",
        "external_macro_scale_future_preserved_cluster_precision_margin_vs_matched_controls",
        "external_macro_scale_future_preserved_long_range_recall_margin_vs_matched_controls",
        "external_macro_scale_future_preserved_cluster_precision_margin_vs_adversarial_controls",
        "external_macro_scale_future_preserved_long_range_recall_margin_vs_adversarial_controls",
        "external_macro_scale_future_preserved_beats_matched_controls",
        "external_macro_scale_future_preserved_beats_adversarial_calibrated_controls",
        "external_macro_scale_future_preserved_claim_allowed",
        "external_multiscale_future_preserved_segment_lengths",
        "external_multiscale_future_preserved_max_events_per_row",
        "external_multiscale_future_preserved_selected_event_count",
        "external_multiscale_future_preserved_false_nucleus_rate",
        "external_multiscale_future_preserved_cluster_precision",
        "external_multiscale_future_preserved_long_range_recall",
        "external_multiscale_future_preserved_long_range_recall_delta_vs_macro",
        "external_multiscale_future_preserved_cluster_precision_margin_vs_matched_controls",
        "external_multiscale_future_preserved_long_range_recall_margin_vs_matched_controls",
        "external_multiscale_future_preserved_cluster_precision_margin_vs_adversarial_controls",
        "external_multiscale_future_preserved_long_range_recall_margin_vs_adversarial_controls",
        "external_multiscale_future_preserved_beats_matched_controls",
        "external_multiscale_future_preserved_beats_adversarial_calibrated_controls",
        "external_multiscale_future_preserved_claim_allowed",
        "external_terminal_bridge_replacement_frontier_count",
        "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum",
        "external_terminal_bridge_replacement_frontier_probe_selected_count",
        "external_terminal_bridge_replacement_frontier_probe_selected_native_long_range_delta_sum",
        "external_terminal_bridge_replacement_frontier_external_count_delta_positive_count",
        "external_terminal_bridge_replacement_frontier_external_confidence_delta_positive_count",
        "external_terminal_bridge_replacement_frontier_external_count_delta_sum",
        "external_terminal_bridge_replacement_frontier_external_confidence_delta_sum",
        "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum_with_external_count_gain",
        "external_terminal_bridge_replacement_frontier_native_long_range_delta_sum_with_external_confidence_gain",
        "external_terminal_bridge_replacement_frontier_claim_allowed",
        "external_persistent_rank_consistent_cluster_gated_recovered_event_count",
        "external_persistent_rank_consistent_cluster_gated_recovered_native_contact_count",
        "external_persistent_rank_consistent_cluster_gated_recovered_native_long_range_contact_count",
        "external_persistent_rank_consistent_cluster_gated_recall_frontier_count",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_row_count",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_native_long_range_contact_count",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_candidate_count",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_row_count",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_matched_control_native_long_range_contact_count",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count_margin_vs_matched_controls",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count_margin_vs_matched_controls",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_matched_controls",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_candidate_count",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_row_count",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_max_adversarial_native_long_range_contact_count",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count_margin_vs_adversarial_controls",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_row_count_margin_vs_adversarial_controls",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_native_long_range_margin_vs_adversarial_controls",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_signal_seen",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_beats_matched_controls",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_repeated_independent_row_beats_adversarial_controls",
        "external_persistent_rank_consistent_cluster_gated_score_margin_expansion_claim_allowed",
        "external_rank_consistent_cluster_gated_native_positive_frontier_count",
        "external_rank_consistent_cluster_gated_frontier_native_contact_count",
        "external_rank_consistent_cluster_gated_frontier_native_long_range_contact_count",
        "external_rank_consistent_cluster_gated_frontier_claim_allowed",
        "hard_adversarial_calibrated_probe_passed",
        "mechanism_discovery_claim_allowed",
        "folding_problem_solved",
    )
    cards = "".join(
        "<div class=\"metric\">"
        f"<span>{label}</span><strong>{report.get(label)}</strong>"
        "</div>"
        for label in labels
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>External Evolutionary Coupling Trace Loop V0</title>
  <style>
    body {{ margin: 0; font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; background: #f8f9f7; color: #202523; }}
    header {{ padding: 32px; background: #26362f; color: #f7fbf7; }}
    main {{ max-width: 1180px; margin: 0 auto; padding: 24px; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 12px; }}
    .metric {{ border: 1px solid #d8ded8; border-radius: 6px; background: #fff; padding: 14px; }}
    .metric span {{ display: block; color: #5d6762; font-size: 12px; overflow-wrap: anywhere; }}
    .metric strong {{ display: block; margin-top: 8px; font-size: 20px; overflow-wrap: anywhere; }}
    section {{ margin: 24px 0; }}
    code {{ background: #edf1ed; padding: 2px 5px; border-radius: 4px; }}
  </style>
</head>
<body>
  <header>
    <h1>External Evolutionary Coupling Trace Loop V0</h1>
    <p>External MSA/DCA coupling signal versus matched controls. Claims remain locked.</p>
  </header>
  <main>
    <section class="metrics">{cards}</section>
    <section>
      <h2>Boundary</h2>
      <p><code>claim_allowed</code>: {report.get('claim_allowed')}</p>
      <p>{report.get('boundary_statement')}</p>
    </section>
  </main>
</body>
</html>
"""


def run_external_evolutionary_coupling_trace_loop_benchmark(
    *,
    benchmark_file: Path,
    external_coupling_file: Path,
    oracle_coupling_file: Path,
    report_path: Path,
    certificate_path: Path,
    selectors_path: Path,
    selected_events_path: Path,
    frontier_path: Path,
    controls_path: Path,
    row_status_path: Path,
    dashboard_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path, Path]:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    import_result = import_external_coupling_dataset(
        rows=rows,
        external_coupling_file=external_coupling_file,
        policy=SERIOUS_EXTERNAL_COUPLING_POLICY,
    )
    oracle_dataset = load_coupling_dataset(oracle_coupling_file)
    external_context = build_coupling_nucleus_context(
        rows=rows,
        coupling_dataset=import_result.dataset,
    )
    physical_context = external_context.physical_context
    external_real = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_real",
        control_kind="external_real",
    )
    external_margin_gated = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_margin_gated",
        selection_mode="coupling_trace_loop_margin_gated",
        control_kind="external_real_margin_gated",
    )
    external_top_rank_gated = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_top_rank_gated",
        selection_mode="coupling_trace_loop_top_rank_gated",
        control_kind="external_real_top_rank_gated",
    )
    external_core_expanded = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_core_expanded",
        selection_mode="coupling_trace_loop_core_expanded",
        control_kind="external_real_core_expanded",
    )
    external_cluster_gated_core_expanded = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_cluster_gated_core_expanded",
        selection_mode="coupling_trace_loop_cluster_gated_core_expanded",
        control_kind="external_real_cluster_gated_core_expanded",
    )
    external_rank_consistent_cluster_gated = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_rank_consistent_cluster_gated",
        selection_mode="coupling_trace_loop_rank_consistent_cluster_gated",
        control_kind="external_real_rank_consistent_cluster_gated",
    )
    external_persistent_rank_consistent_cluster_gated = (
        _run_trace_loop_selector_from_context(
            context=external_context,
            dataset=import_result.dataset,
            selector_name="external_persistent_rank_consistent_cluster_gated",
            selection_mode=(
                "coupling_trace_loop_persistent_rank_consistent_cluster_gated"
            ),
            control_kind="external_real_persistent_rank_consistent_cluster_gated",
        )
    )
    external_score_margin_expanded = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_score_margin_expanded",
        selection_mode="coupling_trace_loop_score_margin_expanded",
        control_kind="external_real_score_margin_expanded",
    )
    external_boundary_continuity_expanded = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_boundary_continuity_expanded",
        selection_mode="coupling_trace_loop_boundary_continuity_expanded",
        control_kind="external_real_boundary_continuity_expanded",
    )
    external_edge_continuity_expanded = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_edge_continuity_expanded",
        selection_mode="coupling_trace_loop_edge_continuity_expanded",
        control_kind="external_real_edge_continuity_expanded",
    )
    external_pressure_release_expanded = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_pressure_release_expanded",
        selection_mode="coupling_trace_loop_pressure_release_expanded",
        control_kind="external_real_pressure_release_expanded",
    )
    external_registry_extension_expanded = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_registry_extension_expanded",
        selection_mode="coupling_trace_loop_registry_extension_expanded",
        control_kind="external_real_registry_extension_expanded",
    )
    external_terminal_bridge_expanded = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_terminal_bridge_expanded",
        selection_mode="coupling_trace_loop_terminal_bridge_expanded",
        control_kind="external_real_terminal_bridge_expanded",
    )
    external_boundary_field_replacement_probe = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="external_boundary_field_replacement_probe",
        selection_mode="coupling_trace_loop_boundary_field_replacement_probe",
        control_kind="external_real_boundary_field_replacement_probe",
    )
    macro_external_context = _build_macro_scale_coupling_context(
        rows=rows,
        dataset=import_result.dataset,
    )
    external_macro_scale_future_preserved = _run_trace_loop_selector_from_context(
        context=macro_external_context,
        dataset=import_result.dataset,
        selector_name="external_macro_scale_future_preserved",
        selection_mode="coupling_trace_loop_macro_scale_future_preserved",
        control_kind="external_real_macro_scale_future_preserved",
    )
    multiscale_physical_contexts = _build_multiscale_physical_contexts(rows)
    external_multiscale_future_preserved = (
        _run_multiscale_future_preserved_selector(
            rows=rows,
            dataset=import_result.dataset,
            selector_name="external_multiscale_future_preserved",
            control_kind="external_real_multiscale_future_preserved",
            physical_contexts=multiscale_physical_contexts,
        )
    )
    physical_baseline = _run_trace_loop_selector_from_context(
        context=external_context,
        dataset=import_result.dataset,
        selector_name="physical_no_coupling_baseline",
        selection_mode="physical_rerank",
        control_kind="physical_no_coupling_baseline",
    )
    controls: Mapping[str, CouplingControlDataset] = (
        generate_external_coupling_negative_controls(
            rows=rows,
            dataset=import_result.dataset,
        )
    )
    adversarial_controls: Mapping[str, CouplingControlDataset] = (
        generate_adversarial_calibrated_external_coupling_controls(
            rows=rows,
            dataset=import_result.dataset,
        )
    )
    control_contexts = {
        name: build_coupling_nucleus_context(
            rows=rows,
            coupling_dataset=control.dataset,
            physical_context=physical_context,
        )
        for name, control in controls.items()
    }
    adversarial_control_contexts = {
        name: build_coupling_nucleus_context(
            rows=rows,
            coupling_dataset=control.dataset,
            physical_context=physical_context,
        )
        for name, control in adversarial_controls.items()
    }
    macro_physical_context = macro_external_context.physical_context
    macro_control_contexts = {
        name: build_coupling_nucleus_context(
            rows=rows,
            coupling_dataset=control.dataset,
            physical_context=macro_physical_context,
        )
        for name, control in controls.items()
    }
    adversarial_macro_control_contexts = {
        name: build_coupling_nucleus_context(
            rows=rows,
            coupling_dataset=control.dataset,
            physical_context=macro_physical_context,
        )
        for name, control in adversarial_controls.items()
    }
    matched_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=name,
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    margin_gated_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_margin_gated_{name}",
            selection_mode="coupling_trace_loop_margin_gated",
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    top_rank_gated_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_top_rank_gated_{name}",
            selection_mode="coupling_trace_loop_top_rank_gated",
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    core_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_core_expanded_{name}",
            selection_mode="coupling_trace_loop_core_expanded",
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    cluster_gated_core_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_cluster_gated_core_expanded_{name}",
            selection_mode="coupling_trace_loop_cluster_gated_core_expanded",
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    rank_consistent_cluster_gated_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_rank_consistent_cluster_gated_{name}",
            selection_mode="coupling_trace_loop_rank_consistent_cluster_gated",
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    persistent_rank_consistent_cluster_gated_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_persistent_rank_consistent_cluster_gated_{name}",
            selection_mode=(
                "coupling_trace_loop_persistent_rank_consistent_cluster_gated"
            ),
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    adversarial_rank_consistent_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=adversarial_control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_rank_consistent_cluster_gated_{name}",
            selection_mode="coupling_trace_loop_rank_consistent_cluster_gated",
            control_kind=control.control_kind,
        )
        for name, control in adversarial_controls.items()
    )
    adversarial_persistent_rank_consistent_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=adversarial_control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_persistent_rank_consistent_cluster_gated_{name}",
            selection_mode=(
                "coupling_trace_loop_persistent_rank_consistent_cluster_gated"
            ),
            control_kind=control.control_kind,
        )
        for name, control in adversarial_controls.items()
    )
    score_margin_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_score_margin_expanded_{name}",
            selection_mode="coupling_trace_loop_score_margin_expanded",
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    adversarial_score_margin_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=adversarial_control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_score_margin_expanded_{name}",
            selection_mode="coupling_trace_loop_score_margin_expanded",
            control_kind=control.control_kind,
        )
        for name, control in adversarial_controls.items()
    )
    boundary_continuity_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_boundary_continuity_expanded_{name}",
            selection_mode="coupling_trace_loop_boundary_continuity_expanded",
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    adversarial_boundary_continuity_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=adversarial_control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_boundary_continuity_expanded_{name}",
            selection_mode="coupling_trace_loop_boundary_continuity_expanded",
            control_kind=control.control_kind,
        )
        for name, control in adversarial_controls.items()
    )
    edge_continuity_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_edge_continuity_expanded_{name}",
            selection_mode="coupling_trace_loop_edge_continuity_expanded",
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    adversarial_edge_continuity_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=adversarial_control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_edge_continuity_expanded_{name}",
            selection_mode="coupling_trace_loop_edge_continuity_expanded",
            control_kind=control.control_kind,
        )
        for name, control in adversarial_controls.items()
    )
    pressure_release_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_pressure_release_expanded_{name}",
            selection_mode="coupling_trace_loop_pressure_release_expanded",
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    adversarial_pressure_release_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=adversarial_control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_pressure_release_expanded_{name}",
            selection_mode="coupling_trace_loop_pressure_release_expanded",
            control_kind=control.control_kind,
        )
        for name, control in adversarial_controls.items()
    )
    registry_extension_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_registry_extension_expanded_{name}",
            selection_mode="coupling_trace_loop_registry_extension_expanded",
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    adversarial_registry_extension_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=adversarial_control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_registry_extension_expanded_{name}",
            selection_mode="coupling_trace_loop_registry_extension_expanded",
            control_kind=control.control_kind,
        )
        for name, control in adversarial_controls.items()
    )
    terminal_bridge_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_terminal_bridge_expanded_{name}",
            selection_mode="coupling_trace_loop_terminal_bridge_expanded",
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    adversarial_terminal_bridge_expanded_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=adversarial_control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_terminal_bridge_expanded_{name}",
            selection_mode="coupling_trace_loop_terminal_bridge_expanded",
            control_kind=control.control_kind,
        )
        for name, control in adversarial_controls.items()
    )
    macro_scale_future_preserved_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=macro_control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_macro_scale_future_preserved_{name}",
            selection_mode="coupling_trace_loop_macro_scale_future_preserved",
            control_kind=control.control_kind,
        )
        for name, control in controls.items()
    )
    adversarial_macro_scale_future_preserved_control_runs = tuple(
        _run_trace_loop_selector_from_context(
            context=adversarial_macro_control_contexts[name],
            dataset=control.dataset,
            selector_name=f"external_macro_scale_future_preserved_{name}",
            selection_mode="coupling_trace_loop_macro_scale_future_preserved",
            control_kind=control.control_kind,
        )
        for name, control in adversarial_controls.items()
    )
    multiscale_future_preserved_control_runs = tuple(
        _run_multiscale_future_preserved_selector(
            rows=rows,
            dataset=control.dataset,
            selector_name=f"external_multiscale_future_preserved_{name}",
            control_kind=control.control_kind,
            physical_contexts=multiscale_physical_contexts,
        )
        for name, control in controls.items()
    )
    adversarial_multiscale_future_preserved_control_runs = tuple(
        _run_multiscale_future_preserved_selector(
            rows=rows,
            dataset=control.dataset,
            selector_name=f"external_multiscale_future_preserved_{name}",
            control_kind=control.control_kind,
            physical_contexts=multiscale_physical_contexts,
        )
        for name, control in adversarial_controls.items()
    )
    oracle_context = build_coupling_nucleus_context(
        rows=rows,
        coupling_dataset=oracle_dataset,
        physical_context=physical_context,
    )
    oracle_positive_control = _run_trace_loop_selector_from_context(
        context=oracle_context,
        dataset=oracle_dataset,
        selector_name="oracle_coordinate_positive_control",
        control_kind="oracle_coordinate_positive_control",
    )
    all_runs = (
        external_real,
        external_margin_gated,
        external_top_rank_gated,
        external_core_expanded,
        external_cluster_gated_core_expanded,
        external_rank_consistent_cluster_gated,
        external_persistent_rank_consistent_cluster_gated,
        external_score_margin_expanded,
        external_boundary_continuity_expanded,
        external_edge_continuity_expanded,
        external_pressure_release_expanded,
        external_registry_extension_expanded,
        external_terminal_bridge_expanded,
        external_boundary_field_replacement_probe,
        external_macro_scale_future_preserved,
        external_multiscale_future_preserved,
        physical_baseline,
        *matched_control_runs,
        *margin_gated_control_runs,
        *top_rank_gated_control_runs,
        *core_expanded_control_runs,
        *cluster_gated_core_expanded_control_runs,
        *rank_consistent_cluster_gated_control_runs,
        *adversarial_rank_consistent_control_runs,
        *persistent_rank_consistent_cluster_gated_control_runs,
        *adversarial_persistent_rank_consistent_control_runs,
        *score_margin_expanded_control_runs,
        *adversarial_score_margin_expanded_control_runs,
        *boundary_continuity_expanded_control_runs,
        *adversarial_boundary_continuity_expanded_control_runs,
        *edge_continuity_expanded_control_runs,
        *adversarial_edge_continuity_expanded_control_runs,
        *pressure_release_expanded_control_runs,
        *adversarial_pressure_release_expanded_control_runs,
        *registry_extension_expanded_control_runs,
        *adversarial_registry_extension_expanded_control_runs,
        *terminal_bridge_expanded_control_runs,
        *adversarial_terminal_bridge_expanded_control_runs,
        *macro_scale_future_preserved_control_runs,
        *adversarial_macro_scale_future_preserved_control_runs,
        *multiscale_future_preserved_control_runs,
        *adversarial_multiscale_future_preserved_control_runs,
        oracle_positive_control,
    )
    frontier_rows = rank_consistent_frontier_rows(
        context=external_context,
        cluster_gated_run=external_cluster_gated_core_expanded,
        rank_consistent_run=external_rank_consistent_cluster_gated,
    )
    recall_frontier_rows = persistent_recall_frontier_rows(
        context=external_context,
        persistent_run=external_persistent_rank_consistent_cluster_gated,
    )
    replacement_frontier_rows = terminal_bridge_replacement_frontier_rows(
        context=external_context,
        terminal_run=external_terminal_bridge_expanded,
        replacement_probe_run=external_boundary_field_replacement_probe,
    )
    persistent_control_run_by_name = {
        run.selector_name.replace(
            "external_persistent_rank_consistent_cluster_gated_",
            "",
            1,
        ): run
        for run in persistent_rank_consistent_cluster_gated_control_runs
    }
    matched_control_recall_frontier_summaries = tuple(
        _score_margin_expansion_summary(
            persistent_recall_frontier_rows(
                context=control_contexts[name],
                persistent_run=persistent_control_run_by_name[name],
            )
        )
        for name in controls
    )
    adversarial_persistent_control_run_by_name = {
        run.selector_name.replace(
            "external_persistent_rank_consistent_cluster_gated_",
            "",
            1,
        ): run
        for run in adversarial_persistent_rank_consistent_control_runs
    }
    adversarial_recall_frontier_summaries = tuple(
        _score_margin_expansion_summary(
            persistent_recall_frontier_rows(
                context=adversarial_control_contexts[name],
                persistent_run=adversarial_persistent_control_run_by_name[name],
            )
        )
        for name in adversarial_controls
    )
    report = _build_report(
        rows=rows,
        import_result=import_result,
        external_real=external_real,
        external_margin_gated=external_margin_gated,
        external_top_rank_gated=external_top_rank_gated,
        external_core_expanded=external_core_expanded,
        external_cluster_gated_core_expanded=external_cluster_gated_core_expanded,
        external_rank_consistent_cluster_gated=(
            external_rank_consistent_cluster_gated
        ),
        external_persistent_rank_consistent_cluster_gated=(
            external_persistent_rank_consistent_cluster_gated
        ),
        external_score_margin_expanded=external_score_margin_expanded,
        external_boundary_continuity_expanded=(
            external_boundary_continuity_expanded
        ),
        external_edge_continuity_expanded=external_edge_continuity_expanded,
        external_pressure_release_expanded=external_pressure_release_expanded,
        external_registry_extension_expanded=(
            external_registry_extension_expanded
        ),
        external_terminal_bridge_expanded=external_terminal_bridge_expanded,
        external_boundary_field_replacement_probe=(
            external_boundary_field_replacement_probe
        ),
        external_macro_scale_future_preserved=external_macro_scale_future_preserved,
        macro_scale_future_preserved_controls=(
            macro_scale_future_preserved_control_runs
        ),
        adversarial_macro_scale_future_preserved_controls=(
            adversarial_macro_scale_future_preserved_control_runs
        ),
        external_multiscale_future_preserved=(
            external_multiscale_future_preserved
        ),
        multiscale_future_preserved_controls=(
            multiscale_future_preserved_control_runs
        ),
        adversarial_multiscale_future_preserved_controls=(
            adversarial_multiscale_future_preserved_control_runs
        ),
        physical_baseline=physical_baseline,
        matched_controls=matched_control_runs,
        margin_gated_controls=margin_gated_control_runs,
        top_rank_gated_controls=top_rank_gated_control_runs,
        core_expanded_controls=core_expanded_control_runs,
        cluster_gated_core_expanded_controls=(
            cluster_gated_core_expanded_control_runs
        ),
        rank_consistent_cluster_gated_controls=(
            rank_consistent_cluster_gated_control_runs
        ),
        adversarial_rank_consistent_controls=(
            adversarial_rank_consistent_control_runs
        ),
        persistent_rank_consistent_controls=(
            persistent_rank_consistent_cluster_gated_control_runs
        ),
        adversarial_persistent_rank_consistent_controls=(
            adversarial_persistent_rank_consistent_control_runs
        ),
        score_margin_expanded_controls=score_margin_expanded_control_runs,
        adversarial_score_margin_expanded_controls=(
            adversarial_score_margin_expanded_control_runs
        ),
        boundary_continuity_expanded_controls=(
            boundary_continuity_expanded_control_runs
        ),
        adversarial_boundary_continuity_expanded_controls=(
            adversarial_boundary_continuity_expanded_control_runs
        ),
        edge_continuity_expanded_controls=edge_continuity_expanded_control_runs,
        adversarial_edge_continuity_expanded_controls=(
            adversarial_edge_continuity_expanded_control_runs
        ),
        pressure_release_expanded_controls=pressure_release_expanded_control_runs,
        adversarial_pressure_release_expanded_controls=(
            adversarial_pressure_release_expanded_control_runs
        ),
        registry_extension_expanded_controls=(
            registry_extension_expanded_control_runs
        ),
        adversarial_registry_extension_expanded_controls=(
            adversarial_registry_extension_expanded_control_runs
        ),
        terminal_bridge_expanded_controls=terminal_bridge_expanded_control_runs,
        adversarial_terminal_bridge_expanded_controls=(
            adversarial_terminal_bridge_expanded_control_runs
        ),
        oracle_positive_control=oracle_positive_control,
        frontier_rows=frontier_rows,
        recall_frontier_rows=recall_frontier_rows,
        replacement_frontier_rows=replacement_frontier_rows,
        matched_control_recall_frontier_summaries=(
            matched_control_recall_frontier_summaries
        ),
        adversarial_recall_frontier_summaries=adversarial_recall_frontier_summaries,
        source_benchmark_file=benchmark_file,
        external_coupling_file=external_coupling_file,
        oracle_coupling_file=oracle_coupling_file,
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    certificate_path.write_text(
        json.dumps(_certificate(report), ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    write_csv_rows([run.metric.to_dict() for run in all_runs], selectors_path)
    selected_rows: list[dict[str, object]] = []
    for run in all_runs:
        selected_rows.extend(run.selected_rows)
    write_csv_rows(selected_rows, selected_events_path)
    write_csv_rows(
        [*replacement_frontier_rows, *frontier_rows, *recall_frontier_rows],
        frontier_path,
    )
    write_csv_rows(_controls_from_runs(all_runs), controls_path)
    write_csv_rows(
        [status.to_dict() for status in import_result.row_statuses],
        row_status_path,
    )
    dashboard_path.write_text(
        render_external_coupling_trace_loop_dashboard(report),
        encoding="utf-8",
    )
    return (
        report_path,
        certificate_path,
        selectors_path,
        selected_events_path,
        frontier_path,
        controls_path,
        row_status_path,
        dashboard_path,
    )
