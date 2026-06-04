from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import mean
from typing import Mapping, Sequence

from pharmacotopology.artifact_io import write_csv_rows
from pharmacotopology.folding_coupling_decoy_falsification import (
    COUPLING_DECOY_FALSIFICATION_KIND,
    coupling_decoy_comparisons,
    real_beats_decoy_coupling_score_rate,
    real_vs_decoy_coupling_enrichment_ratio,
)
from pharmacotopology.folding_evolutionary_constraints import (
    COUPLING_ASSESSMENT_KIND,
    EVOLUTIONARY_COUPLING_LAYER_KIND,
    CouplingClosureAssessment,
    CouplingConstraint,
    CouplingDataset,
    assess_coupling_closures,
    compatible_future_event,
    load_coupling_dataset,
)
from pharmacotopology.folding_external_coupling_sources import (
    ACCEPTED_EXTERNAL_COUPLING_SOURCE_KINDS,
    EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID,
)
from pharmacotopology.folding_nucleus_closure_search import (
    FOLDING_NUCLEUS_CLOSURE_EVENT_KIND,
    NucleusClosureEvent,
)
from pharmacotopology.folding_nucleus_decoy_falsification import (
    decoy_distance,
    matched_decoys_for_selected_events,
)
from pharmacotopology.folding_physical_selection import (
    ACTIVE_PHYSICAL_SELECTION_SCORE_KIND,
    ActivePhysicalContext,
    active_physical_score,
    build_active_physical_context,
    select_events as select_active_physical_events,
)
from pharmacotopology.folding_physical_state import (
    PHYSICAL_CLOSURE_STATE_KIND,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


COUPLING_NUCLEUS_SELECTOR_REPORT_KIND = (
    "coupling_preserved_nucleus_selector_benchmark_v1"
)
COUPLING_NUCLEUS_SELECTOR_CERTIFICATE_KIND = (
    "coupling_preserved_nucleus_selector_certificate"
)
COUPLING_NUCLEUS_SELECTOR_SCORE_KIND = (
    "coupling_future_physical_decoy_margin_score_v1"
)
EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_KIND = (
    EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_ID
)

SELECTED_EVENTS_PER_ROW = 40
COUPLING_DIRECT_SUPPORT_MIN = 0.22
COUPLING_FUTURE_PRESERVATION_MIN = 0.62
COUPLING_BLOCKED_FUTURE_MAX = 0.16
COUPLING_DECOY_MARGIN_MIN = 0.04
TRACE_LOOP_MARGIN_GATE_MIN = 0.0
TRACE_LOOP_MARGIN_GATE_BLOCKED_FUTURE_MAX = 0.16

SURVIVAL_FALSE_RATE_MAX = 0.25
SURVIVAL_CLUSTER_PRECISION_MIN = 0.08
SURVIVAL_LONG_RANGE_RECALL_MIN = 0.35
SURVIVAL_DECOY_ENRICHMENT_MIN = 1.35

ROOT_OUTPUT_NAMES = (
    "coupling_nucleus_selector_report.json",
    "coupling_nucleus_selector_selectors.csv",
    "coupling_nucleus_selector_selected_events.csv",
    "coupling_nucleus_selector_assessments.csv",
    "coupling_nucleus_selector_decoys.csv",
    "coupling_nucleus_selector_dashboard.html",
    "coupling_nucleus_selector_certificate.json",
)
EXTERNAL_EVOLUTIONARY_COUPLING_SOURCE_KINDS = ACCEPTED_EXTERNAL_COUPLING_SOURCE_KINDS


@dataclass(frozen=True)
class CouplingNucleusContext:
    rows: tuple[RealCoordinateVisualRow, ...]
    physical_context: ActivePhysicalContext
    coupling_dataset: CouplingDataset
    assessments: tuple[CouplingClosureAssessment, ...]
    assessment_by_event_id: Mapping[str, CouplingClosureAssessment]
    coupling_decoy_margin_by_event_id: Mapping[str, float]

    @property
    def competitive_events(self) -> tuple[NucleusClosureEvent, ...]:
        return self.physical_context.competitive_events

    @property
    def graph_events(self) -> tuple[NucleusClosureEvent, ...]:
        return self.physical_context.graph_events

    @property
    def event_by_id(self) -> Mapping[str, NucleusClosureEvent]:
        return self.physical_context.event_by_id


@dataclass(frozen=True)
class CouplingSelectorMetric:
    selector_name: str
    selected_event_count: int
    false_nucleus_rate: float
    contact_cluster_precision: float
    long_range_contact_recall: float
    coupling_constraint_recall: float
    real_vs_decoy_coupling_enrichment_ratio: float
    real_beats_decoy_coupling_score_rate: float
    mean_coupling_nucleus_score: float
    survives_targets: bool
    coordinate_truth_used_to_build_constraints: bool
    native_truth_used_before_coupling_selection: bool
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(value, 6)


def _events_by_row(events: Sequence[NucleusClosureEvent]) -> dict[str, list[NucleusClosureEvent]]:
    output: dict[str, list[NucleusClosureEvent]] = defaultdict(list)
    for event in events:
        output[event.row_id].append(event)
    return output


def _region_union(events: Sequence[NucleusClosureEvent]) -> set[tuple[int, int]]:
    pairs: set[tuple[int, int]] = set()
    for event in events:
        pairs.update(event.candidate_region_pairs())
    return pairs


def _constraint_confidence_by_pair(
    constraints: Sequence[CouplingConstraint],
) -> dict[tuple[int, int], float]:
    return {constraint.pair(): constraint.confidence for constraint in constraints}


def build_coupling_nucleus_context(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    coupling_dataset: CouplingDataset,
) -> CouplingNucleusContext:
    physical_context = build_active_physical_context(rows)
    assessments = assess_coupling_closures(
        rows=rows,
        events=physical_context.competitive_events,
        dataset=coupling_dataset,
    )
    assessment_by_event = {assessment.event_id: assessment for assessment in assessments}
    by_row = _events_by_row(physical_context.competitive_events)
    margins: dict[str, float] = {}
    for row_events in by_row.values():
        for event in row_events:
            decoy_candidates = [
                candidate
                for candidate in row_events
                if candidate.event_id != event.event_id
            ]
            if not decoy_candidates:
                margins[event.event_id] = 0.0
                continue
            decoy = min(decoy_candidates, key=lambda item: decoy_distance(event, item))
            margins[event.event_id] = _rounded(
                assessment_by_event[event.event_id].coupling_selectivity_score
                - assessment_by_event[decoy.event_id].coupling_selectivity_score
            )
    return CouplingNucleusContext(
        rows=tuple(rows),
        physical_context=physical_context,
        coupling_dataset=coupling_dataset,
        assessments=assessments,
        assessment_by_event_id=assessment_by_event,
        coupling_decoy_margin_by_event_id=margins,
    )


def coupling_nucleus_score(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> float:
    assessment = context.assessment_by_event_id[event.event_id]
    state = context.physical_context.state_by_event_id[event.event_id]
    return _rounded(
        0.62 * assessment.coupling_selectivity_score
        + 0.18 * active_physical_score(event, context.physical_context)
        + 0.12 * state.burial_gain
        + 0.08 * event.registry_support
        + 0.08 * event.contact_cluster_gain
        + 0.10 * context.coupling_decoy_margin_by_event_id[event.event_id]
        - 0.14 * state.unsatisfied_polar_penalty
        - 0.10 * state.steric_clash_score
        - 0.12 * assessment.blocked_future_pressure
    )


def _passes_coupling_future_gate(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> bool:
    assessment = context.assessment_by_event_id[event.event_id]
    return (
        assessment.direct_support_score >= COUPLING_DIRECT_SUPPORT_MIN
        and assessment.future_preservation_score >= COUPLING_FUTURE_PRESERVATION_MIN
        and assessment.blocked_future_pressure <= COUPLING_BLOCKED_FUTURE_MAX
    )


def select_coupling_events(
    context: CouplingNucleusContext,
    *,
    selector_name: str,
) -> tuple[NucleusClosureEvent, ...]:
    if selector_name == "graph_only":
        return context.graph_events
    if selector_name == "physical_rerank":
        return select_active_physical_events(
            context.physical_context,
            selector_name="physical_rerank",
        )
    if selector_name == "coupling_trace_loop":
        return select_coupling_trace_loop_events(context)
    if selector_name == "coupling_trace_loop_margin_gated":
        return select_coupling_trace_loop_events(
            context,
            min_coupling_decoy_margin=TRACE_LOOP_MARGIN_GATE_MIN,
            max_blocked_future_pressure=TRACE_LOOP_MARGIN_GATE_BLOCKED_FUTURE_MAX,
        )

    selected: list[NucleusClosureEvent] = []
    competitive_by_row = _events_by_row(context.competitive_events)
    for row in context.rows:
        candidates = tuple(competitive_by_row.get(row.row_id, ()))
        if selector_name in {"coupling_future_gate", "coupling_decoy_falsifier"}:
            candidates = tuple(
                event
                for event in candidates
                if _passes_coupling_future_gate(event, context)
            )
        if selector_name == "coupling_decoy_falsifier":
            candidates = tuple(
                event
                for event in candidates
                if context.coupling_decoy_margin_by_event_id[event.event_id]
                >= COUPLING_DECOY_MARGIN_MIN
            )
        selected.extend(
            sorted(
                candidates,
                key=lambda event: (
                    -coupling_nucleus_score(event, context),
                    -context.assessment_by_event_id[
                        event.event_id
                    ].direct_support_score,
                    event.segment_a_start,
                    event.segment_b_start,
                    event.event_id,
                ),
            )[:SELECTED_EVENTS_PER_ROW]
        )
    return tuple(selected)


def select_coupling_trace_loop_events(
    context: CouplingNucleusContext,
    *,
    min_coupling_decoy_margin: float | None = None,
    max_blocked_future_pressure: float | None = None,
) -> tuple[NucleusClosureEvent, ...]:
    constraints_by_row = context.coupling_dataset.constraints_by_row_id()
    competitive_by_row = _events_by_row(context.competitive_events)
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        row_constraints = tuple(constraints_by_row.get(row.row_id, ()))
        confidence_by_pair = _constraint_confidence_by_pair(row_constraints)
        uncovered = set(confidence_by_pair)
        row_selected: list[NucleusClosureEvent] = []
        row_candidates = tuple(competitive_by_row.get(row.row_id, ()))
        while uncovered and len(row_selected) < SELECTED_EVENTS_PER_ROW:
            scored_candidates: list[tuple[float, NucleusClosureEvent]] = []
            for event in row_candidates:
                if event.event_id in {item.event_id for item in row_selected}:
                    continue
                if any(
                    not compatible_future_event(selected_event, event)
                    for selected_event in row_selected
                ):
                    continue
                event_pairs = set(event.candidate_region_pairs())
                newly_covered = event_pairs & uncovered
                if not newly_covered:
                    continue
                assessment = context.assessment_by_event_id[event.event_id]
                if (
                    max_blocked_future_pressure is not None
                    and assessment.blocked_future_pressure
                    > max_blocked_future_pressure
                ):
                    continue
                coupling_decoy_margin = context.coupling_decoy_margin_by_event_id[
                    event.event_id
                ]
                if (
                    min_coupling_decoy_margin is not None
                    and coupling_decoy_margin < min_coupling_decoy_margin
                ):
                    continue
                covered_confidence = sum(
                    confidence_by_pair[pair] for pair in newly_covered
                )
                trace_score = _rounded(
                    0.48 * min(1.0, covered_confidence / 3.0)
                    + 0.28 * coupling_nucleus_score(event, context)
                    + 0.16 * assessment.future_preservation_score
                    - 0.10 * assessment.blocked_future_pressure
                )
                scored_candidates.append((trace_score, event))
            if not scored_candidates:
                break
            _, chosen = max(
                scored_candidates,
                key=lambda item: (
                    item[0],
                    context.assessment_by_event_id[
                        item[1].event_id
                    ].direct_support_score,
                    -item[1].segment_a_start,
                    -item[1].segment_b_start,
                    item[1].event_id,
                ),
            )
            row_selected.append(chosen)
            uncovered -= set(chosen.candidate_region_pairs())
        selected.extend(row_selected)
    return tuple(selected)


def selector_metrics(
    context: CouplingNucleusContext,
    *,
    selector_name: str,
    selected_events: Sequence[NucleusClosureEvent],
) -> CouplingSelectorMetric:
    by_row = _events_by_row(selected_events)
    constraints_by_row = context.coupling_dataset.constraints_by_row_id()
    false_rates: list[float] = []
    precisions: list[float] = []
    long_recalls: list[float] = []
    constraint_recalls: list[float] = []
    for row in context.rows:
        row_events = tuple(by_row.get(row.row_id, ()))
        native_pairs = set(row.native_contact_pairs())
        native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
        constraint_pairs = {
            constraint.pair()
            for constraint in constraints_by_row.get(row.row_id, ())
        }
        region = _region_union(row_events)
        native_hit_count = sum(
            event.native_contact_count_after_scoring for event in row_events
        )
        possible_region_pair_count = len(row_events) * 64
        false_count = sum(
            1 for event in row_events if event.native_contact_count_after_scoring == 0
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
    matches = matched_decoys_for_selected_events(
        selected_events=selected_events,
        candidate_events=context.competitive_events,
    )
    comparisons = coupling_decoy_comparisons(
        matches=matches,
        assessments=context.assessments,
    )
    false_rate = _rounded(mean(false_rates) if false_rates else 0.0)
    precision = _rounded(mean(precisions) if precisions else 0.0)
    long_recall = _rounded(mean(long_recalls) if long_recalls else 0.0)
    enrichment = real_vs_decoy_coupling_enrichment_ratio(comparisons)
    beats = real_beats_decoy_coupling_score_rate(comparisons)
    survives = (
        false_rate < SURVIVAL_FALSE_RATE_MAX
        and precision > SURVIVAL_CLUSTER_PRECISION_MIN
        and long_recall > SURVIVAL_LONG_RANGE_RECALL_MIN
        and enrichment > SURVIVAL_DECOY_ENRICHMENT_MIN
    )
    return CouplingSelectorMetric(
        selector_name=selector_name,
        selected_event_count=len(selected_events),
        false_nucleus_rate=false_rate,
        contact_cluster_precision=precision,
        long_range_contact_recall=long_recall,
        coupling_constraint_recall=_rounded(
            mean(constraint_recalls) if constraint_recalls else 0.0
        ),
        real_vs_decoy_coupling_enrichment_ratio=enrichment,
        real_beats_decoy_coupling_score_rate=beats,
        mean_coupling_nucleus_score=_rounded(
            mean([coupling_nucleus_score(event, context) for event in selected_events])
            if selected_events
            else 0.0
        ),
        survives_targets=survives,
        coordinate_truth_used_to_build_constraints=(
            context.coupling_dataset.coordinate_truth_tainted
        ),
        native_truth_used_before_coupling_selection=(
            context.coupling_dataset.native_truth_used_before_coupling_selection
        ),
    )


def selected_event_rows(
    context: CouplingNucleusContext,
    selections: Mapping[str, Sequence[NucleusClosureEvent]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for selector_name in sorted(selections):
        for rank, event in enumerate(
            sorted(
                selections[selector_name],
                key=lambda item: (
                    item.row_id,
                    -coupling_nucleus_score(item, context),
                    item.segment_a_start,
                    item.segment_b_start,
                    item.event_id,
                ),
            ),
            start=1,
        ):
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
                }
            )
    return rows


def decoy_rows(
    context: CouplingNucleusContext,
    selected_events: Sequence[NucleusClosureEvent],
) -> list[dict[str, object]]:
    matches = matched_decoys_for_selected_events(
        selected_events=selected_events,
        candidate_events=context.competitive_events,
    )
    return [
        row.to_dict()
        for row in coupling_decoy_comparisons(
            matches=matches,
            assessments=context.assessments,
        )
    ]


def coupling_claim_mode_validation_failures(
    dataset: CouplingDataset,
) -> tuple[str, ...]:
    failures: list[str] = []
    if not dataset.external_evolutionary_couplings_used:
        failures.append("external_evolutionary_couplings_used=false")
    if dataset.coupling_source_kind not in EXTERNAL_EVOLUTIONARY_COUPLING_SOURCE_KINDS:
        failures.append(f"coupling_source_kind={dataset.coupling_source_kind}")
    if dataset.coordinate_truth_tainted:
        failures.append("coordinate_truth_used_to_build_constraints=true")
    if dataset.native_truth_tainted:
        failures.append("native_truth_used_before_coupling_selection=true")
    if dataset.structure_model_tainted:
        failures.append("structure_model_used=true")
    if dataset.oracle_constraint_control:
        failures.append("oracle_constraint_control=true")
    return tuple(failures)


def build_coupling_nucleus_selector_report(
    *,
    context: CouplingNucleusContext,
    selector_rows: Sequence[CouplingSelectorMetric],
    source_benchmark_file: Path,
    coupling_file: Path,
) -> dict[str, object]:
    selector_lookup = {row.selector_name: row for row in selector_rows}
    coupling_target_survives = any(row.survives_targets for row in selector_rows)
    claim_mode_failures = coupling_claim_mode_validation_failures(
        context.coupling_dataset
    )
    claim_mode_validation_passed = not claim_mode_failures
    claim_allowed = coupling_target_survives and claim_mode_validation_passed
    return {
        "report_kind": COUPLING_NUCLEUS_SELECTOR_REPORT_KIND,
        "source_benchmark_kind": REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "coupling_layer_kind": EVOLUTIONARY_COUPLING_LAYER_KIND,
        "coupling_file": str(coupling_file),
        "external_batch_kind": EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_BATCH_KIND,
        "source_event_kind": FOLDING_NUCLEUS_CLOSURE_EVENT_KIND,
        "physical_state_kind": PHYSICAL_CLOSURE_STATE_KIND,
        "active_physical_selection_score_kind": ACTIVE_PHYSICAL_SELECTION_SCORE_KIND,
        "coupling_assessment_kind": COUPLING_ASSESSMENT_KIND,
        "coupling_decoy_falsification_kind": COUPLING_DECOY_FALSIFICATION_KIND,
        "coupling_nucleus_selector_score_kind": COUPLING_NUCLEUS_SELECTOR_SCORE_KIND,
        "benchmark_size": len(context.rows),
        "candidate_event_count": len(context.competitive_events),
        "coupling_constraint_count": len(context.coupling_dataset.constraints),
        "graph_only_false_nucleus_rate": selector_lookup[
            "graph_only"
        ].false_nucleus_rate,
        "physical_rerank_false_nucleus_rate": selector_lookup[
            "physical_rerank"
        ].false_nucleus_rate,
        "coupling_rerank_false_nucleus_rate": selector_lookup[
            "coupling_rerank"
        ].false_nucleus_rate,
        "coupling_trace_loop_false_nucleus_rate": selector_lookup[
            "coupling_trace_loop"
        ].false_nucleus_rate,
        "coupling_future_gate_false_nucleus_rate": selector_lookup[
            "coupling_future_gate"
        ].false_nucleus_rate,
        "coupling_decoy_falsifier_false_nucleus_rate": selector_lookup[
            "coupling_decoy_falsifier"
        ].false_nucleus_rate,
        "graph_only_cluster_precision": selector_lookup[
            "graph_only"
        ].contact_cluster_precision,
        "physical_rerank_cluster_precision": selector_lookup[
            "physical_rerank"
        ].contact_cluster_precision,
        "coupling_rerank_cluster_precision": selector_lookup[
            "coupling_rerank"
        ].contact_cluster_precision,
        "coupling_trace_loop_cluster_precision": selector_lookup[
            "coupling_trace_loop"
        ].contact_cluster_precision,
        "coupling_future_gate_cluster_precision": selector_lookup[
            "coupling_future_gate"
        ].contact_cluster_precision,
        "coupling_decoy_falsifier_cluster_precision": selector_lookup[
            "coupling_decoy_falsifier"
        ].contact_cluster_precision,
        "graph_only_long_range_recall": selector_lookup[
            "graph_only"
        ].long_range_contact_recall,
        "physical_rerank_long_range_recall": selector_lookup[
            "physical_rerank"
        ].long_range_contact_recall,
        "coupling_rerank_long_range_recall": selector_lookup[
            "coupling_rerank"
        ].long_range_contact_recall,
        "coupling_rerank_constraint_recall": selector_lookup[
            "coupling_rerank"
        ].coupling_constraint_recall,
        "coupling_trace_loop_long_range_recall": selector_lookup[
            "coupling_trace_loop"
        ].long_range_contact_recall,
        "coupling_trace_loop_constraint_recall": selector_lookup[
            "coupling_trace_loop"
        ].coupling_constraint_recall,
        "coupling_future_gate_long_range_recall": selector_lookup[
            "coupling_future_gate"
        ].long_range_contact_recall,
        "coupling_future_gate_constraint_recall": selector_lookup[
            "coupling_future_gate"
        ].coupling_constraint_recall,
        "coupling_decoy_falsifier_long_range_recall": selector_lookup[
            "coupling_decoy_falsifier"
        ].long_range_contact_recall,
        "coupling_decoy_falsifier_constraint_recall": selector_lookup[
            "coupling_decoy_falsifier"
        ].coupling_constraint_recall,
        "coupling_rerank_real_vs_decoy_enrichment_ratio": selector_lookup[
            "coupling_rerank"
        ].real_vs_decoy_coupling_enrichment_ratio,
        "coupling_trace_loop_real_vs_decoy_enrichment_ratio": selector_lookup[
            "coupling_trace_loop"
        ].real_vs_decoy_coupling_enrichment_ratio,
        "coupling_future_gate_real_vs_decoy_enrichment_ratio": selector_lookup[
            "coupling_future_gate"
        ].real_vs_decoy_coupling_enrichment_ratio,
        "coupling_decoy_falsifier_real_vs_decoy_enrichment_ratio": selector_lookup[
            "coupling_decoy_falsifier"
        ].real_vs_decoy_coupling_enrichment_ratio,
        "coupling_selector_targets_met": coupling_target_survives,
        "external_evolutionary_couplings_used": (
            context.coupling_dataset.external_evolutionary_couplings_used
        ),
        "coupling_source_kind": context.coupling_dataset.coupling_source_kind,
        "per_constraint_coordinate_truth_used": (
            context.coupling_dataset.per_constraint_coordinate_truth_used
        ),
        "coordinate_truth_used_to_build_constraints": (
            context.coupling_dataset.coordinate_truth_tainted
        ),
        "native_truth_used_before_coupling_selection": (
            context.coupling_dataset.native_truth_used_before_coupling_selection
        ),
        "oracle_constraint_control": context.coupling_dataset.oracle_constraint_control,
        "claim_mode_validation_passed": claim_mode_validation_passed,
        "claim_mode_validation_failures": claim_mode_failures,
        "mechanism_discovery_claim_allowed": claim_allowed,
        "mechanism_discovery_claim_created": False,
        "global_folding_claim_allowed": claim_allowed,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "claim_allowed": claim_allowed,
        "boundary_statement": (
            "Coupling-supported selection tests the missing native-selective "
            "constraint channel by asking whether direct coupling support, "
            "future closure preservation, and matched-decoy coupling margins "
            "can separate plausible nuclei from fake burial. The checked-in "
            "constraint file is an oracle control unless replaced by external "
            "MSA/DCA couplings, so it cannot create a folding mechanism claim."
        ),
        "selectors": [row.to_dict() for row in selector_rows],
    }


def build_coupling_nucleus_selector_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": COUPLING_NUCLEUS_SELECTOR_CERTIFICATE_KIND,
        "report_kind": report["report_kind"],
        "coupling_selector_targets_met": report["coupling_selector_targets_met"],
        "coupling_rerank_false_nucleus_rate": report[
            "coupling_rerank_false_nucleus_rate"
        ],
        "coupling_trace_loop_false_nucleus_rate": report[
            "coupling_trace_loop_false_nucleus_rate"
        ],
        "coupling_rerank_cluster_precision": report[
            "coupling_rerank_cluster_precision"
        ],
        "coupling_trace_loop_cluster_precision": report[
            "coupling_trace_loop_cluster_precision"
        ],
        "coupling_rerank_long_range_recall": report[
            "coupling_rerank_long_range_recall"
        ],
        "coupling_trace_loop_long_range_recall": report[
            "coupling_trace_loop_long_range_recall"
        ],
        "coupling_trace_loop_constraint_recall": report[
            "coupling_trace_loop_constraint_recall"
        ],
        "coupling_rerank_real_vs_decoy_enrichment_ratio": report[
            "coupling_rerank_real_vs_decoy_enrichment_ratio"
        ],
        "coupling_trace_loop_real_vs_decoy_enrichment_ratio": report[
            "coupling_trace_loop_real_vs_decoy_enrichment_ratio"
        ],
        "external_evolutionary_couplings_used": report[
            "external_evolutionary_couplings_used"
        ],
        "coordinate_truth_used_to_build_constraints": report[
            "coordinate_truth_used_to_build_constraints"
        ],
        "native_truth_used_before_coupling_selection": report[
            "native_truth_used_before_coupling_selection"
        ],
        "per_constraint_coordinate_truth_used": report[
            "per_constraint_coordinate_truth_used"
        ],
        "oracle_constraint_control": report["oracle_constraint_control"],
        "claim_mode_validation_passed": report["claim_mode_validation_passed"],
        "claim_mode_validation_failures": tuple(
            report["claim_mode_validation_failures"]  # type: ignore[index]
        ),
        "mechanism_discovery_claim_allowed": report[
            "mechanism_discovery_claim_allowed"
        ],
        "global_folding_claim_allowed": report["global_folding_claim_allowed"],
        "folding_problem_solved": report["folding_problem_solved"],
        "claim_allowed": report["claim_allowed"],
        "output_artifacts": tuple(output_names),
    }


def write_coupling_nucleus_selector_outputs(
    *,
    report: Mapping[str, object],
    selector_rows: Sequence[CouplingSelectorMetric],
    selected_rows: Sequence[Mapping[str, object]],
    assessments: Sequence[CouplingClosureAssessment],
    decoys: Sequence[Mapping[str, object]],
    report_path: Path,
    selectors_path: Path,
    selected_events_path: Path,
    assessments_path: Path,
    decoys_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_csv_rows([row.to_dict() for row in selector_rows], selectors_path)
    write_csv_rows(selected_rows, selected_events_path)
    write_csv_rows([row.to_dict() for row in assessments], assessments_path)
    write_csv_rows(decoys, decoys_path)
    dashboard_path.write_text(
        render_coupling_nucleus_selector_dashboard(report),
        encoding="utf-8",
    )
    certificate = build_coupling_nucleus_selector_certificate(
        report,
        output_names=ROOT_OUTPUT_NAMES,
    )
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return (
        report_path,
        selectors_path,
        selected_events_path,
        assessments_path,
        decoys_path,
        dashboard_path,
        certificate_path,
    )


def run_coupling_nucleus_selector_benchmark(
    *,
    benchmark_file: Path,
    coupling_file: Path,
    report_path: Path,
    selectors_path: Path,
    selected_events_path: Path,
    assessments_path: Path,
    decoys_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path]:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    coupling_dataset = load_coupling_dataset(coupling_file)
    context = build_coupling_nucleus_context(
        rows=rows,
        coupling_dataset=coupling_dataset,
    )
    selections = {
        "graph_only": select_coupling_events(context, selector_name="graph_only"),
        "physical_rerank": select_coupling_events(
            context,
            selector_name="physical_rerank",
        ),
        "coupling_rerank": select_coupling_events(
            context,
            selector_name="coupling_rerank",
        ),
        "coupling_trace_loop": select_coupling_events(
            context,
            selector_name="coupling_trace_loop",
        ),
        "coupling_future_gate": select_coupling_events(
            context,
            selector_name="coupling_future_gate",
        ),
        "coupling_decoy_falsifier": select_coupling_events(
            context,
            selector_name="coupling_decoy_falsifier",
        ),
    }
    selector_rows = tuple(
        selector_metrics(
            context,
            selector_name=name,
            selected_events=events,
        )
        for name, events in selections.items()
    )
    selected_rows = selected_event_rows(context, selections)
    decoys = decoy_rows(context, selections["coupling_rerank"])
    report = build_coupling_nucleus_selector_report(
        context=context,
        selector_rows=selector_rows,
        source_benchmark_file=benchmark_file,
        coupling_file=coupling_file,
    )
    return write_coupling_nucleus_selector_outputs(
        report=report,
        selector_rows=selector_rows,
        selected_rows=selected_rows,
        assessments=context.assessments,
        decoys=decoys,
        report_path=report_path,
        selectors_path=selectors_path,
        selected_events_path=selected_events_path,
        assessments_path=assessments_path,
        decoys_path=decoys_path,
        dashboard_path=dashboard_path,
        certificate_path=certificate_path,
    )


def _escape(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _metric_cards(report: Mapping[str, object]) -> str:
    labels = (
        "coupling_rerank_false_nucleus_rate",
        "coupling_rerank_cluster_precision",
        "coupling_rerank_long_range_recall",
        "coupling_rerank_real_vs_decoy_enrichment_ratio",
        "coupling_trace_loop_false_nucleus_rate",
        "coupling_trace_loop_long_range_recall",
        "coupling_future_gate_false_nucleus_rate",
        "coupling_future_gate_long_range_recall",
        "coupling_decoy_falsifier_false_nucleus_rate",
        "coupling_decoy_falsifier_real_vs_decoy_enrichment_ratio",
        "coupling_selector_targets_met",
        "claim_mode_validation_passed",
        "oracle_constraint_control",
        "per_constraint_coordinate_truth_used",
        "mechanism_discovery_claim_allowed",
        "folding_problem_solved",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
    )


def _selector_table(report: Mapping[str, object]) -> str:
    rows = report.get("selectors", [])
    if not isinstance(rows, list):
        rows = []
    body = "".join(
        "<tr>"
        f"<td>{_escape(row.get('selector_name'))}</td>"
        f"<td>{_escape(row.get('selected_event_count'))}</td>"
        f"<td>{_escape(row.get('false_nucleus_rate'))}</td>"
        f"<td>{_escape(row.get('contact_cluster_precision'))}</td>"
        f"<td>{_escape(row.get('long_range_contact_recall'))}</td>"
        f"<td>{_escape(row.get('coupling_constraint_recall'))}</td>"
        f"<td>{_escape(row.get('real_vs_decoy_coupling_enrichment_ratio'))}</td>"
        f"<td>{_escape(row.get('real_beats_decoy_coupling_score_rate'))}</td>"
        "</tr>"
        for row in rows
        if isinstance(row, Mapping)
    )
    return (
        "<section><h2>Selector Comparison</h2>"
        "<table><thead><tr><th>selector</th><th>events</th><th>false</th>"
        "<th>precision</th><th>long-range</th><th>constraint recall</th><th>decoy ratio</th>"
        "<th>beats decoy</th></tr></thead><tbody>"
        + body
        + "</tbody></table></section>"
    )


def render_coupling_nucleus_selector_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Coupling-Preserved Nucleus Selector</title>
  <style>
    body {{
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f7f8f5;
      color: #1f2522;
    }}
    header {{
      padding: 34px;
      background: #20332b;
      color: #f8fbf5;
    }}
    main {{
      max-width: 1240px;
      margin: 0 auto;
      padding: 24px;
    }}
    h1, h2, h3 {{
      margin: 0 0 12px;
      letter-spacing: 0;
    }}
    section {{
      margin: 24px 0;
    }}
    .metrics {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
      gap: 12px;
    }}
    .metric {{
      background: #ffffff;
      border: 1px solid #d7ded8;
      border-radius: 6px;
      padding: 14px;
    }}
    .metric span {{
      display: block;
      color: #5c6661;
      font-size: 12px;
      overflow-wrap: anywhere;
    }}
    .metric strong {{
      display: block;
      margin-top: 8px;
      font-size: 20px;
      overflow-wrap: anywhere;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: #ffffff;
      border: 1px solid #d7ded8;
    }}
    th, td {{
      padding: 9px 10px;
      border-bottom: 1px solid #e3e9e4;
      text-align: left;
      font-size: 13px;
      vertical-align: top;
      overflow-wrap: anywhere;
    }}
  </style>
</head>
<body>
  <header>
    <h1>Coupling-Preserved Nucleus Selector</h1>
    <p>Closure candidates are ranked by direct coupling support, future closure preservation, physical state, and matched-decoy margin.</p>
  </header>
  <main>
    <section class="metrics">{_metric_cards(report)}</section>
    {_selector_table(report)}
    <section>
      <h2>Claim Boundary</h2>
      <p>{_escape(report.get("boundary_statement", ""))}</p>
    </section>
  </main>
</body>
</html>
"""
