from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Mapping, Optional, Sequence

from pharmacotopology.artifact_io import write_csv_rows
from pharmacotopology.folding_coupling_negative_controls import (
    CouplingControlDataset,
    generate_external_coupling_negative_controls,
)
from pharmacotopology.folding_coupling_nucleus_selector import (
    CouplingNucleusContext,
    CouplingSelectorMetric,
    build_coupling_nucleus_context,
    coupling_claim_mode_validation_failures,
    select_coupling_events,
    selected_event_rows,
    selector_metrics,
)
from pharmacotopology.folding_evolutionary_constraints import (
    CouplingDataset,
    load_coupling_dataset,
)
from pharmacotopology.folding_external_coupling_importer import (
    ExternalCouplingImportResult,
    import_external_coupling_dataset,
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


@dataclass(frozen=True)
class TraceLoopRun:
    selector_name: str
    dataset: CouplingDataset
    metric: CouplingSelectorMetric
    selected_rows: tuple[dict[str, object], ...]
    constraint_count: int
    control_kind: str


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
        selected_rows=tuple(rows_out),
        constraint_count=len(dataset.constraints),
        control_kind=control_kind,
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
        )
        rows.append(row)
    return rows


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
    physical_baseline: TraceLoopRun,
    matched_controls: Sequence[TraceLoopRun],
    margin_gated_controls: Sequence[TraceLoopRun],
    top_rank_gated_controls: Sequence[TraceLoopRun],
    core_expanded_controls: Sequence[TraceLoopRun],
    cluster_gated_core_expanded_controls: Sequence[TraceLoopRun],
    rank_consistent_cluster_gated_controls: Sequence[TraceLoopRun],
    oracle_positive_control: TraceLoopRun,
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
    external_probe_passed = (
        acceptance_criteria_met or rank_consistent_acceptance_criteria_met
    )
    result = classify_external_probe_result(
        available_rows=available_rows,
        external_constraint_count=len(import_result.dataset.constraints),
        external_real_beats_physical=(
            external_real_beats_physical or rank_consistent_beats_physical
        ),
        external_real_beats_matched_controls=(
            external_real_beats_matched_controls
            or rank_consistent_beats_matched_controls
        ),
    )
    reason = (
        "provenance-calibrated external couplings beat physical and matched controls"
        if rank_consistent_acceptance_criteria_met
        else _external_probe_reason(
            result=result,
            selected_event_count=external_selected_event_count,
            external_constraint_count=len(import_result.dataset.constraints),
        )
    )
    external_metric_defined = external_selected_event_count > 0
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
        "external_rank_consistent_cluster_gated_meets_oracle_recall_floor": (
            rank_consistent_meets_oracle_recall_floor
        ),
        "external_rank_consistent_cluster_gated_claim_allowed": False,
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
        "external_rank_consistent_cluster_gated_mean_selected_coupling_selectivity_score",
        "external_rank_consistent_cluster_gated_max_control_mean_selected_coupling_selectivity_score",
        "external_rank_consistent_cluster_gated_mean_coupling_decoy_selectivity_margin",
        "external_rank_consistent_cluster_gated_max_control_mean_coupling_decoy_selectivity_margin",
        "external_rank_consistent_cluster_gated_beats_matched_controls",
        "external_rank_consistent_cluster_gated_probe_passed",
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
    controls_path: Path,
    row_status_path: Path,
    dashboard_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path]:
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
    control_contexts = {
        name: build_coupling_nucleus_context(
            rows=rows,
            coupling_dataset=control.dataset,
            physical_context=physical_context,
        )
        for name, control in controls.items()
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
        physical_baseline,
        *matched_control_runs,
        *margin_gated_control_runs,
        *top_rank_gated_control_runs,
        *core_expanded_control_runs,
        *cluster_gated_core_expanded_control_runs,
        *rank_consistent_cluster_gated_control_runs,
        oracle_positive_control,
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
        oracle_positive_control=oracle_positive_control,
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
        controls_path,
        row_status_path,
        dashboard_path,
    )
