from __future__ import annotations

import csv
import json
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

from pharmacotopology.folding_closure_geometry import (
    CLOSURE_GEOMETRY_KIND,
    closure_compatibility,
    compatibility_rows,
)
from pharmacotopology.folding_contact_law_features import contact_law_feature_rows
from pharmacotopology.folding_frustration_filter import (
    FRUSTRATION_FILTER_KIND,
    FrustrationAssessment,
    assess_event_frustration,
)
from pharmacotopology.folding_nucleus_closure_search import (
    FOLDING_NUCLEUS_CLOSURE_EVENT_KIND,
    NucleusClosureEvent,
    accepted_events,
    build_folding_nucleus_closure_report,
    nucleus_closure_events,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


COMPETITIVE_NUCLEUS_SELECTION_REPORT_KIND = "competitive_nucleus_selection_v1"
COMPETITIVE_NUCLEUS_SELECTION_CERTIFICATE_KIND = (
    "competitive_nucleus_selection_certificate"
)

PRE_COMPETITION_THRESHOLD = 0.30
MAX_SELECTED_EVENTS_PER_ROW = 100
MAX_SELECTED_EVENTS_PER_100_RESIDUES = 95
TARGET_POST_SELECTED_EVENTS_MAX = 800
TARGET_FALSE_NUCLEUS_RATE_MAX = 0.45
TARGET_CONTACT_CLUSTER_PRECISION_MIN = 0.08
TARGET_LONG_RANGE_RECALL_MIN = 0.45

ROOT_OUTPUT_NAMES = (
    "competitive_nucleus_selection_report.json",
    "competitive_nucleus_selection_selected_events.csv",
    "competitive_nucleus_selection_rejections.csv",
    "competitive_nucleus_selection_compatibility.csv",
    "competitive_nucleus_selection_trajectory.csv",
    "competitive_nucleus_selection_metrics.csv",
    "competitive_nucleus_selection_dashboard.html",
    "competitive_nucleus_selection_certificate.json",
)


@dataclass(frozen=True)
class CompetitiveEventDecision:
    row_id: str
    source_accession: str
    event_id: str
    competition_score: float
    false_contact_risk_proxy: float
    cluster_precision_proxy: float
    selected: bool
    selected_rank: int
    rejection_reason: str
    compatibility_rejection_label: str
    native_contact_count_after_scoring: int
    native_long_range_contact_count_after_scoring: int
    native_label_attached_after_selection: bool
    native_truth_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class CompetitiveRowMetrics:
    row_id: str
    source_accession: str
    pre_competition_event_count: int
    post_competition_selected_event_count: int
    event_reduction_ratio: float
    post_false_nucleus_rate: float
    post_contact_cluster_precision: float
    post_native_nucleus_recall: float
    post_long_range_contact_recall: float
    frustration_rejection_count: int
    geometry_rejection_count: int
    competition_rejection_count: int
    overlap_rejection_count: int
    trap_rejection_count: int
    long_range_target_met: bool
    false_rate_target_met: bool
    precision_target_met: bool

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(value, 6)


def _mean(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)


def _region_union(events: Sequence[NucleusClosureEvent]) -> set[tuple[int, int]]:
    pairs: set[tuple[int, int]] = set()
    for event in events:
        pairs.update(event.candidate_region_pairs())
    return pairs


def competition_score(
    event: NucleusClosureEvent,
    assessment: FrustrationAssessment | None = None,
) -> float:
    precision_proxy = (
        assessment.cluster_precision_proxy
        if assessment is not None
        else 0.0
    )
    return _rounded(
        0.54 * event.nucleus_score
        + 0.24 * event.closure_event_stability
        + 0.34 * event.contact_cluster_gain
        + 0.22 * event.registry_support
        + 0.10 * event.hydrophobic_burial_gain
        + 0.08 * event.secondary_structure_compatibility
        + 0.14 * precision_proxy
        - 0.24 * event.loop_entropy_cost
        - 0.12 * event.isolation_penalty
        - 0.08 * event.geometry_violation_cost
        - 0.06 * event.frustration_cost
    )


def selection_budget_for_row(row: RealCoordinateVisualRow) -> int:
    per_length_budget = round(
        row.sequence_length * MAX_SELECTED_EVENTS_PER_100_RESIDUES / 100
    )
    return max(8, min(MAX_SELECTED_EVENTS_PER_ROW, per_length_budget))


def _compatibility_rejection_label(
    event: NucleusClosureEvent,
    selected_events: Sequence[NucleusClosureEvent],
) -> str:
    labels = Counter(
        closure_compatibility(event, selected).compatibility_label
        for selected in selected_events
        if selected.event_id != event.event_id
    )
    for label in (
        "overlapping",
        "topologically_conflicting",
        "sterically_risky",
        "domain_inconsistent",
        "competing",
    ):
        if labels[label]:
            return label
    return "lower_rank_no_direct_conflict"


def select_competitive_events_for_row(
    row: RealCoordinateVisualRow,
    row_events: Sequence[NucleusClosureEvent],
) -> tuple[tuple[NucleusClosureEvent, ...], tuple[CompetitiveEventDecision, ...]]:
    assessments = {event.event_id: assess_event_frustration(event) for event in row_events}
    scored = sorted(
        row_events,
        key=lambda event: (
            -competition_score(event, assessments[event.event_id]),
            event.segment_a_start,
            event.segment_b_start,
            event.event_id,
        ),
    )
    budget = selection_budget_for_row(row)
    eligible = [
        event
        for event in scored
        if assessments[event.event_id].passed_filter
    ]
    selected = tuple(eligible[:budget])
    selected_ids = {event.event_id for event in selected}
    selected_ranks = {event.event_id: index for index, event in enumerate(selected, start=1)}
    decisions: list[CompetitiveEventDecision] = []
    for event in scored:
        assessment = assessments[event.event_id]
        selected_event = event.event_id in selected_ids
        rejection_reason = "selected" if selected_event else "competition_rejection"
        compatibility_label = "selected"
        if not selected_event:
            if not assessment.passed_filter:
                rejection_reason = assessment.primary_rejection_reason
                compatibility_label = "filter_rejection"
            else:
                compatibility_label = _compatibility_rejection_label(event, selected)
                if compatibility_label == "overlapping":
                    rejection_reason = "overlap_rejection"
                elif compatibility_label in {
                    "topologically_conflicting",
                    "sterically_risky",
                    "domain_inconsistent",
                }:
                    rejection_reason = "geometry_rejection"
        decisions.append(
            CompetitiveEventDecision(
                row_id=event.row_id,
                source_accession=event.source_accession,
                event_id=event.event_id,
                competition_score=competition_score(event, assessment),
                false_contact_risk_proxy=assessment.false_contact_risk_proxy,
                cluster_precision_proxy=assessment.cluster_precision_proxy,
                selected=selected_event,
                selected_rank=selected_ranks.get(event.event_id, 0),
                rejection_reason=rejection_reason,
                compatibility_rejection_label=compatibility_label,
                native_contact_count_after_scoring=(
                    event.native_contact_count_after_scoring
                ),
                native_long_range_contact_count_after_scoring=(
                    event.native_long_range_contact_count_after_scoring
                ),
                native_label_attached_after_selection=True,
            )
        )
    return selected, tuple(decisions)


def select_competitive_events(
    rows: Sequence[RealCoordinateVisualRow],
    events: Sequence[NucleusClosureEvent],
) -> tuple[tuple[NucleusClosureEvent, ...], tuple[CompetitiveEventDecision, ...]]:
    by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in events:
        by_row.setdefault(event.row_id, []).append(event)
    selected: list[NucleusClosureEvent] = []
    decisions: list[CompetitiveEventDecision] = []
    for row in rows:
        row_selected, row_decisions = select_competitive_events_for_row(
            row,
            tuple(by_row.get(row.row_id, ())),
        )
        selected.extend(row_selected)
        decisions.extend(row_decisions)
    return tuple(selected), tuple(decisions)


def row_metrics(
    row: RealCoordinateVisualRow,
    selected: Sequence[NucleusClosureEvent],
    decisions: Sequence[CompetitiveEventDecision],
) -> CompetitiveRowMetrics:
    native_pairs = set(row.native_contact_pairs())
    native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
    region_pairs = _region_union(selected)
    native_covered = region_pairs & native_pairs
    native_long_covered = region_pairs & native_long
    possible_region_pair_count = sum(event.candidate_contact_count for event in selected)
    native_hit_count = sum(
        event.native_contact_count_after_scoring for event in selected
    )
    selected_false = [
        event for event in selected if event.native_contact_count_after_scoring == 0
    ]
    reason_counts = Counter(decision.rejection_reason for decision in decisions)
    pre_count = len(decisions)
    post_count = len(selected)
    post_false = len(selected_false) / post_count if post_count else 0.0
    precision = (
        native_hit_count / possible_region_pair_count
        if possible_region_pair_count
        else 0.0
    )
    long_range = (
        len(native_long_covered) / len(native_long) if native_long else 1.0
    )
    return CompetitiveRowMetrics(
        row_id=row.row_id,
        source_accession=row.source_accession,
        pre_competition_event_count=pre_count,
        post_competition_selected_event_count=post_count,
        event_reduction_ratio=_rounded(1.0 - (post_count / pre_count if pre_count else 0.0)),
        post_false_nucleus_rate=_rounded(post_false),
        post_contact_cluster_precision=_rounded(precision),
        post_native_nucleus_recall=_rounded(
            len(native_covered) / len(native_pairs) if native_pairs else 0.0
        ),
        post_long_range_contact_recall=_rounded(long_range),
        frustration_rejection_count=reason_counts["frustration_rejection"],
        geometry_rejection_count=reason_counts["geometry_rejection"],
        competition_rejection_count=reason_counts["competition_rejection"],
        overlap_rejection_count=reason_counts["overlap_rejection"],
        trap_rejection_count=reason_counts["trap_rejection"],
        long_range_target_met=long_range > TARGET_LONG_RANGE_RECALL_MIN,
        false_rate_target_met=post_false < TARGET_FALSE_NUCLEUS_RATE_MAX,
        precision_target_met=precision > TARGET_CONTACT_CLUSTER_PRECISION_MIN,
    )


def competitive_row_metrics(
    rows: Sequence[RealCoordinateVisualRow],
    selected_events: Sequence[NucleusClosureEvent],
    decisions: Sequence[CompetitiveEventDecision],
) -> tuple[CompetitiveRowMetrics, ...]:
    selected_by_row: dict[str, list[NucleusClosureEvent]] = {}
    decisions_by_row: dict[str, list[CompetitiveEventDecision]] = {}
    for event in selected_events:
        selected_by_row.setdefault(event.row_id, []).append(event)
    for decision in decisions:
        decisions_by_row.setdefault(decision.row_id, []).append(decision)
    return tuple(
        row_metrics(
            row,
            tuple(selected_by_row.get(row.row_id, ())),
            tuple(decisions_by_row.get(row.row_id, ())),
        )
        for row in rows
    )


def aggregate_row_metrics(
    metrics: Sequence[CompetitiveRowMetrics],
) -> dict[str, object]:
    post_count = sum(metric.post_competition_selected_event_count for metric in metrics)
    pre_count = sum(metric.pre_competition_event_count for metric in metrics)
    return {
        "post_competition_selected_event_count": post_count,
        "event_reduction_ratio": _rounded(
            1.0 - (post_count / pre_count if pre_count else 0.0)
        ),
        "post_false_nucleus_rate": _rounded(
            _mean([metric.post_false_nucleus_rate for metric in metrics])
        ),
        "post_contact_cluster_precision": _rounded(
            _mean([metric.post_contact_cluster_precision for metric in metrics])
        ),
        "post_native_nucleus_recall": _rounded(
            _mean([metric.post_native_nucleus_recall for metric in metrics])
        ),
        "post_long_range_contact_recall": _rounded(
            _mean([metric.post_long_range_contact_recall for metric in metrics])
        ),
        "frustration_rejection_count": sum(
            metric.frustration_rejection_count for metric in metrics
        ),
        "geometry_rejection_count": sum(
            metric.geometry_rejection_count for metric in metrics
        ),
        "competition_rejection_count": sum(
            metric.competition_rejection_count for metric in metrics
        ),
        "overlap_rejection_count": sum(
            metric.overlap_rejection_count for metric in metrics
        ),
        "trap_rejection_count": sum(metric.trap_rejection_count for metric in metrics),
    }


def selected_event_rows(
    selected_events: Sequence[NucleusClosureEvent],
    decisions: Sequence[CompetitiveEventDecision],
) -> list[dict[str, object]]:
    by_event = {decision.event_id: decision for decision in decisions}
    rows: list[dict[str, object]] = []
    for event in sorted(
        selected_events,
        key=lambda item: (
            item.row_id,
            by_event[item.event_id].selected_rank,
            item.event_id,
        ),
    ):
        decision = by_event[event.event_id]
        rows.append(
            {
                "row_id": event.row_id,
                "source_accession": event.source_accession,
                "event_id": event.event_id,
                "selected_rank": decision.selected_rank,
                "competition_score": decision.competition_score,
                "false_contact_risk_proxy": decision.false_contact_risk_proxy,
                "cluster_precision_proxy": decision.cluster_precision_proxy,
                "segment_a_start": event.segment_a_start,
                "segment_a_end": event.segment_a_end,
                "segment_b_start": event.segment_b_start,
                "segment_b_end": event.segment_b_end,
                "normalized_span": event.normalized_span,
                "contact_cluster_gain": event.contact_cluster_gain,
                "registry_support": event.registry_support,
                "loop_entropy_cost": event.loop_entropy_cost,
                "geometry_violation_cost": event.geometry_violation_cost,
                "frustration_cost": event.frustration_cost,
                "native_contact_count_after_scoring": (
                    event.native_contact_count_after_scoring
                ),
                "native_long_range_contact_count_after_scoring": (
                    event.native_long_range_contact_count_after_scoring
                ),
                "native_truth_used_before_selection": False,
                "raw_sequence_exposed": False,
            }
        )
    return rows


def rejection_rows(
    decisions: Sequence[CompetitiveEventDecision],
) -> list[dict[str, object]]:
    return [
        decision.to_dict()
        for decision in decisions
        if not decision.selected
    ]


def metric_rows_from_report(report: Mapping[str, object]) -> list[dict[str, object]]:
    rows = report.get("metrics", [])
    if not isinstance(rows, list):
        return []
    return [dict(row) for row in rows if isinstance(row, Mapping)]


def trajectory_rows(
    rows: Sequence[RealCoordinateVisualRow],
    selected_events: Sequence[NucleusClosureEvent],
    decisions: Sequence[CompetitiveEventDecision],
) -> list[dict[str, object]]:
    selected_by_row: dict[str, list[NucleusClosureEvent]] = {}
    decision_by_event = {decision.event_id: decision for decision in decisions}
    for event in selected_events:
        selected_by_row.setdefault(event.row_id, []).append(event)
    output: list[dict[str, object]] = []
    for row in rows:
        native_pairs = set(row.native_contact_pairs())
        native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
        cumulative: set[tuple[int, int]] = set()
        ordered = sorted(
            selected_by_row.get(row.row_id, ()),
            key=lambda event: (
                decision_by_event[event.event_id].selected_rank,
                event.event_id,
            ),
        )
        false_count = 0
        native_hit_count = 0
        possible_region_pair_count = 0
        for step, event in enumerate(ordered, start=1):
            cumulative.update(event.candidate_region_pairs())
            if event.native_contact_count_after_scoring == 0:
                false_count += 1
            native_hit_count += event.native_contact_count_after_scoring
            possible_region_pair_count += event.candidate_contact_count
            output.append(
                {
                    "row_id": row.row_id,
                    "source_accession": row.source_accession,
                    "step": step,
                    "event_id": event.event_id,
                    "competition_score": (
                        decision_by_event[event.event_id].competition_score
                    ),
                    "cumulative_false_nucleus_rate": _rounded(false_count / step),
                    "cumulative_contact_cluster_precision": _rounded(
                        native_hit_count / possible_region_pair_count
                        if possible_region_pair_count
                        else 0.0
                    ),
                    "cumulative_native_nucleus_recall": _rounded(
                        len(cumulative & native_pairs) / len(native_pairs)
                        if native_pairs
                        else 0.0
                    ),
                    "cumulative_long_range_contact_recall": _rounded(
                        len(cumulative & native_long) / len(native_long)
                        if native_long
                        else 1.0
                    ),
                    "native_truth_used_before_selection": False,
                    "raw_sequence_exposed": False,
                }
            )
    return output


def build_competitive_nucleus_selection_report(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    events: Sequence[NucleusClosureEvent],
    selected_events: Sequence[NucleusClosureEvent],
    decisions: Sequence[CompetitiveEventDecision],
    source_benchmark_file: Path,
    pre_competition_report: Mapping[str, object],
) -> dict[str, object]:
    metrics = competitive_row_metrics(rows, selected_events, decisions)
    aggregate = aggregate_row_metrics(metrics)
    targets = {
        "selected_event_count_target_met": (
            int(aggregate["post_competition_selected_event_count"])
            < TARGET_POST_SELECTED_EVENTS_MAX
        ),
        "false_nucleus_rate_target_met": (
            float(aggregate["post_false_nucleus_rate"])
            < TARGET_FALSE_NUCLEUS_RATE_MAX
        ),
        "contact_cluster_precision_target_met": (
            float(aggregate["post_contact_cluster_precision"])
            > TARGET_CONTACT_CLUSTER_PRECISION_MIN
        ),
        "long_range_contact_recall_target_met": (
            float(aggregate["post_long_range_contact_recall"])
            > TARGET_LONG_RANGE_RECALL_MIN
        ),
    }
    law_survives = all(bool(value) for value in targets.values())
    return {
        "report_kind": COMPETITIVE_NUCLEUS_SELECTION_REPORT_KIND,
        "source_benchmark_kind": REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "source_event_kind": FOLDING_NUCLEUS_CLOSURE_EVENT_KIND,
        "closure_geometry_kind": CLOSURE_GEOMETRY_KIND,
        "frustration_filter_kind": FRUSTRATION_FILTER_KIND,
        "benchmark_size": len(rows),
        "pre_competition_threshold": PRE_COMPETITION_THRESHOLD,
        "selection_budget_max_events_per_row": MAX_SELECTED_EVENTS_PER_ROW,
        "selection_budget_max_events_per_100_residues": (
            MAX_SELECTED_EVENTS_PER_100_RESIDUES
        ),
        "candidate_closure_event_count": len(events),
        "pre_competition_event_count": pre_competition_report[
            "accepted_event_count"
        ],
        "pre_false_nucleus_rate": pre_competition_report["false_nucleus_rate"],
        "pre_contact_cluster_precision": pre_competition_report[
            "contact_cluster_precision"
        ],
        "pre_long_range_contact_recall": pre_competition_report[
            "long_range_contact_recall_after_nucleus"
        ],
        **aggregate,
        **targets,
        "nucleus_competition_law_survives": law_survives,
        "competition_reduces_event_flood": (
            int(aggregate["post_competition_selected_event_count"])
            < int(pre_competition_report["accepted_event_count"])
        ),
        "native_truth_used_before_event_generation": False,
        "native_truth_used_before_selection": False,
        "native_label_attached_after_selection": True,
        "row_specific_nucleus_thresholds_forbidden": True,
        "mechanism_discovery_claim_allowed": False,
        "mechanism_discovery_claim_created": False,
        "global_folding_claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "claim_allowed": False,
        "boundary_statement": (
            "Competitive nucleus selection reduces the closure-event flood "
            "and audits geometry, overlap, frustration, and trap pressure. "
            "This batch does not discover a folding law: false closures and "
            "cluster precision remain explicit survival gates."
        ),
        "metrics": [metric.to_dict() for metric in metrics],
    }


def build_competitive_nucleus_selection_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": COMPETITIVE_NUCLEUS_SELECTION_CERTIFICATE_KIND,
        "report_kind": report["report_kind"],
        "pre_competition_event_count": report["pre_competition_event_count"],
        "post_competition_selected_event_count": report[
            "post_competition_selected_event_count"
        ],
        "event_reduction_ratio": report["event_reduction_ratio"],
        "pre_false_nucleus_rate": report["pre_false_nucleus_rate"],
        "post_false_nucleus_rate": report["post_false_nucleus_rate"],
        "pre_contact_cluster_precision": report["pre_contact_cluster_precision"],
        "post_contact_cluster_precision": report["post_contact_cluster_precision"],
        "pre_long_range_contact_recall": report["pre_long_range_contact_recall"],
        "post_long_range_contact_recall": report["post_long_range_contact_recall"],
        "nucleus_competition_law_survives": report[
            "nucleus_competition_law_survives"
        ],
        "native_truth_used_before_event_generation": report[
            "native_truth_used_before_event_generation"
        ],
        "native_truth_used_before_selection": report[
            "native_truth_used_before_selection"
        ],
        "mechanism_discovery_claim_allowed": report[
            "mechanism_discovery_claim_allowed"
        ],
        "global_folding_claim_allowed": report["global_folding_claim_allowed"],
        "folding_problem_solved": report["folding_problem_solved"],
        "claim_allowed": report["claim_allowed"],
        "output_artifacts": tuple(output_names),
    }


def write_competitive_nucleus_selection_outputs(
    *,
    report: Mapping[str, object],
    selected_events: Sequence[NucleusClosureEvent],
    decisions: Sequence[CompetitiveEventDecision],
    trajectories: Sequence[Mapping[str, object]],
    report_path: Path,
    selected_events_path: Path,
    rejections_path: Path,
    compatibility_path: Path,
    trajectory_path: Path,
    metrics_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(selected_event_rows(selected_events, decisions), selected_events_path)
    _write_csv_rows(rejection_rows(decisions), rejections_path)
    _write_csv_rows(compatibility_rows(selected_events), compatibility_path)
    _write_csv_rows(trajectories, trajectory_path)
    _write_csv_rows(metric_rows_from_report(report), metrics_path)
    dashboard_path.write_text(
        render_competitive_nucleus_selection_dashboard(report),
        encoding="utf-8",
    )
    certificate = build_competitive_nucleus_selection_certificate(
        report,
        output_names=ROOT_OUTPUT_NAMES,
    )
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return (
        report_path,
        selected_events_path,
        rejections_path,
        compatibility_path,
        trajectory_path,
        metrics_path,
        dashboard_path,
        certificate_path,
    )


def run_competitive_nucleus_selection(
    *,
    benchmark_file: Path,
    report_path: Path,
    selected_events_path: Path,
    rejections_path: Path,
    compatibility_path: Path,
    trajectory_path: Path,
    metrics_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path, Path]:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    feature_rows = contact_law_feature_rows(rows)
    events = nucleus_closure_events(rows, feature_rows)
    pre_report = build_folding_nucleus_closure_report(
        rows=rows,
        feature_rows=feature_rows,
        events=events,
        source_benchmark_file=benchmark_file,
    )
    pre_selected = accepted_events(events, threshold=PRE_COMPETITION_THRESHOLD)
    selected, decisions = select_competitive_events(rows, pre_selected)
    report = build_competitive_nucleus_selection_report(
        rows=rows,
        events=events,
        selected_events=selected,
        decisions=decisions,
        source_benchmark_file=benchmark_file,
        pre_competition_report=pre_report,
    )
    trajectories = trajectory_rows(rows, selected, decisions)
    return write_competitive_nucleus_selection_outputs(
        report=report,
        selected_events=selected,
        decisions=decisions,
        trajectories=trajectories,
        report_path=report_path,
        selected_events_path=selected_events_path,
        rejections_path=rejections_path,
        compatibility_path=compatibility_path,
        trajectory_path=trajectory_path,
        metrics_path=metrics_path,
        dashboard_path=dashboard_path,
        certificate_path=certificate_path,
    )


def _write_csv_rows(rows: Sequence[Mapping[str, object]], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as file:
        if not rows:
            return path
        fieldnames = list(rows[0])
        writer = csv.DictWriter(file, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    return path


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
        "pre_competition_event_count",
        "post_competition_selected_event_count",
        "event_reduction_ratio",
        "pre_false_nucleus_rate",
        "post_false_nucleus_rate",
        "pre_contact_cluster_precision",
        "post_contact_cluster_precision",
        "pre_long_range_contact_recall",
        "post_long_range_contact_recall",
        "nucleus_competition_law_survives",
        "mechanism_discovery_claim_allowed",
        "folding_problem_solved",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
    )


def _rule_cards() -> str:
    rules = (
        (
            "Competition Filters Event Floods",
            "The selector applies a fixed budget and ranks sequence-only closure events before native scoring.",
        ),
        (
            "Frustration Is A Real Rejection Surface",
            "Geometry, trap, overlap, and competition rejections are counted instead of hidden.",
        ),
        (
            "False Nucleus Rate Remains A Survival Gate",
            "Lower event count alone cannot claim a folding mechanism.",
        ),
        (
            "No Mechanism Discovery Claim",
            "Native labels are attached only after selection for audit metrics.",
        ),
    )
    return "".join(
        "<div class=\"rule\">"
        f"<h3>{_escape(title)}</h3><p>{_escape(body)}</p>"
        "</div>"
        for title, body in rules
    )


def _target_table(report: Mapping[str, object]) -> str:
    rows = (
        ("selected_event_count_target_met", TARGET_POST_SELECTED_EVENTS_MAX),
        ("false_nucleus_rate_target_met", TARGET_FALSE_NUCLEUS_RATE_MAX),
        ("contact_cluster_precision_target_met", TARGET_CONTACT_CLUSTER_PRECISION_MIN),
        ("long_range_contact_recall_target_met", TARGET_LONG_RANGE_RECALL_MIN),
    )
    body = "".join(
        "<tr>"
        f"<td>{_escape(name)}</td>"
        f"<td>{_escape(report.get(name))}</td>"
        f"<td>{_escape(target)}</td>"
        "</tr>"
        for name, target in rows
    )
    return (
        "<section><h2>Survival Targets</h2>"
        "<table><thead><tr><th>target</th><th>met</th><th>threshold</th>"
        "</tr></thead><tbody>"
        + body
        + "</tbody></table></section>"
    )


def _metric_table(report: Mapping[str, object]) -> str:
    body = []
    for row in metric_rows_from_report(report):
        body.append(
            "<tr>"
            f"<td>{_escape(row['source_accession'])}</td>"
            f"<td>{_escape(row['post_competition_selected_event_count'])}</td>"
            f"<td>{_escape(row['post_false_nucleus_rate'])}</td>"
            f"<td>{_escape(row['post_contact_cluster_precision'])}</td>"
            f"<td>{_escape(row['post_long_range_contact_recall'])}</td>"
            f"<td>{_escape(row['event_reduction_ratio'])}</td>"
            "</tr>"
        )
    return (
        "<section><h2>Row Metrics</h2>"
        "<table><thead><tr>"
        "<th>source</th><th>selected</th><th>false rate</th>"
        "<th>cluster precision</th><th>long-range recall</th><th>reduction</th>"
        "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></section>"
    )


def render_competitive_nucleus_selection_dashboard(
    report: Mapping[str, object],
) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Competitive Nucleus Selection</title>
  <style>
    body {{
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f4f6f3;
      color: #202522;
    }}
    header {{
      padding: 34px;
      background: #26332f;
      color: #f7f8f4;
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
    .metrics, .rules {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
      gap: 12px;
    }}
    .metric, .rule {{
      background: #ffffff;
      border: 1px solid #d3ddd5;
      border-radius: 6px;
      padding: 14px;
    }}
    .metric span {{
      display: block;
      color: #5b6660;
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
      border: 1px solid #d3ddd5;
    }}
    th, td {{
      padding: 9px 10px;
      border-bottom: 1px solid #e2e8e3;
      text-align: left;
      font-size: 13px;
      vertical-align: top;
      overflow-wrap: anywhere;
    }}
  </style>
</head>
<body>
  <header>
    <h1>Competitive Nucleus Selection</h1>
    <p>Sequence-only closure candidates are ranked, filtered, and scored against native coordinates only after selection.</p>
  </header>
  <main>
    <section class="metrics">{_metric_cards(report)}</section>
    <section><h2>Boundary Rules</h2><div class="rules">{_rule_cards()}</div></section>
    {_target_table(report)}
    {_metric_table(report)}
    <section>
      <h2>Claim Boundary</h2>
      <p>{_escape(report.get("boundary_statement", ""))}</p>
    </section>
  </main>
</body>
</html>
"""
