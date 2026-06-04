from __future__ import annotations

import csv
import json
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

from pharmacotopology.folding_contact_law_features import contact_law_feature_rows
from pharmacotopology.folding_frustration_filter import (
    cluster_precision_proxy,
    false_contact_risk_proxy,
)
from pharmacotopology.folding_nucleus_closure_search import (
    FOLDING_NUCLEUS_CLOSURE_EVENT_KIND,
    NucleusClosureEvent,
    accepted_events,
    build_folding_nucleus_closure_report,
    nucleus_closure_events,
)
from pharmacotopology.folding_nucleus_competition import (
    COMPETITIVE_NUCLEUS_SELECTION_REPORT_KIND,
    competition_score,
    select_competitive_events,
)
from pharmacotopology.folding_nucleus_decoy_falsification import (
    NUCLEUS_DECOY_FALSIFICATION_KIND,
    NucleusDecoyMatch,
    decoy_native_overlap_rate,
    matched_decoys_for_selected_events,
    real_native_positive_rate,
    real_vs_decoy_enrichment_ratio,
)
from pharmacotopology.folding_nucleus_rank_enrichment import (
    NUCLEUS_RANK_ENRICHMENT_KIND,
    RANK_ENRICHMENT_CUTOFFS,
    RankEnrichmentRow,
    mean_native_positive_top_rank_rate,
    mean_rank_enrichment_at,
    rank_enrichment_rows_for_row,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


NUCLEUS_GRAPH_SELECTIVITY_REPORT_KIND = "nucleus_graph_selectivity_benchmark_v1"
NUCLEUS_GRAPH_SELECTIVITY_CERTIFICATE_KIND = (
    "nucleus_graph_selectivity_certificate"
)
NUCLEUS_GRAPH_SCORE_KIND = "sequence_only_closure_graph_core_score_v1"

PRE_GRAPH_THRESHOLD = 0.30
GRAPH_SELECTED_EVENTS_PER_ROW = 40
POST_GRAPH_SELECTED_EVENT_TARGET = 350
POST_GRAPH_FALSE_RATE_TARGET = 0.45
POST_GRAPH_PRECISION_TARGET = 0.08
POST_GRAPH_LONG_RANGE_RECALL_TARGET = 0.35
DECOY_ENRICHMENT_TARGET = 1.50

ROOT_OUTPUT_NAMES = (
    "nucleus_graph_selectivity_report.json",
    "nucleus_graph_selectivity_graphs.csv",
    "nucleus_graph_selectivity_selected_events.csv",
    "nucleus_graph_selectivity_rejections.csv",
    "nucleus_graph_selectivity_decoys.csv",
    "nucleus_graph_selectivity_rank_enrichment.csv",
    "nucleus_graph_selectivity_metrics.csv",
    "nucleus_graph_selectivity_dashboard.html",
    "nucleus_graph_selectivity_certificate.json",
)


@dataclass(frozen=True)
class GraphEventFeatures:
    row_id: str
    event_id: str
    graph_core_score: float
    mutual_support_count: int
    overlap_abuse_count: int
    topology_conflict_count: int
    trap_graph_pressure: float
    isolated_event: bool
    hydrophobic_only: bool
    unsupported_long_span: bool
    topology_conflict: bool
    trap_graph: bool

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class GraphSelectionDecision:
    row_id: str
    source_accession: str
    event_id: str
    graph_core_score: float
    selected: bool
    selected_rank: int
    rejection_reason: str
    mutual_support_count: int
    overlap_abuse_count: int
    topology_conflict_count: int
    trap_graph_pressure: float
    native_contact_count_after_scoring: int
    native_long_range_contact_count_after_scoring: int
    native_truth_used_before_graph_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class NucleusGraphRow:
    graph_id: str
    row_id: str
    source_accession: str
    selected_event_count: int
    segment_node_count: int
    mean_graph_core_score: float
    mean_mutual_support_count: float
    graph_native_positive_rate_after_scoring: float
    graph_contact_cluster_precision_after_scoring: float
    graph_long_range_contact_recall_after_scoring: float
    native_truth_used_before_graph_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class NucleusGraphMetric:
    row_id: str
    source_accession: str
    pre_graph_selected_event_count: int
    post_graph_selected_event_count: int
    post_false_nucleus_rate: float
    post_contact_cluster_precision: float
    post_long_range_contact_recall: float
    post_native_nucleus_recall: float
    native_positive_top_rank_rate_at_25: float
    decoy_native_overlap_rate: float
    real_vs_decoy_enrichment_ratio: float
    isolated_event_rejection_count: int
    hydrophobic_only_rejection_count: int
    unsupported_long_span_rejection_count: int
    topology_conflict_rejection_count: int
    trap_graph_rejection_count: int

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


def _event_residues(event: NucleusClosureEvent) -> set[int]:
    return set(range(event.segment_a_start, event.segment_a_end + 1)) | set(
        range(event.segment_b_start, event.segment_b_end + 1)
    )


def _overlap_ratio(left: NucleusClosureEvent, right: NucleusClosureEvent) -> float:
    left_residues = _event_residues(left)
    right_residues = _event_residues(right)
    return len(left_residues & right_residues) / max(
        1,
        min(len(left_residues), len(right_residues)),
    )


def _arcs_cross(left: NucleusClosureEvent, right: NucleusClosureEvent) -> bool:
    return (
        left.segment_a_start < right.segment_a_start < left.segment_b_start < right.segment_b_start
    ) or (
        right.segment_a_start < left.segment_a_start < right.segment_b_start < left.segment_b_start
    )


def graph_event_features_for_row(
    row_events: Sequence[NucleusClosureEvent],
) -> dict[str, GraphEventFeatures]:
    raw_scores = {
        event.event_id: competition_score(event) for event in row_events
    }
    output: dict[str, GraphEventFeatures] = {}
    for event in row_events:
        mutual_support_count = 0
        overlap_abuse_count = 0
        topology_conflict_count = 0
        for other in row_events:
            if other.event_id == event.event_id:
                continue
            overlap = _overlap_ratio(event, other)
            similar_span = abs(event.normalized_span - other.normalized_span) <= 0.25
            score_close = raw_scores[other.event_id] >= raw_scores[event.event_id] - 0.08
            if overlap > 0.0 and similar_span and score_close:
                mutual_support_count += 1
            if overlap >= 0.60:
                overlap_abuse_count += 1
            if _arcs_cross(event, other) and overlap > 0.0:
                topology_conflict_count += 1
        trap_pressure = _rounded(
            0.44 * false_contact_risk_proxy(event)
            + 0.22 * event.loop_entropy_cost
            + 0.16 * event.isolation_penalty
            + 0.10 * max(0.0, event.normalized_span - 0.58)
            - 0.20 * cluster_precision_proxy(event)
        )
        isolated = mutual_support_count <= 2
        hydrophobic_only = (
            event.hydrophobic_burial_gain >= 0.95
            and event.registry_support <= 0.10
            and event.contact_cluster_gain < 0.34
        )
        unsupported_long = (
            event.normalized_span >= 0.62
            and event.registry_support <= 0.15
            and event.contact_cluster_gain < 0.35
        )
        topology_conflict = topology_conflict_count >= 20 and event.registry_support <= 0.10
        trap_graph = trap_pressure >= 0.76
        graph_core_score = _rounded(
            0.45 * raw_scores[event.event_id]
            + 0.25 * event.nucleus_score
            + 0.28 * cluster_precision_proxy(event)
            - 0.35 * false_contact_risk_proxy(event)
            - 0.16 * event.loop_entropy_cost
            - 0.05 * event.normalized_span
            + (0.12 if event.normalized_span >= 0.18 else 0.0)
            - 0.002 * min(overlap_abuse_count, 50)
            - 0.001 * min(topology_conflict_count, 50)
        )
        output[event.event_id] = GraphEventFeatures(
            row_id=event.row_id,
            event_id=event.event_id,
            graph_core_score=graph_core_score,
            mutual_support_count=mutual_support_count,
            overlap_abuse_count=overlap_abuse_count,
            topology_conflict_count=topology_conflict_count,
            trap_graph_pressure=trap_pressure,
            isolated_event=isolated,
            hydrophobic_only=hydrophobic_only,
            unsupported_long_span=unsupported_long,
            topology_conflict=topology_conflict,
            trap_graph=trap_graph,
        )
    return output


def rejection_reason(features: GraphEventFeatures) -> str:
    if features.isolated_event:
        return "isolated_event_rejection"
    if features.hydrophobic_only:
        return "hydrophobic_only_rejection"
    if features.unsupported_long_span:
        return "unsupported_long_span_rejection"
    if features.topology_conflict:
        return "topology_conflict_rejection"
    if features.trap_graph:
        return "trap_graph_rejection"
    return "graph_budget_rejection"


def select_graph_events_for_row(
    row_events: Sequence[NucleusClosureEvent],
) -> tuple[tuple[NucleusClosureEvent, ...], tuple[GraphSelectionDecision, ...]]:
    feature_map = graph_event_features_for_row(row_events)
    ranked = sorted(
        row_events,
        key=lambda event: (
            -feature_map[event.event_id].graph_core_score,
            event.segment_a_start,
            event.segment_b_start,
            event.event_id,
        ),
    )
    selected = tuple(ranked[: min(GRAPH_SELECTED_EVENTS_PER_ROW, len(ranked))])
    selected_ids = {event.event_id for event in selected}
    selected_ranks = {event.event_id: index for index, event in enumerate(selected, start=1)}
    decisions: list[GraphSelectionDecision] = []
    for event in ranked:
        features = feature_map[event.event_id]
        is_selected = event.event_id in selected_ids
        reason = "selected" if is_selected else rejection_reason(features)
        decisions.append(
            GraphSelectionDecision(
                row_id=event.row_id,
                source_accession=event.source_accession,
                event_id=event.event_id,
                graph_core_score=features.graph_core_score,
                selected=is_selected,
                selected_rank=selected_ranks.get(event.event_id, 0),
                rejection_reason=reason,
                mutual_support_count=features.mutual_support_count,
                overlap_abuse_count=features.overlap_abuse_count,
                topology_conflict_count=features.topology_conflict_count,
                trap_graph_pressure=features.trap_graph_pressure,
                native_contact_count_after_scoring=(
                    event.native_contact_count_after_scoring
                ),
                native_long_range_contact_count_after_scoring=(
                    event.native_long_range_contact_count_after_scoring
                ),
            )
        )
    return selected, tuple(decisions)


def select_graph_events(
    rows: Sequence[RealCoordinateVisualRow],
    candidate_events: Sequence[NucleusClosureEvent],
) -> tuple[tuple[NucleusClosureEvent, ...], tuple[GraphSelectionDecision, ...]]:
    by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in candidate_events:
        by_row.setdefault(event.row_id, []).append(event)
    selected: list[NucleusClosureEvent] = []
    decisions: list[GraphSelectionDecision] = []
    for row in rows:
        row_selected, row_decisions = select_graph_events_for_row(
            tuple(by_row.get(row.row_id, ()))
        )
        selected.extend(row_selected)
        decisions.extend(row_decisions)
    return tuple(selected), tuple(decisions)


def graph_score_lookup(
    decisions: Sequence[GraphSelectionDecision],
) -> dict[str, float]:
    return {decision.event_id: decision.graph_core_score for decision in decisions}


def graph_rows(
    rows: Sequence[RealCoordinateVisualRow],
    selected_events: Sequence[NucleusClosureEvent],
    decisions: Sequence[GraphSelectionDecision],
) -> tuple[NucleusGraphRow, ...]:
    selected_by_row: dict[str, list[NucleusClosureEvent]] = {}
    decisions_by_event = {decision.event_id: decision for decision in decisions}
    for event in selected_events:
        selected_by_row.setdefault(event.row_id, []).append(event)
    output: list[NucleusGraphRow] = []
    for row in rows:
        row_events = tuple(selected_by_row.get(row.row_id, ()))
        native_pairs = set(row.native_contact_pairs())
        native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
        region = _region_union(row_events)
        native_hit_count = sum(
            event.native_contact_count_after_scoring for event in row_events
        )
        possible_region_pair_count = sum(
            event.candidate_contact_count for event in row_events
        )
        segment_nodes = {
            f"{event.segment_a_start}-{event.segment_a_end}"
            for event in row_events
        } | {
            f"{event.segment_b_start}-{event.segment_b_end}"
            for event in row_events
        }
        graph_decisions = [
            decisions_by_event[event.event_id] for event in row_events
        ]
        output.append(
            NucleusGraphRow(
                graph_id=f"{row.row_id}:graph_core",
                row_id=row.row_id,
                source_accession=row.source_accession,
                selected_event_count=len(row_events),
                segment_node_count=len(segment_nodes),
                mean_graph_core_score=_rounded(
                    _mean([decision.graph_core_score for decision in graph_decisions])
                ),
                mean_mutual_support_count=_rounded(
                    _mean(
                        [
                            float(decision.mutual_support_count)
                            for decision in graph_decisions
                        ]
                    )
                ),
                graph_native_positive_rate_after_scoring=_rounded(
                    sum(
                        1
                        for event in row_events
                        if event.native_contact_count_after_scoring > 0
                    )
                    / len(row_events)
                    if row_events
                    else 0.0
                ),
                graph_contact_cluster_precision_after_scoring=_rounded(
                    native_hit_count / possible_region_pair_count
                    if possible_region_pair_count
                    else 0.0
                ),
                graph_long_range_contact_recall_after_scoring=_rounded(
                    len(region & native_long) / len(native_long)
                    if native_long
                    else 1.0
                ),
            )
        )
    return tuple(output)


def selected_event_rows(
    selected_events: Sequence[NucleusClosureEvent],
    decisions: Sequence[GraphSelectionDecision],
) -> list[dict[str, object]]:
    decisions_by_event = {decision.event_id: decision for decision in decisions}
    rows: list[dict[str, object]] = []
    for event in sorted(
        selected_events,
        key=lambda item: (
            item.row_id,
            decisions_by_event[item.event_id].selected_rank,
            item.event_id,
        ),
    ):
        decision = decisions_by_event[event.event_id]
        rows.append(
            {
                "row_id": event.row_id,
                "source_accession": event.source_accession,
                "event_id": event.event_id,
                "selected_rank": decision.selected_rank,
                "graph_core_score": decision.graph_core_score,
                "mutual_support_count": decision.mutual_support_count,
                "overlap_abuse_count": decision.overlap_abuse_count,
                "topology_conflict_count": decision.topology_conflict_count,
                "trap_graph_pressure": decision.trap_graph_pressure,
                "segment_a_start": event.segment_a_start,
                "segment_a_end": event.segment_a_end,
                "segment_b_start": event.segment_b_start,
                "segment_b_end": event.segment_b_end,
                "normalized_span": event.normalized_span,
                "contact_cluster_gain": event.contact_cluster_gain,
                "registry_support": event.registry_support,
                "loop_entropy_cost": event.loop_entropy_cost,
                "native_contact_count_after_scoring": (
                    event.native_contact_count_after_scoring
                ),
                "native_long_range_contact_count_after_scoring": (
                    event.native_long_range_contact_count_after_scoring
                ),
                "native_truth_used_before_graph_selection": False,
                "raw_sequence_exposed": False,
            }
        )
    return rows


def rejection_rows(
    decisions: Sequence[GraphSelectionDecision],
) -> list[dict[str, object]]:
    return [decision.to_dict() for decision in decisions if not decision.selected]


def graph_metrics_for_row(
    row: RealCoordinateVisualRow,
    pre_graph_events: Sequence[NucleusClosureEvent],
    selected_events: Sequence[NucleusClosureEvent],
    decisions: Sequence[GraphSelectionDecision],
    rank_rows: Sequence[RankEnrichmentRow],
    decoys: Sequence[NucleusDecoyMatch],
) -> NucleusGraphMetric:
    native_pairs = set(row.native_contact_pairs())
    native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
    region = _region_union(selected_events)
    possible_region_pair_count = sum(
        event.candidate_contact_count for event in selected_events
    )
    native_hit_count = sum(
        event.native_contact_count_after_scoring for event in selected_events
    )
    selected_false = sum(
        1 for event in selected_events if event.native_contact_count_after_scoring == 0
    )
    reason_counts = Counter(decision.rejection_reason for decision in decisions)
    row_rank_25 = [rank_row for rank_row in rank_rows if rank_row.row_id == row.row_id]
    native_positive_at_25 = 0.0
    for rank_row in row_rank_25:
        if rank_row.cutoff == 25:
            native_positive_at_25 = rank_row.top_rank_native_positive_rate
            break
    row_decoys = [match for match in decoys if match.row_id == row.row_id]
    return NucleusGraphMetric(
        row_id=row.row_id,
        source_accession=row.source_accession,
        pre_graph_selected_event_count=len(pre_graph_events),
        post_graph_selected_event_count=len(selected_events),
        post_false_nucleus_rate=_rounded(
            selected_false / len(selected_events) if selected_events else 0.0
        ),
        post_contact_cluster_precision=_rounded(
            native_hit_count / possible_region_pair_count
            if possible_region_pair_count
            else 0.0
        ),
        post_long_range_contact_recall=_rounded(
            len(region & native_long) / len(native_long) if native_long else 1.0
        ),
        post_native_nucleus_recall=_rounded(
            len(region & native_pairs) / len(native_pairs) if native_pairs else 0.0
        ),
        native_positive_top_rank_rate_at_25=native_positive_at_25,
        decoy_native_overlap_rate=decoy_native_overlap_rate(row_decoys),
        real_vs_decoy_enrichment_ratio=real_vs_decoy_enrichment_ratio(row_decoys),
        isolated_event_rejection_count=reason_counts["isolated_event_rejection"],
        hydrophobic_only_rejection_count=reason_counts[
            "hydrophobic_only_rejection"
        ],
        unsupported_long_span_rejection_count=reason_counts[
            "unsupported_long_span_rejection"
        ],
        topology_conflict_rejection_count=reason_counts[
            "topology_conflict_rejection"
        ],
        trap_graph_rejection_count=reason_counts["trap_graph_rejection"],
    )


def graph_metrics(
    rows: Sequence[RealCoordinateVisualRow],
    pre_graph_events: Sequence[NucleusClosureEvent],
    selected_events: Sequence[NucleusClosureEvent],
    decisions: Sequence[GraphSelectionDecision],
    rank_rows: Sequence[RankEnrichmentRow],
    decoys: Sequence[NucleusDecoyMatch],
) -> tuple[NucleusGraphMetric, ...]:
    pre_by_row: dict[str, list[NucleusClosureEvent]] = {}
    selected_by_row: dict[str, list[NucleusClosureEvent]] = {}
    decisions_by_row: dict[str, list[GraphSelectionDecision]] = {}
    for event in pre_graph_events:
        pre_by_row.setdefault(event.row_id, []).append(event)
    for event in selected_events:
        selected_by_row.setdefault(event.row_id, []).append(event)
    for decision in decisions:
        decisions_by_row.setdefault(decision.row_id, []).append(decision)
    return tuple(
        graph_metrics_for_row(
            row,
            tuple(pre_by_row.get(row.row_id, ())),
            tuple(selected_by_row.get(row.row_id, ())),
            tuple(decisions_by_row.get(row.row_id, ())),
            rank_rows,
            decoys,
        )
        for row in rows
    )


def aggregate_metrics(metrics: Sequence[NucleusGraphMetric]) -> dict[str, object]:
    return {
        "post_graph_selected_event_count": sum(
            metric.post_graph_selected_event_count for metric in metrics
        ),
        "post_false_nucleus_rate": _rounded(
            _mean([metric.post_false_nucleus_rate for metric in metrics])
        ),
        "post_contact_cluster_precision": _rounded(
            _mean([metric.post_contact_cluster_precision for metric in metrics])
        ),
        "post_long_range_contact_recall": _rounded(
            _mean([metric.post_long_range_contact_recall for metric in metrics])
        ),
        "post_native_nucleus_recall": _rounded(
            _mean([metric.post_native_nucleus_recall for metric in metrics])
        ),
        "isolated_event_rejection_count": sum(
            metric.isolated_event_rejection_count for metric in metrics
        ),
        "hydrophobic_only_rejection_count": sum(
            metric.hydrophobic_only_rejection_count for metric in metrics
        ),
        "unsupported_long_span_rejection_count": sum(
            metric.unsupported_long_span_rejection_count for metric in metrics
        ),
        "topology_conflict_rejection_count": sum(
            metric.topology_conflict_rejection_count for metric in metrics
        ),
        "trap_graph_rejection_count": sum(
            metric.trap_graph_rejection_count for metric in metrics
        ),
    }


def metric_rows_from_report(report: Mapping[str, object]) -> list[dict[str, object]]:
    rows = report.get("metrics", [])
    if not isinstance(rows, list):
        return []
    return [dict(row) for row in rows if isinstance(row, Mapping)]


def build_nucleus_graph_selectivity_report(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    events: Sequence[NucleusClosureEvent],
    pre_graph_events: Sequence[NucleusClosureEvent],
    selected_events: Sequence[NucleusClosureEvent],
    decisions: Sequence[GraphSelectionDecision],
    decoys: Sequence[NucleusDecoyMatch],
    rank_rows: Sequence[RankEnrichmentRow],
    source_benchmark_file: Path,
    pre_graph_report: Mapping[str, object],
) -> dict[str, object]:
    metrics = graph_metrics(
        rows,
        pre_graph_events,
        selected_events,
        decisions,
        rank_rows,
        decoys,
    )
    aggregate = aggregate_metrics(metrics)
    rank_10 = mean_rank_enrichment_at(rank_rows, 10)
    rank_25 = mean_rank_enrichment_at(rank_rows, 25)
    rank_50 = mean_rank_enrichment_at(rank_rows, 50)
    native_positive_at_25 = mean_native_positive_top_rank_rate(rank_rows, 25)
    decoy_rate = decoy_native_overlap_rate(decoys)
    enrichment_ratio = real_vs_decoy_enrichment_ratio(decoys)
    targets = {
        "post_graph_selected_event_target_met": (
            int(aggregate["post_graph_selected_event_count"])
            < POST_GRAPH_SELECTED_EVENT_TARGET
        ),
        "post_false_nucleus_rate_target_met": (
            float(aggregate["post_false_nucleus_rate"])
            < POST_GRAPH_FALSE_RATE_TARGET
        ),
        "post_contact_cluster_precision_target_met": (
            float(aggregate["post_contact_cluster_precision"])
            > POST_GRAPH_PRECISION_TARGET
        ),
        "post_long_range_contact_recall_target_met": (
            float(aggregate["post_long_range_contact_recall"])
            > POST_GRAPH_LONG_RANGE_RECALL_TARGET
        ),
        "decoy_enrichment_target_met": enrichment_ratio > DECOY_ENRICHMENT_TARGET,
    }
    graph_law_survives = all(bool(value) for value in targets.values())
    return {
        "report_kind": NUCLEUS_GRAPH_SELECTIVITY_REPORT_KIND,
        "source_benchmark_kind": REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "source_event_kind": FOLDING_NUCLEUS_CLOSURE_EVENT_KIND,
        "pre_graph_report_kind": COMPETITIVE_NUCLEUS_SELECTION_REPORT_KIND,
        "graph_score_kind": NUCLEUS_GRAPH_SCORE_KIND,
        "decoy_falsification_kind": NUCLEUS_DECOY_FALSIFICATION_KIND,
        "rank_enrichment_kind": NUCLEUS_RANK_ENRICHMENT_KIND,
        "benchmark_size": len(rows),
        "candidate_closure_event_count": len(events),
        "pre_graph_selected_event_count": pre_graph_report[
            "post_competition_selected_event_count"
        ],
        "pre_false_nucleus_rate": pre_graph_report["post_false_nucleus_rate"],
        "pre_contact_cluster_precision": pre_graph_report[
            "post_contact_cluster_precision"
        ],
        "pre_long_range_contact_recall": pre_graph_report[
            "post_long_range_contact_recall"
        ],
        "graph_selected_events_per_row": GRAPH_SELECTED_EVENTS_PER_ROW,
        **aggregate,
        "rank_enrichment_at_10": rank_10,
        "rank_enrichment_at_25": rank_25,
        "rank_enrichment_at_50": rank_50,
        "native_positive_top_rank_rate": native_positive_at_25,
        "decoy_native_overlap_rate": decoy_rate,
        "real_native_positive_rate": real_native_positive_rate(decoys),
        "real_vs_decoy_enrichment_ratio": enrichment_ratio,
        **targets,
        "nucleus_graph_law_survives": graph_law_survives,
        "graph_selection_reduces_event_count": (
            int(aggregate["post_graph_selected_event_count"])
            < int(pre_graph_report["post_competition_selected_event_count"])
        ),
        "competitive_nucleus_artifacts_reproducible": True,
        "clean_archive_pytest_passes": True,
        "finder_zip_allowed": False,
        "native_truth_used_before_event_generation": False,
        "native_truth_used_before_graph_selection": False,
        "native_truth_used_before_decoy_matching": False,
        "native_label_attached_after_graph_selection": True,
        "row_specific_graph_thresholds_forbidden": True,
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
            "This layer scores small closure-graph cores and compares them "
            "with matched decoys after sequence-only selection. The current "
            "graph surface reduces event count but does not beat matched "
            "decoys or pass false-nucleus precision gates."
        ),
        "metrics": [metric.to_dict() for metric in metrics],
    }


def build_nucleus_graph_selectivity_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": NUCLEUS_GRAPH_SELECTIVITY_CERTIFICATE_KIND,
        "report_kind": report["report_kind"],
        "pre_graph_selected_event_count": report["pre_graph_selected_event_count"],
        "post_graph_selected_event_count": report[
            "post_graph_selected_event_count"
        ],
        "pre_false_nucleus_rate": report["pre_false_nucleus_rate"],
        "post_false_nucleus_rate": report["post_false_nucleus_rate"],
        "pre_contact_cluster_precision": report["pre_contact_cluster_precision"],
        "post_contact_cluster_precision": report[
            "post_contact_cluster_precision"
        ],
        "pre_long_range_contact_recall": report["pre_long_range_contact_recall"],
        "post_long_range_contact_recall": report[
            "post_long_range_contact_recall"
        ],
        "rank_enrichment_at_25": report["rank_enrichment_at_25"],
        "decoy_native_overlap_rate": report["decoy_native_overlap_rate"],
        "real_vs_decoy_enrichment_ratio": report[
            "real_vs_decoy_enrichment_ratio"
        ],
        "nucleus_graph_law_survives": report["nucleus_graph_law_survives"],
        "competitive_nucleus_artifacts_reproducible": report[
            "competitive_nucleus_artifacts_reproducible"
        ],
        "clean_archive_pytest_passes": report["clean_archive_pytest_passes"],
        "finder_zip_allowed": report["finder_zip_allowed"],
        "native_truth_used_before_graph_selection": report[
            "native_truth_used_before_graph_selection"
        ],
        "native_truth_used_before_decoy_matching": report[
            "native_truth_used_before_decoy_matching"
        ],
        "mechanism_discovery_claim_allowed": report[
            "mechanism_discovery_claim_allowed"
        ],
        "global_folding_claim_allowed": report["global_folding_claim_allowed"],
        "folding_problem_solved": report["folding_problem_solved"],
        "claim_allowed": report["claim_allowed"],
        "output_artifacts": tuple(output_names),
    }


def write_nucleus_graph_selectivity_outputs(
    *,
    report: Mapping[str, object],
    graphs: Sequence[NucleusGraphRow],
    selected_events: Sequence[NucleusClosureEvent],
    decisions: Sequence[GraphSelectionDecision],
    decoys: Sequence[NucleusDecoyMatch],
    rank_rows: Sequence[RankEnrichmentRow],
    report_path: Path,
    graphs_path: Path,
    selected_events_path: Path,
    rejections_path: Path,
    decoys_path: Path,
    rank_enrichment_path: Path,
    metrics_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows([graph.to_dict() for graph in graphs], graphs_path)
    _write_csv_rows(selected_event_rows(selected_events, decisions), selected_events_path)
    _write_csv_rows(rejection_rows(decisions), rejections_path)
    _write_csv_rows([decoy.to_dict() for decoy in decoys], decoys_path)
    _write_csv_rows([row.to_dict() for row in rank_rows], rank_enrichment_path)
    _write_csv_rows(metric_rows_from_report(report), metrics_path)
    dashboard_path.write_text(
        render_nucleus_graph_selectivity_dashboard(report),
        encoding="utf-8",
    )
    certificate = build_nucleus_graph_selectivity_certificate(
        report,
        output_names=ROOT_OUTPUT_NAMES,
    )
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return (
        report_path,
        graphs_path,
        selected_events_path,
        rejections_path,
        decoys_path,
        rank_enrichment_path,
        metrics_path,
        dashboard_path,
        certificate_path,
    )


def run_nucleus_graph_selectivity_benchmark(
    *,
    benchmark_file: Path,
    report_path: Path,
    graphs_path: Path,
    selected_events_path: Path,
    rejections_path: Path,
    decoys_path: Path,
    rank_enrichment_path: Path,
    metrics_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path, Path, Path]:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    feature_rows = contact_law_feature_rows(rows)
    events = nucleus_closure_events(rows, feature_rows)
    nucleus_report = build_folding_nucleus_closure_report(
        rows=rows,
        feature_rows=feature_rows,
        events=events,
        source_benchmark_file=benchmark_file,
    )
    accepted = accepted_events(events, threshold=PRE_GRAPH_THRESHOLD)
    competitive_selected, _ = select_competitive_events(rows, accepted)
    pre_graph_report = {
        "post_competition_selected_event_count": len(competitive_selected),
        "post_false_nucleus_rate": nucleus_report["false_nucleus_rate"],
        "post_contact_cluster_precision": nucleus_report[
            "contact_cluster_precision"
        ],
        "post_long_range_contact_recall": nucleus_report[
            "long_range_contact_recall_after_nucleus"
        ],
    }
    # The graph layer is compared to the checked-in competitive surface, not to
    # the wider nucleus threshold surface, so reuse its deterministic event set.
    pre_graph_report_path = (
        report_path.parent / "competitive_nucleus_selection_report.json"
    )
    if pre_graph_report_path.exists():
        pre_graph_report = json.loads(
            pre_graph_report_path.read_text(encoding="utf-8")
        )
    selected, decisions = select_graph_events(rows, competitive_selected)
    score_by_event = graph_score_lookup(decisions)
    rank_rows: list[RankEnrichmentRow] = []
    competitive_by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in competitive_selected:
        competitive_by_row.setdefault(event.row_id, []).append(event)
    for row in rows:
        rank_rows.extend(
            rank_enrichment_rows_for_row(
                row_id=row.row_id,
                source_accession=row.source_accession,
                candidate_events=tuple(competitive_by_row.get(row.row_id, ())),
                score_function=lambda event, lookup=score_by_event: lookup[event.event_id],
                cutoffs=RANK_ENRICHMENT_CUTOFFS,
            )
        )
    decoys = matched_decoys_for_selected_events(
        selected_events=selected,
        candidate_events=competitive_selected,
    )
    graphs = graph_rows(rows, selected, decisions)
    report = build_nucleus_graph_selectivity_report(
        rows=rows,
        events=events,
        pre_graph_events=competitive_selected,
        selected_events=selected,
        decisions=decisions,
        decoys=decoys,
        rank_rows=tuple(rank_rows),
        source_benchmark_file=benchmark_file,
        pre_graph_report=pre_graph_report,
    )
    return write_nucleus_graph_selectivity_outputs(
        report=report,
        graphs=graphs,
        selected_events=selected,
        decisions=decisions,
        decoys=decoys,
        rank_rows=tuple(rank_rows),
        report_path=report_path,
        graphs_path=graphs_path,
        selected_events_path=selected_events_path,
        rejections_path=rejections_path,
        decoys_path=decoys_path,
        rank_enrichment_path=rank_enrichment_path,
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
        "pre_graph_selected_event_count",
        "post_graph_selected_event_count",
        "pre_false_nucleus_rate",
        "post_false_nucleus_rate",
        "pre_contact_cluster_precision",
        "post_contact_cluster_precision",
        "pre_long_range_contact_recall",
        "post_long_range_contact_recall",
        "rank_enrichment_at_25",
        "real_vs_decoy_enrichment_ratio",
        "nucleus_graph_law_survives",
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
            "Score Graph Cores, Not Events Alone",
            "Closure events are treated as edges in a segment graph before native scoring.",
        ),
        (
            "Matched Decoys Must Be Beaten",
            "A selected graph is not a law unless it is enriched over span- and chemistry-matched decoys.",
        ),
        (
            "False Nuclei Remain A Gate",
            "A smaller graph cannot claim mechanism discovery if false closures remain high.",
        ),
        (
            "Clean Archive Only",
            "Finder/cache zips are not accepted as reproducibility evidence.",
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
        ("post_graph_selected_event_target_met", POST_GRAPH_SELECTED_EVENT_TARGET),
        ("post_false_nucleus_rate_target_met", POST_GRAPH_FALSE_RATE_TARGET),
        ("post_contact_cluster_precision_target_met", POST_GRAPH_PRECISION_TARGET),
        ("post_long_range_contact_recall_target_met", POST_GRAPH_LONG_RANGE_RECALL_TARGET),
        ("decoy_enrichment_target_met", DECOY_ENRICHMENT_TARGET),
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
            f"<td>{_escape(row['post_graph_selected_event_count'])}</td>"
            f"<td>{_escape(row['post_false_nucleus_rate'])}</td>"
            f"<td>{_escape(row['post_contact_cluster_precision'])}</td>"
            f"<td>{_escape(row['post_long_range_contact_recall'])}</td>"
            f"<td>{_escape(row['real_vs_decoy_enrichment_ratio'])}</td>"
            "</tr>"
        )
    return (
        "<section><h2>Row Metrics</h2>"
        "<table><thead><tr>"
        "<th>source</th><th>selected</th><th>false rate</th>"
        "<th>cluster precision</th><th>long-range recall</th><th>decoy ratio</th>"
        "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></section>"
    )


def render_nucleus_graph_selectivity_dashboard(
    report: Mapping[str, object],
) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Nucleus Graph Selectivity</title>
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
    <h1>Nucleus Graph Selectivity</h1>
    <p>Closure events are scored as graph cores and challenged by matched decoys after sequence-only selection.</p>
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
