from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

from pharmacotopology.folding_contact_law_features import (
    ContactLawFeatureRow,
    contact_law_feature_rows,
    feature_rows_by_row_id,
)
from pharmacotopology.folding_native_contact_eval import evaluate_contact_prediction
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


FOLDING_NUCLEUS_CLOSURE_REPORT_KIND = "folding_nucleus_closure_search_v1"
FOLDING_NUCLEUS_CLOSURE_CERTIFICATE_KIND = (
    "folding_nucleus_closure_search_certificate"
)
FOLDING_NUCLEUS_CLOSURE_EVENT_KIND = "sequence_only_segment_closure_event_v1"

SEGMENT_LENGTH = 8
SEGMENT_STRIDE = 4
MIN_SEGMENT_GAP = 4
PAIR_LEVEL_REFERENCE_MODEL = "pair_plus_entropy_score"
PAIR_LEVEL_REFERENCE_THRESHOLD = 0.44
FRUSTRATION_REJECTION_LIMIT = 0.75
NUCLEUS_THRESHOLD_MIN = 0.30
NUCLEUS_THRESHOLDS = tuple(round(index / 100, 2) for index in range(30, 101))

ROOT_OUTPUT_NAMES = (
    "folding_nucleus_closure_report.json",
    "folding_nucleus_closure_events.csv",
    "folding_nucleus_closure_trajectory.csv",
    "folding_nucleus_closure_metrics.csv",
    "folding_nucleus_closure_failures.csv",
    "folding_nucleus_closure_dashboard.html",
    "folding_nucleus_closure_certificate.json",
)


@dataclass(frozen=True)
class NucleusClosureEvent:
    row_id: str
    source_accession: str
    sequence_hash: str
    sequence_length: int
    event_id: str
    segment_a_start: int
    segment_a_end: int
    segment_b_start: int
    segment_b_end: int
    sequence_span: int
    normalized_span: float
    candidate_contact_count: int
    contact_cluster_gain: float
    secondary_structure_compatibility: float
    hydrophobic_burial_gain: float
    registry_support: float
    loop_entropy_cost: float
    geometry_violation_cost: float
    frustration_cost: float
    isolation_penalty: float
    nucleus_score: float
    closure_event_stability: float
    native_contact_count_after_scoring: int
    native_long_range_contact_count_after_scoring: int
    native_label_attached_after_event_generation: bool
    native_truth_used_before_event_generation: bool = False
    raw_sequence_exposed: bool = False

    def to_safe_dict(self) -> dict[str, object]:
        return asdict(self)

    def candidate_region_pairs(self) -> tuple[tuple[int, int], ...]:
        return tuple(
            (i, j)
            for i in range(self.segment_a_start, self.segment_a_end + 1)
            for j in range(self.segment_b_start, self.segment_b_end + 1)
            if j - i >= 3
        )


@dataclass(frozen=True)
class NucleusRowMetrics:
    row_id: str
    source_accession: str
    selected_threshold: float
    candidate_event_count: int
    accepted_event_count: int
    native_nucleus_recall: float
    false_nucleus_rate: float
    long_range_contact_recall_after_nucleus: float
    pair_level_long_range_contact_recall: float
    long_range_recall_delta_vs_pair_level: float
    contact_cluster_precision: float
    closure_event_stability: float
    trap_event_count: int
    frustration_rejection_count: int
    trajectory_native_gain: float
    nucleus_level_long_range_beats_pair_level: bool

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _mean(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)


def _hash_event(
    *,
    row_id: str,
    segment_a_start: int,
    segment_a_end: int,
    segment_b_start: int,
    segment_b_end: int,
) -> str:
    encoded = (
        f"{row_id}:{segment_a_start}-{segment_a_end}:"
        f"{segment_b_start}-{segment_b_end}"
    )
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def _top_features(features: Sequence[ContactLawFeatureRow]) -> tuple[ContactLawFeatureRow, ...]:
    return tuple(
        sorted(
            features,
            key=lambda feature: (
                -feature.pair_plus_cluster_score,
                -feature.cluster_neighbor_support,
                feature.i,
                feature.j,
            ),
        )[:10]
    )


def _event_score(
    *,
    contact_cluster_gain: float,
    secondary_structure_compatibility: float,
    hydrophobic_burial_gain: float,
    registry_support: float,
    loop_entropy_cost: float,
    geometry_violation_cost: float,
    frustration_cost: float,
    isolation_penalty: float,
) -> tuple[float, float]:
    score = (
        0.38 * contact_cluster_gain
        + 0.12 * secondary_structure_compatibility
        + 0.20 * hydrophobic_burial_gain
        + 0.14 * registry_support
        - 0.16 * loop_entropy_cost
        - 0.08 * geometry_violation_cost
        - 0.12 * frustration_cost
        - 0.06 * isolation_penalty
        + 0.12
    )
    nucleus_score = _rounded(score)
    stability = _rounded(
        nucleus_score
        * (1.0 - min(frustration_cost, 1.0))
        * (1.0 - min(geometry_violation_cost, 1.0))
    )
    return nucleus_score, stability


def nucleus_closure_events_for_row(
    row: RealCoordinateVisualRow,
    features: Sequence[ContactLawFeatureRow],
    *,
    segment_length: int = SEGMENT_LENGTH,
    segment_stride: int = SEGMENT_STRIDE,
) -> tuple[NucleusClosureEvent, ...]:
    feature_lookup = {(feature.i, feature.j): feature for feature in features}
    native_pairs = set(row.native_contact_pairs())
    events: list[NucleusClosureEvent] = []
    for segment_a_start in range(
        1,
        row.sequence_length - segment_length + 2,
        segment_stride,
    ):
        segment_a_end = segment_a_start + segment_length - 1
        first_b = segment_a_end + MIN_SEGMENT_GAP + 1
        for segment_b_start in range(
            first_b,
            row.sequence_length - segment_length + 2,
            segment_stride,
        ):
            segment_b_end = segment_b_start + segment_length - 1
            region_features = tuple(
                feature_lookup[(i, j)]
                for i in range(segment_a_start, segment_a_end + 1)
                for j in range(segment_b_start, segment_b_end + 1)
                if (i, j) in feature_lookup and j - i >= 3
            )
            if not region_features:
                continue
            top = _top_features(region_features)
            contact_cluster_gain = _rounded(
                _mean([feature.pair_plus_cluster_score for feature in top])
            )
            secondary_structure_compatibility = _rounded(
                _mean(
                    [
                        max(
                            feature.helix_window_support,
                            feature.beta_window_support,
                        )
                        for feature in top
                    ]
                )
            )
            hydrophobic_burial_gain = _rounded(
                _mean(
                    [
                        max(
                            feature.hydrophobic_pair_support,
                            feature.aromatic_anchor_support,
                        )
                        for feature in top
                    ]
                )
            )
            registry_support = _rounded(
                _mean([feature.parallel_contact_support for feature in top])
            )
            loop_entropy_cost = _rounded(
                _mean([feature.loop_entropy_cost for feature in top])
            )
            center_a = (segment_a_start + segment_a_end) / 2
            center_b = (segment_b_start + segment_b_end) / 2
            normalized_span = round(abs(center_b - center_a) / row.sequence_length, 6)
            geometry_violation_cost = _rounded(max(0.0, normalized_span - 0.65) * 1.3)
            frustration_cost = _rounded(
                _mean(
                    [
                        feature.same_charge_penalty + feature.breaker_penalty
                        for feature in top
                    ]
                )
            )
            isolation_penalty = _rounded(
                _mean([feature.isolation_penalty for feature in top])
            )
            nucleus_score, closure_event_stability = _event_score(
                contact_cluster_gain=contact_cluster_gain,
                secondary_structure_compatibility=secondary_structure_compatibility,
                hydrophobic_burial_gain=hydrophobic_burial_gain,
                registry_support=registry_support,
                loop_entropy_cost=loop_entropy_cost,
                geometry_violation_cost=geometry_violation_cost,
                frustration_cost=frustration_cost,
                isolation_penalty=isolation_penalty,
            )
            region_pairs = tuple(
                (i, j)
                for i in range(segment_a_start, segment_a_end + 1)
                for j in range(segment_b_start, segment_b_end + 1)
                if j - i >= 3
            )
            native_hits = native_pairs & set(region_pairs)
            events.append(
                NucleusClosureEvent(
                    row_id=row.row_id,
                    source_accession=row.source_accession,
                    sequence_hash=row.sequence_sha256,
                    sequence_length=row.sequence_length,
                    event_id=_hash_event(
                        row_id=row.row_id,
                        segment_a_start=segment_a_start,
                        segment_a_end=segment_a_end,
                        segment_b_start=segment_b_start,
                        segment_b_end=segment_b_end,
                    ),
                    segment_a_start=segment_a_start,
                    segment_a_end=segment_a_end,
                    segment_b_start=segment_b_start,
                    segment_b_end=segment_b_end,
                    sequence_span=segment_b_end - segment_a_start + 1,
                    normalized_span=normalized_span,
                    candidate_contact_count=len(region_pairs),
                    contact_cluster_gain=contact_cluster_gain,
                    secondary_structure_compatibility=(
                        secondary_structure_compatibility
                    ),
                    hydrophobic_burial_gain=hydrophobic_burial_gain,
                    registry_support=registry_support,
                    loop_entropy_cost=loop_entropy_cost,
                    geometry_violation_cost=geometry_violation_cost,
                    frustration_cost=frustration_cost,
                    isolation_penalty=isolation_penalty,
                    nucleus_score=nucleus_score,
                    closure_event_stability=closure_event_stability,
                    native_contact_count_after_scoring=len(native_hits),
                    native_long_range_contact_count_after_scoring=sum(
                        1 for left, right in native_hits if right - left >= 24
                    ),
                    native_label_attached_after_event_generation=True,
                )
            )
    return tuple(events)


def nucleus_closure_events(
    rows: Sequence[RealCoordinateVisualRow],
    feature_rows: Sequence[ContactLawFeatureRow],
) -> tuple[NucleusClosureEvent, ...]:
    features_by_row = feature_rows_by_row_id(feature_rows)
    output: list[NucleusClosureEvent] = []
    for row in rows:
        output.extend(
            nucleus_closure_events_for_row(row, features_by_row[row.row_id])
        )
    return tuple(output)


def accepted_events(
    events: Sequence[NucleusClosureEvent],
    *,
    threshold: float,
) -> tuple[NucleusClosureEvent, ...]:
    return tuple(
        event
        for event in events
        if event.nucleus_score >= threshold
        and event.frustration_cost < FRUSTRATION_REJECTION_LIMIT
    )


def _region_union(events: Sequence[NucleusClosureEvent]) -> set[tuple[int, int]]:
    output: set[tuple[int, int]] = set()
    for event in events:
        output.update(event.candidate_region_pairs())
    return output


def _pair_level_reference_metrics(
    features: Sequence[ContactLawFeatureRow],
) -> tuple[float, float]:
    native_pairs = tuple(
        (feature.i, feature.j)
        for feature in features
        if feature.native_contact
    )
    predicted_pairs = tuple(
        (feature.i, feature.j)
        for feature in features
        if feature.pair_plus_entropy_score >= PAIR_LEVEL_REFERENCE_THRESHOLD
    )
    metrics = evaluate_contact_prediction(
        native_pairs=native_pairs,
        predicted_pairs=predicted_pairs,
    )
    return metrics.long_range_contact_recall, metrics.contact_map_f1


def row_metrics_for_threshold(
    *,
    row: RealCoordinateVisualRow,
    row_features: Sequence[ContactLawFeatureRow],
    row_events: Sequence[NucleusClosureEvent],
    threshold: float,
) -> NucleusRowMetrics:
    native_pairs = set(row.native_contact_pairs())
    native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
    selected = accepted_events(row_events, threshold=threshold)
    selected_region = _region_union(selected)
    native_covered = selected_region & native_pairs
    native_long_covered = selected_region & native_long
    possible_region_pair_count = sum(event.candidate_contact_count for event in selected)
    native_hit_count = sum(event.native_contact_count_after_scoring for event in selected)
    false_events = [
        event
        for event in selected
        if event.native_contact_count_after_scoring == 0
    ]
    frustrated_rejections = [
        event
        for event in row_events
        if event.nucleus_score >= threshold
        and event.frustration_cost >= FRUSTRATION_REJECTION_LIMIT
    ]
    pair_long_range, _ = _pair_level_reference_metrics(row_features)
    trajectory = trajectory_rows_for_row(
        row=row,
        row_events=row_events,
        threshold=threshold,
    )
    first_gain = float(trajectory[0]["cumulative_native_nucleus_recall"]) if trajectory else 0.0
    final_gain = (
        float(trajectory[-1]["cumulative_native_nucleus_recall"])
        if trajectory
        else 0.0
    )
    long_range_recall = (
        len(native_long_covered) / len(native_long)
        if native_long
        else 1.0
    )
    return NucleusRowMetrics(
        row_id=row.row_id,
        source_accession=row.source_accession,
        selected_threshold=threshold,
        candidate_event_count=len(row_events),
        accepted_event_count=len(selected),
        native_nucleus_recall=_rounded(
            len(native_covered) / len(native_pairs) if native_pairs else 0.0
        ),
        false_nucleus_rate=_rounded(
            len(false_events) / len(selected) if selected else 0.0
        ),
        long_range_contact_recall_after_nucleus=_rounded(long_range_recall),
        pair_level_long_range_contact_recall=pair_long_range,
        long_range_recall_delta_vs_pair_level=_rounded(
            long_range_recall - pair_long_range
        ),
        contact_cluster_precision=_rounded(
            native_hit_count / possible_region_pair_count
            if possible_region_pair_count
            else 0.0
        ),
        closure_event_stability=_rounded(
            _mean([event.closure_event_stability for event in selected])
        ),
        trap_event_count=len(false_events),
        frustration_rejection_count=len(frustrated_rejections),
        trajectory_native_gain=_rounded(final_gain - first_gain),
        nucleus_level_long_range_beats_pair_level=(
            long_range_recall > pair_long_range
        ),
    )


def metrics_for_threshold(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    features_by_row: Mapping[str, Sequence[ContactLawFeatureRow]],
    events_by_row: Mapping[str, Sequence[NucleusClosureEvent]],
    threshold: float,
) -> list[NucleusRowMetrics]:
    return [
        row_metrics_for_threshold(
            row=row,
            row_features=features_by_row[row.row_id],
            row_events=events_by_row[row.row_id],
            threshold=threshold,
        )
        for row in rows
    ]


def aggregate_metrics(metrics: Sequence[NucleusRowMetrics]) -> dict[str, object]:
    return {
        "native_nucleus_recall": _rounded(
            _mean([metric.native_nucleus_recall for metric in metrics])
        ),
        "false_nucleus_rate": _rounded(
            _mean([metric.false_nucleus_rate for metric in metrics])
        ),
        "long_range_contact_recall_after_nucleus": _rounded(
            _mean(
                [
                    metric.long_range_contact_recall_after_nucleus
                    for metric in metrics
                ]
            )
        ),
        "pair_level_long_range_contact_recall": _rounded(
            _mean(
                [
                    metric.pair_level_long_range_contact_recall
                    for metric in metrics
                ]
            )
        ),
        "long_range_recall_delta_vs_pair_level": _rounded(
            _mean(
                [
                    metric.long_range_recall_delta_vs_pair_level
                    for metric in metrics
                ]
            )
        ),
        "contact_cluster_precision": _rounded(
            _mean([metric.contact_cluster_precision for metric in metrics])
        ),
        "closure_event_stability": _rounded(
            _mean([metric.closure_event_stability for metric in metrics])
        ),
        "trap_event_count": sum(metric.trap_event_count for metric in metrics),
        "frustration_rejection_count": sum(
            metric.frustration_rejection_count for metric in metrics
        ),
        "trajectory_native_gain": _rounded(
            _mean([metric.trajectory_native_gain for metric in metrics])
        ),
        "nucleus_level_long_range_beats_pair_level": all(
            metric.nucleus_level_long_range_beats_pair_level
            for metric in metrics
        ),
        "accepted_event_count": sum(metric.accepted_event_count for metric in metrics),
    }


def closure_quality_score(aggregate: Mapping[str, object]) -> float:
    return round(
        0.45 * float(aggregate["long_range_contact_recall_after_nucleus"])
        + 0.25 * float(aggregate["native_nucleus_recall"])
        + 0.15 * float(aggregate["contact_cluster_precision"])
        + 0.10 * float(aggregate["closure_event_stability"])
        - 0.25 * float(aggregate["false_nucleus_rate"]),
        6,
    )


def select_threshold(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    features_by_row: Mapping[str, Sequence[ContactLawFeatureRow]],
    events_by_row: Mapping[str, Sequence[NucleusClosureEvent]],
) -> tuple[float, list[NucleusRowMetrics], dict[str, object]]:
    best_threshold = 0.0
    best_metrics: list[NucleusRowMetrics] = []
    best_aggregate: dict[str, object] = {}
    best_score = -999.0
    for threshold in NUCLEUS_THRESHOLDS:
        metrics = metrics_for_threshold(
            rows=rows,
            features_by_row=features_by_row,
            events_by_row=events_by_row,
            threshold=threshold,
        )
        aggregate = aggregate_metrics(metrics)
        score = closure_quality_score(aggregate)
        if score > best_score:
            best_threshold = threshold
            best_metrics = metrics
            best_aggregate = {**aggregate, "closure_quality_score": score}
            best_score = score
    return best_threshold, best_metrics, best_aggregate


def trajectory_rows_for_row(
    *,
    row: RealCoordinateVisualRow,
    row_events: Sequence[NucleusClosureEvent],
    threshold: float,
) -> list[dict[str, object]]:
    native_pairs = set(row.native_contact_pairs())
    native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
    cumulative_region: set[tuple[int, int]] = set()
    rows = []
    selected = sorted(
        accepted_events(row_events, threshold=threshold),
        key=lambda event: (-event.nucleus_score, event.segment_a_start, event.segment_b_start),
    )
    for step, event in enumerate(selected, start=1):
        cumulative_region.update(event.candidate_region_pairs())
        native_covered = cumulative_region & native_pairs
        native_long_covered = cumulative_region & native_long
        rows.append(
            {
                "row_id": row.row_id,
                "source_accession": row.source_accession,
                "step": step,
                "event_id": event.event_id,
                "nucleus_score": event.nucleus_score,
                "closure_event_stability": event.closure_event_stability,
                "native_contact_count_after_scoring": (
                    event.native_contact_count_after_scoring
                ),
                "cumulative_native_nucleus_recall": _rounded(
                    len(native_covered) / len(native_pairs)
                    if native_pairs
                    else 0.0
                ),
                "cumulative_long_range_contact_recall": _rounded(
                    len(native_long_covered) / len(native_long)
                    if native_long
                    else 1.0
                ),
                "raw_sequence_exposed": False,
            }
        )
    return rows


def build_folding_nucleus_closure_report(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    feature_rows: Sequence[ContactLawFeatureRow],
    events: Sequence[NucleusClosureEvent],
    source_benchmark_file: Path,
) -> dict[str, object]:
    features_by_row = feature_rows_by_row_id(feature_rows)
    events_by_row: dict[str, tuple[NucleusClosureEvent, ...]] = {}
    for row in rows:
        events_by_row[row.row_id] = tuple(
            event for event in events if event.row_id == row.row_id
        )
    threshold, metrics, aggregate = select_threshold(
        rows=rows,
        features_by_row=features_by_row,
        events_by_row=events_by_row,
    )
    pair_long_range_values = [
        metric.pair_level_long_range_contact_recall for metric in metrics
    ]
    nucleus_law_survives = (
        bool(aggregate["nucleus_level_long_range_beats_pair_level"])
        and float(aggregate["false_nucleus_rate"]) < 0.35
        and float(aggregate["contact_cluster_precision"]) >= 0.10
    )
    failure_count = sum(
        1
        for metric in metrics
        if metric.false_nucleus_rate >= 0.50
        or metric.long_range_contact_recall_after_nucleus <= (
            metric.pair_level_long_range_contact_recall
        )
    )
    return {
        "report_kind": FOLDING_NUCLEUS_CLOSURE_REPORT_KIND,
        "source_benchmark_kind": REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "event_kind": FOLDING_NUCLEUS_CLOSURE_EVENT_KIND,
        "benchmark_size": len(rows),
        "segment_length": SEGMENT_LENGTH,
        "segment_stride": SEGMENT_STRIDE,
        "nucleus_threshold_min": NUCLEUS_THRESHOLD_MIN,
        "candidate_closure_event_count": len(events),
        "selected_threshold": threshold,
        **aggregate,
        "pair_level_reference_model": PAIR_LEVEL_REFERENCE_MODEL,
        "pair_level_reference_threshold": PAIR_LEVEL_REFERENCE_THRESHOLD,
        "pair_level_mean_long_range_contact_recall": _rounded(
            _mean(pair_long_range_values)
        ),
        "nucleus_level_scoring_completed": True,
        "nucleus_level_long_range_beats_pair_level": bool(
            aggregate["nucleus_level_long_range_beats_pair_level"]
        ),
        "nucleus_law_survives": nucleus_law_survives,
        "cooperative_closure_supported": bool(
            aggregate["nucleus_level_long_range_beats_pair_level"]
        ),
        "native_truth_used_before_event_generation": False,
        "native_label_attached_after_event_generation": True,
        "row_specific_nucleus_thresholds_forbidden": True,
        "candidate_law_failure_count": failure_count,
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
            "This benchmark tests segment-level cooperative closure events "
            "against coordinate-native contact regions. It shows whether "
            "closure nuclei recover long-range native regions better than "
            "pair-level thresholds, but it does not discover a folding law "
            "or solve protein folding."
        ),
        "metrics": [metric.to_dict() for metric in metrics],
    }


def metric_rows_from_report(report: Mapping[str, object]) -> list[dict[str, object]]:
    rows = report.get("metrics", [])
    if not isinstance(rows, list):
        return []
    return [dict(row) for row in rows if isinstance(row, Mapping)]


def event_rows(events: Sequence[NucleusClosureEvent]) -> list[dict[str, object]]:
    return [event.to_safe_dict() for event in events]


def trajectory_rows(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    events: Sequence[NucleusClosureEvent],
    threshold: float,
) -> list[dict[str, object]]:
    output = []
    for row in rows:
        row_events = tuple(event for event in events if event.row_id == row.row_id)
        output.extend(
            trajectory_rows_for_row(
                row=row,
                row_events=row_events,
                threshold=threshold,
            )
        )
    return output


def failure_rows_from_metrics(metrics: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    rows = []
    for metric in metrics:
        false_rate = float(metric["false_nucleus_rate"])
        long_delta = float(metric["long_range_recall_delta_vs_pair_level"])
        if false_rate < 0.50 and long_delta > 0:
            continue
        rows.append(
            {
                "row_id": metric["row_id"],
                "source_accession": metric["source_accession"],
                "false_nucleus_rate": metric["false_nucleus_rate"],
                "long_range_contact_recall_after_nucleus": metric[
                    "long_range_contact_recall_after_nucleus"
                ],
                "pair_level_long_range_contact_recall": metric[
                    "pair_level_long_range_contact_recall"
                ],
                "failure_reason": (
                    "false_nucleus_rate_high"
                    if false_rate >= 0.50
                    else "no_long_range_recall_gain"
                ),
                "mechanism_discovery_claim_allowed": False,
            }
        )
    return rows


def build_folding_nucleus_closure_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": FOLDING_NUCLEUS_CLOSURE_CERTIFICATE_KIND,
        "report_kind": report["report_kind"],
        "nucleus_level_scoring_completed": report[
            "nucleus_level_scoring_completed"
        ],
        "nucleus_level_long_range_beats_pair_level": report[
            "nucleus_level_long_range_beats_pair_level"
        ],
        "nucleus_law_survives": report["nucleus_law_survives"],
        "native_truth_used_before_event_generation": report[
            "native_truth_used_before_event_generation"
        ],
        "row_specific_nucleus_thresholds_forbidden": report[
            "row_specific_nucleus_thresholds_forbidden"
        ],
        "mechanism_discovery_claim_allowed": report[
            "mechanism_discovery_claim_allowed"
        ],
        "global_folding_claim_allowed": report["global_folding_claim_allowed"],
        "folding_problem_solved": report["folding_problem_solved"],
        "claim_allowed": report["claim_allowed"],
        "output_artifacts": tuple(output_names),
    }


def write_folding_nucleus_closure_outputs(
    *,
    report: Mapping[str, object],
    events: Sequence[NucleusClosureEvent],
    trajectories: Sequence[Mapping[str, object]],
    report_path: Path,
    events_path: Path,
    trajectory_path: Path,
    metrics_path: Path,
    failures_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    metric_rows = metric_rows_from_report(report)
    _write_csv_rows(event_rows(events), events_path)
    _write_csv_rows(trajectories, trajectory_path)
    _write_csv_rows(metric_rows, metrics_path)
    _write_csv_rows(failure_rows_from_metrics(metric_rows), failures_path)
    dashboard_path.write_text(
        render_folding_nucleus_closure_dashboard(report),
        encoding="utf-8",
    )
    certificate = build_folding_nucleus_closure_certificate(
        report,
        output_names=ROOT_OUTPUT_NAMES,
    )
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return (
        report_path,
        events_path,
        trajectory_path,
        metrics_path,
        failures_path,
        dashboard_path,
        certificate_path,
    )


def run_folding_nucleus_closure_search(
    *,
    benchmark_file: Path,
    report_path: Path,
    events_path: Path,
    trajectory_path: Path,
    metrics_path: Path,
    failures_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path, Path]:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    feature_rows = contact_law_feature_rows(rows)
    events = nucleus_closure_events(rows, feature_rows)
    report = build_folding_nucleus_closure_report(
        rows=rows,
        feature_rows=feature_rows,
        events=events,
        source_benchmark_file=benchmark_file,
    )
    trajectories = trajectory_rows(
        rows=rows,
        events=events,
        threshold=float(report["selected_threshold"]),
    )
    return write_folding_nucleus_closure_outputs(
        report=report,
        events=events,
        trajectories=trajectories,
        report_path=report_path,
        events_path=events_path,
        trajectory_path=trajectory_path,
        metrics_path=metrics_path,
        failures_path=failures_path,
        dashboard_path=dashboard_path,
        certificate_path=certificate_path,
    )


def _write_csv_rows(
    rows: Sequence[Mapping[str, object]],
    path: Path,
) -> Path:
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
        "candidate_closure_event_count",
        "selected_threshold",
        "native_nucleus_recall",
        "long_range_contact_recall_after_nucleus",
        "pair_level_mean_long_range_contact_recall",
        "false_nucleus_rate",
        "contact_cluster_precision",
        "nucleus_law_survives",
        "mechanism_discovery_claim_allowed",
        "folding_problem_solved",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
    )


def _metric_table(report: Mapping[str, object]) -> str:
    rows = metric_rows_from_report(report)
    body = []
    for row in rows:
        body.append(
            "<tr>"
            f"<td>{_escape(row['source_accession'])}</td>"
            f"<td>{_escape(row['accepted_event_count'])}</td>"
            f"<td>{_escape(row['native_nucleus_recall'])}</td>"
            f"<td>{_escape(row['long_range_contact_recall_after_nucleus'])}</td>"
            f"<td>{_escape(row['false_nucleus_rate'])}</td>"
            f"<td>{_escape(row['contact_cluster_precision'])}</td>"
            "</tr>"
        )
    return (
        "<section><h2>Row Metrics</h2>"
        "<table><thead><tr>"
        "<th>source</th><th>events</th><th>native recall</th>"
        "<th>long-range recall</th><th>false rate</th><th>cluster precision</th>"
        "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></section>"
    )


def _rule_cards() -> str:
    rules = (
        (
            "Closure Events, Not Pair Events",
            "Each candidate closes two sequence segments and is scored as a cooperative region.",
        ),
        (
            "Native Labels After Event Generation",
            "Coordinate-native contacts are used only after sequence-only event scoring.",
        ),
        (
            "Long-Range Recovery Is The First Target",
            "The benchmark asks whether segment closures recover long-range native regions better than pair thresholds.",
        ),
        (
            "No Folding Law Claim",
            "High false closure rates keep the nucleus law unclaimed.",
        ),
    )
    return "".join(
        "<div class=\"rule\">"
        f"<h3>{_escape(title)}</h3><p>{_escape(body)}</p>"
        "</div>"
        for title, body in rules
    )


def render_folding_nucleus_closure_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Folding Nucleus Closure Search</title>
  <style>
    body {{
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f5f6f1;
      color: #202623;
    }}
    header {{
      padding: 34px;
      background: #24302c;
      color: #f6f7f2;
    }}
    main {{
      max-width: 1220px;
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
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
    }}
    .metric, .rule {{
      background: #ffffff;
      border: 1px solid #d4ddd6;
      border-radius: 6px;
      padding: 14px;
    }}
    .metric span {{
      display: block;
      color: #58635e;
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
      border: 1px solid #d4ddd6;
    }}
    th, td {{
      padding: 9px 10px;
      border-bottom: 1px solid #e3e8e3;
      text-align: left;
      font-size: 13px;
      vertical-align: top;
      overflow-wrap: anywhere;
    }}
  </style>
</head>
<body>
  <header>
    <h1>Folding Nucleus Closure Search</h1>
    <p>Sequence-only segment closure events are scored against coordinate-native contact regions after event generation.</p>
  </header>
  <main>
    <section class="metrics">{_metric_cards(report)}</section>
    <section><h2>Boundary Rules</h2><div class="rules">{_rule_cards()}</div></section>
    {_metric_table(report)}
    <section>
      <h2>Claim Boundary</h2>
      <p>{_escape(report.get("boundary_statement", ""))}</p>
    </section>
  </main>
</body>
</html>
"""
