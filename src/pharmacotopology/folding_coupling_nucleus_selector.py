from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from math import log2, sqrt
from pathlib import Path
from statistics import mean
from typing import Mapping, Sequence

from pharmacotopology.artifact_io import write_csv_rows
from pharmacotopology.folding_contact_law_features import (
    contact_law_feature_rows,
    feature_rows_by_row_id,
)
from pharmacotopology.folding_coupling_decoy_falsification import (
    COUPLING_DECOY_FALSIFICATION_KIND,
    coupling_decoy_comparisons,
    real_beats_decoy_coupling_score_rate,
    real_vs_decoy_coupling_enrichment_ratio,
)
from pharmacotopology.folding_event_region_contact_collapse import (
    DEFAULT_BALANCED_PAIRS_PER_EVENT,
    EVENT_REGION_CONTACT_COLLAPSE_BOUNDARY,
    EVENT_REGION_CONTACT_COLLAPSE_KIND,
    SELF_DECIDING_STRATEGY_NAME,
    RowCollapseResult,
    collapse_event_region_contacts,
    collapse_row_event_regions,
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
from pharmacotopology.folding_native_contact_eval import contact_map_hash
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

# Inter-lobe contact detection thresholds (signal-based, protein-agnostic)
# If a protein contains any contact matching this strict signature,
# the entire protein is treated as multi-domain and gets adaptive gating.
COUPLING_INTERLOBE_NORMALIZED_SPAN_MIN = 0.25
COUPLING_INTERLOBE_DIRECT_SUPPORT_MAX = 0.15
COUPLING_INTERLOBE_CLUSTER_GAIN_MAX = 0.35

# Adaptive threshold safety rails for inter-lobe / low-signal contacts.
# The default remains conservative, but the effective floor can now move
# per row from phase shape + sequence complexity + MSA depth/coverage.
COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR_DEFAULT = 0.42
COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR_MIN = 0.35
COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR_MAX = 0.48
# Backwards-compatible name used by older tests/docs.
COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR = (
    COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR_DEFAULT
)
COUPLING_INTERLOBE_BLOCKED_FUTURE_CEILING = 0.25
COUPLING_ADAPTIVE_LOW_SIGNAL_FUTURE_MIN = 0.24
COUPLING_ADAPTIVE_LOW_SIGNAL_DIRECT_SUPPORT_MIN = 0.50
COUPLING_ADAPTIVE_LOW_SIGNAL_CLUSTER_MIN = 0.40
COUPLING_ADAPTIVE_LOW_SIGNAL_BLOCKED_FUTURE_CEILING = 0.35
COUPLING_ADAPTIVE_LOW_SIGNAL_PHYSICAL_SCORE_FLOOR_DEFAULT = 0.42
COUPLING_ADAPTIVE_LOW_SIGNAL_PHYSICAL_SCORE_FLOOR_MIN = 0.35
COUPLING_ADAPTIVE_LOW_SIGNAL_PHYSICAL_SCORE_FLOOR_MAX = 0.48
TRACE_LOOP_MARGIN_GATE_MIN = 0.0
TRACE_LOOP_MARGIN_GATE_BLOCKED_FUTURE_MAX = 0.16
TRACE_LOOP_TOP_RANK_FRACTION = 0.30
TRACE_LOOP_TOP_RANK_MIN_NEW_PAIRS = 2
TRACE_LOOP_CORE_TOP_RANK_FRACTION = 0.10
TRACE_LOOP_EXPANSION_TOP_RANK_FRACTION = 0.30
TRACE_LOOP_CLUSTER_GATE_MIN = 0.32
TRACE_LOOP_RANK_CONSISTENT_CLUSTER_GATE_MIN = 0.34
TRACE_LOOP_RANK_CONSISTENT_RECOVERY_CLUSTER_GATE_MIN = 0.32
TRACE_LOOP_RANK_CONSISTENT_RECOVERY_DIRECT_SUPPORT_MIN = 0.64
TRACE_LOOP_RANK_CONSISTENT_RECOVERY_FUTURE_PRESERVATION_MIN = 0.85
TRACE_LOOP_RANK_CONSISTENT_RECOVERY_DECOY_MARGIN_MIN = 0.30
TRACE_LOOP_RANK_CONSISTENT_RECOVERY_BLOCKED_FUTURE_MAX = 0.10
TRACE_LOOP_RANK_CONFIDENCE_CONSISTENCY_MIN = 0.95
TRACE_LOOP_SCORE_CONFIDENCE_CALIBRATION_MIN = 0.98
TRACE_LOOP_RANK_LENGTH_CALIBRATION_MIN = 0.999
TRACE_LOOP_PERSISTENT_RECOVERY_CLUSTER_GATE_MIN = 0.32
TRACE_LOOP_PERSISTENT_RECOVERY_DIRECT_SUPPORT_MIN = 0.64
TRACE_LOOP_PERSISTENT_RECOVERY_FUTURE_PRESERVATION_MIN = 0.78
TRACE_LOOP_PERSISTENT_RECOVERY_BLOCKED_FUTURE_MAX = 0.10
TRACE_LOOP_PERSISTENT_RECOVERY_SCORE_MIN = 0.70
TRACE_LOOP_PERSISTENT_RECOVERY_NEIGHBOR_COUNT_MIN = 2
TRACE_LOOP_PERSISTENT_NEIGHBOR_CLUSTER_MIN = 0.30
TRACE_LOOP_PERSISTENT_NEIGHBOR_DIRECT_SUPPORT_MIN = 0.45
TRACE_LOOP_PERSISTENT_NEIGHBOR_BLOCKED_FUTURE_MAX = 0.16
TRACE_LOOP_PERSISTENT_NEIGHBOR_WINDOW = 32
TRACE_LOOP_PERSISTENT_LOCAL_NEIGHBOR_WINDOW = 16
TRACE_LOOP_PERSISTENT_NEIGHBOR_LIMIT = 4
TRACE_LOOP_SCORE_MARGIN_EXPANSION_SCORE_MIN = 0.44
TRACE_LOOP_SCORE_MARGIN_EXPANSION_DECOY_MARGIN_MIN = 0.15
TRACE_LOOP_SCORE_MARGIN_EXPANSION_CLUSTER_MIN = 0.46
TRACE_LOOP_SCORE_MARGIN_EXPANSION_DIRECT_SUPPORT_MIN = 0.10
TRACE_LOOP_SCORE_MARGIN_EXPANSION_DIRECT_CONSTRAINT_COUNT_MIN = 2
TRACE_LOOP_SCORE_MARGIN_EXPANSION_FUTURE_PRESERVATION_MIN = 0.18
TRACE_LOOP_SCORE_MARGIN_EXPANSION_BLOCKED_FUTURE_MAX = 0.08
TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_SCORE_MIN = 0.45
TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_SCORE_MAX = 0.60
TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_DECOY_MARGIN_MIN = 0.18
TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_CLUSTER_MIN = 0.30
TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_DIRECT_SUPPORT_MIN = 0.25
TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_FUTURE_PRESERVATION_MIN = 0.50
TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_BLOCKED_FUTURE_MAX = 0.08
TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_SECONDARY_STRUCTURE_MIN = 0.53
TRACE_LOOP_EDGE_CONTINUITY_RESCUE_SCORE_MIN = 0.42
TRACE_LOOP_EDGE_CONTINUITY_RESCUE_SCORE_MAX = 0.50
TRACE_LOOP_EDGE_CONTINUITY_RESCUE_DECOY_MARGIN_MIN = 0.09
TRACE_LOOP_EDGE_CONTINUITY_RESCUE_CLUSTER_MIN = 0.43
TRACE_LOOP_EDGE_CONTINUITY_RESCUE_DIRECT_SUPPORT_MIN = 0.25
TRACE_LOOP_EDGE_CONTINUITY_RESCUE_FUTURE_PRESERVATION_MIN = 0.54
TRACE_LOOP_EDGE_CONTINUITY_RESCUE_BLOCKED_FUTURE_MAX = 0.08
TRACE_LOOP_EDGE_CONTINUITY_RESCUE_SECONDARY_STRUCTURE_MIN = 0.58
TRACE_LOOP_PRESSURE_RELEASE_RESCUE_SCORE_MIN = 0.42
TRACE_LOOP_PRESSURE_RELEASE_RESCUE_SCORE_MAX = 0.46
TRACE_LOOP_PRESSURE_RELEASE_RESCUE_DECOY_MARGIN_MIN = 0.14
TRACE_LOOP_PRESSURE_RELEASE_RESCUE_CLUSTER_MIN = 0.38
TRACE_LOOP_PRESSURE_RELEASE_RESCUE_DIRECT_SUPPORT_MIN = 0.40
TRACE_LOOP_PRESSURE_RELEASE_RESCUE_FUTURE_PRESERVATION_MIN = 0.30
TRACE_LOOP_PRESSURE_RELEASE_RESCUE_BLOCKED_FUTURE_MIN = 0.12
TRACE_LOOP_PRESSURE_RELEASE_RESCUE_BLOCKED_FUTURE_MAX = 0.22
TRACE_LOOP_PRESSURE_RELEASE_RESCUE_SECONDARY_STRUCTURE_MIN = 0.50
TRACE_LOOP_PRESSURE_RELEASE_RESCUE_DIRECT_CONSTRAINT_COUNT_MIN = 3
TRACE_LOOP_PRESSURE_RELEASE_RESCUE_DIRECT_CONFIDENCE_SUM_MIN = 1.20
TRACE_LOOP_PRESSURE_RELEASE_RESCUE_DIRECT_TOP_RANK_COUNT_MIN = 1
TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_SCORE_MIN = 0.40
TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_SCORE_MAX = 0.50
TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_DECOY_MARGIN_MIN = 0.0
TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_CLUSTER_MIN = 0.48
TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_DIRECT_SUPPORT_MIN = 0.10
TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_FUTURE_PRESERVATION_MIN = 0.18
TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_BLOCKED_FUTURE_MAX = 0.12
TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_SECONDARY_STRUCTURE_MIN = 0.50
TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_ANCHOR_WINDOW = 24
TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_ANCHOR_SCORE_MIN = 0.50
TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_ANCHOR_DIRECT_SUPPORT_MIN = 0.20
TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_SCORE_MIN = 0.40
TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_SCORE_MAX = 0.45
TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_DECOY_MARGIN_MIN = 0.04
TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_CLUSTER_MIN = 0.42
TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_DIRECT_SUPPORT_MIN = 0.25
TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_FUTURE_PRESERVATION_MIN = 0.18
TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_FUTURE_PRESERVATION_MAX = 0.25
TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_BLOCKED_FUTURE_MAX = 0.08
TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_COUPLING_MARGIN_MIN = 0.15
TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_SECONDARY_STRUCTURE_MIN = 0.53
TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_DIRECT_CONSTRAINT_COUNT_MIN = 4
TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_DIRECT_CONFIDENCE_SUM_MIN = 0.80
TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_SCORE_MIN = 0.30
TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_SCORE_MAX = 0.34
TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_DECOY_MARGIN_MIN = -0.01
TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_CLUSTER_MIN = 0.36
TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_DIRECT_SUPPORT_MIN = 0.20
TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_FUTURE_PRESERVATION_MIN = 0.22
TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_BLOCKED_FUTURE_MAX = 0.08
TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_COUPLING_MARGIN_MIN = 0.05
TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_SECONDARY_STRUCTURE_MIN = 0.65
TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_DIRECT_CONSTRAINT_COUNT_MIN = 3
TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_DIRECT_CONFIDENCE_SUM_MIN = 0.70
TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_SCORE_MIN = 0.35
TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_SCORE_MAX = 0.38
TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_DECOY_MARGIN_MIN = 0.05
TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_CLUSTER_MIN = 0.40
TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_DIRECT_SUPPORT_MIN = 0.15
TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_FUTURE_PRESERVATION_MIN = 0.30
TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_BLOCKED_FUTURE_MIN = 0.25
TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_BLOCKED_FUTURE_MAX = 0.35
TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_COUPLING_MARGIN_MIN = 0.04
TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_SECONDARY_STRUCTURE_MIN = 0.53
TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_DIRECT_CONSTRAINT_COUNT_MIN = 3
TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_DIRECT_CONFIDENCE_SUM_MIN = 0.50
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_DIRECT_CONSTRAINT_COUNT_MIN = 6
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_DIRECT_CONFIDENCE_SUM_MIN = 1.80
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_PRESERVATION_MIN = 0.70
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_BLOCKED_FUTURE_MAX = 0.10
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_SECONDARY_STRUCTURE_MIN = 0.50
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_CLUSTER_MIN = 0.36
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_SCORE_MIN = 0.50
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_COUPLING_MARGIN_MIN = -0.10
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_DIRECT_CONFIDENCE_DELTA_MIN = 0.20
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_DIRECT_SUPPORT_DROP_MAX = 0.05
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_PRESERVATION_DROP_MAX = 0.08
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_CLUSTER_DROP_MAX = 0.04
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_SECONDARY_STRUCTURE_DROP_MAX = 0.08
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_SPAN_CONTINUITY_DELTA_MAX = 8
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_NORMALIZED_SPAN_MAX = 0.55
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_SCORE_MIN = 0.49
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_PHYSICAL_MIN = 0.10
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_PHYSICAL_MAX = 0.50
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_CLUSTER_MIN = 0.34
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_SUPPORT_MIN = 0.40
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_PRESERVATION_MIN = 0.69
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_BLOCKED_MAX = 0.13
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_COUNT_MIN = 3
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_CONFIDENCE_MIN = 1.00
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_SCORE_MIN = 0.40
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_STATE_MIN = 0.65
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_BURIAL_MIN = 0.67
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_CLUSTER_MIN = 0.46
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_SECONDARY_STRUCTURE_MIN = 0.56
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_FUTURE_MIN = 0.17
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_BLOCKED_MAX = 0.11
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_DIRECT_COUNT_MIN = 1
TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_LIMIT_PER_ROW = 3
TRACE_LOOP_MACRO_SCALE_FUTURE_PRESERVATION_MIN = 0.32

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
    "coupling_nucleus_selector_contact_collapse_rows.csv",
    "coupling_nucleus_selector_collapsed_contacts.csv",
    "coupling_nucleus_selector_contact_collapse_events.csv",
    "coupling_nucleus_selector_dashboard.html",
    "coupling_nucleus_selector_certificate.json",
)
EXTERNAL_EVOLUTIONARY_COUPLING_SOURCE_KINDS = ACCEPTED_EXTERNAL_COUPLING_SOURCE_KINDS
PRIMARY_CONTACT_COLLAPSE_SELECTOR_NAME = "coupling_trace_loop"
PRIMARY_CONTACT_COLLAPSE_STRATEGY = SELF_DECIDING_STRATEGY_NAME

# Frontier expansion is now verified by the same native-free collapse layer that
# will later score the added region.  The controller only opens new frontier
# regions when the row already contains a broad, low-score ridge seed and the
# candidate survives an internal collapse-confidence distribution.
SELF_VERIFIED_EXPANSION_PROFILE = "direct_ridge_trace"
SELF_VERIFIED_EXPANSION_TIER_MODE = "direct_ridge_trace_low_score_tier"


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
    mean_selected_coupling_selectivity_score: float
    mean_decoy_coupling_selectivity_score: float
    mean_coupling_decoy_selectivity_margin: float
    mean_coupling_nucleus_score: float
    mean_decoy_coupling_nucleus_score: float
    mean_coupling_nucleus_decoy_margin: float
    real_vs_decoy_coupling_nucleus_enrichment_ratio: float
    real_beats_decoy_coupling_nucleus_score_rate: float
    survives_targets: bool
    coordinate_truth_used_to_build_constraints: bool
    native_truth_used_before_coupling_selection: bool
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class TraceLoopPersistenceEvidence:
    trace_loop_persistence_score: float
    persistent_neighbor_count: int
    mean_neighbor_direct_support: float
    mean_neighbor_future_preservation: float
    local_neighbor_fraction: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class AdaptiveCouplingFloorProfile:
    row_id: str
    source_accession: str
    phase_mode: str
    sequence_complexity: float
    coupling_depth_over_length: float
    target_coverage: float
    row_future_preservation_ceiling: float
    future_preservation_floor: float
    physical_score_floor: float
    adaptive_gate_enabled: bool
    low_signal_rescue_enabled: bool
    signal_reason: str

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


def _bounded_float(value: float, *, minimum: float, maximum: float) -> float:
    return max(minimum, min(maximum, value))


def _row_for_event(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> RealCoordinateVisualRow | None:
    for row in context.rows:
        if row.row_id == event.row_id:
            return row
    return None


def _constraints_for_row(
    row_id: str,
    context: CouplingNucleusContext,
) -> tuple[CouplingConstraint, ...]:
    return tuple(context.coupling_dataset.constraints_by_row_id().get(row_id, ()))


def _row_sequence_complexity(row: RealCoordinateVisualRow | None) -> float:
    if row is None or not row.sequence:
        return 0.0
    counts = Counter(residue for residue in row.sequence if residue)
    total = sum(counts.values())
    if total <= 0 or len(counts) <= 1:
        return 0.0
    entropy = -sum(
        (count / total) * log2(count / total)
        for count in counts.values()
        if count > 0
    )
    return _rounded(_bounded_float(entropy / log2(20), minimum=0.0, maximum=1.0))


def _population_stddev(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    center = mean(values)
    return sqrt(mean((value - center) ** 2 for value in values))


def _pair_correlation(pairs: Sequence[tuple[int, int]]) -> float:
    if len(pairs) <= 1:
        return 0.0
    left = [float(pair[0]) for pair in pairs]
    right = [float(pair[1]) for pair in pairs]
    left_center = mean(left)
    right_center = mean(right)
    left_ss = sum((value - left_center) ** 2 for value in left)
    right_ss = sum((value - right_center) ** 2 for value in right)
    if left_ss == 0.0 or right_ss == 0.0:
        return 0.0
    covariance = sum(
        (left_value - left_center) * (right_value - right_center)
        for left_value, right_value in zip(left, right)
    )
    return covariance / sqrt(left_ss * right_ss)


def _nearest_phase_votes(
    pairs: Sequence[tuple[int, int]],
) -> tuple[int, int]:
    diagonal_votes = 0
    square_votes = 0
    for pair in pairs:
        nearest: tuple[int, int, int] | None = None
        for other in pairs:
            if other == pair:
                continue
            delta_left = other[0] - pair[0]
            delta_right = other[1] - pair[1]
            distance = max(abs(delta_left), abs(delta_right))
            candidate = (distance, delta_left, delta_right)
            if nearest is None or candidate < nearest:
                nearest = candidate
        if nearest is None:
            continue
        distance, delta_left, delta_right = nearest
        if abs(delta_left - delta_right) <= 1:
            diagonal_votes += 1
        if distance <= 2:
            square_votes += 1
    return diagonal_votes, square_votes


def _phase_mode_from_pairs(pairs: Sequence[tuple[int, int]]) -> str:
    if len(pairs) <= 1:
        return "point"
    offsets = [float(right - left) for left, right in pairs]
    sums = [float(right + left) for left, right in pairs]
    offset_spread = _population_stddev(offsets)
    sum_spread = _population_stddev(sums)
    correlation = max(0.0, _pair_correlation(pairs))
    _, square_votes = _nearest_phase_votes(pairs)
    diagonal_score = correlation * sum_spread / (offset_spread + 1.0)
    square_score = (square_votes / len(pairs)) * sqrt(offset_spread + 1.0)
    if diagonal_score > square_score:
        return "diagonal"
    return "square"


def _row_phase_mode(
    row: RealCoordinateVisualRow | None,
    context: CouplingNucleusContext,
) -> str:
    if row is None:
        return "unknown"
    constraints = sorted(
        _constraints_for_row(row.row_id, context),
        key=lambda constraint: (
            -constraint.confidence,
            constraint.i,
            constraint.j,
            constraint.constraint_id,
        ),
    )
    top_count = min(24, len(constraints))
    if top_count <= 0:
        return "none"
    return _phase_mode_from_pairs(
        tuple(constraint.pair() for constraint in constraints[:top_count])
    )


def _row_coupling_quality(
    row_id: str,
    context: CouplingNucleusContext,
) -> tuple[float, float]:
    constraints = _constraints_for_row(row_id, context)
    depth_values = [
        float(constraint.effective_sequence_count_over_length)
        for constraint in constraints
        if constraint.effective_sequence_count_over_length > 0.0
    ]
    coverage_values = [
        float(constraint.target_coverage)
        for constraint in constraints
        if constraint.target_coverage > 0.0
    ]
    if depth_values:
        depth_over_length = mean(depth_values)
    else:
        row = next((candidate for candidate in context.rows if candidate.row_id == row_id), None)
        depth_over_length = (
            len(constraints) / max(1, row.sequence_length)
            if row is not None
            else 0.0
        )
    target_coverage = mean(coverage_values) if coverage_values else 0.0
    return _rounded(depth_over_length), _rounded(target_coverage)


def _row_future_preservation_ceiling(
    row_id: str,
    context: CouplingNucleusContext,
) -> float:
    values = [
        context.assessment_by_event_id[event.event_id].future_preservation_score
        for event in context.competitive_events
        if event.row_id == row_id
    ]
    return _rounded(max(values, default=0.0))


def _adaptive_coupling_floor_profile(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> AdaptiveCouplingFloorProfile:
    row = _row_for_event(event, context)
    source_accession = row.source_accession if row is not None else event.source_accession
    phase_mode = _row_phase_mode(row, context)
    sequence_complexity = _row_sequence_complexity(row)
    depth_over_length, target_coverage = _row_coupling_quality(event.row_id, context)
    future_ceiling = _row_future_preservation_ceiling(event.row_id, context)

    future_floor = COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR_DEFAULT
    reasons: list[str] = []
    if phase_mode == "diagonal":
        future_floor -= 0.04
        reasons.append("diagonal_phase")
    elif phase_mode == "square":
        future_floor -= 0.03
        reasons.append("square_phase")
    elif phase_mode == "point":
        future_floor += 0.03
        reasons.append("point_phase_guard")
    else:
        reasons.append(f"phase={phase_mode}")

    if sequence_complexity >= 0.90:
        future_floor -= 0.015
        reasons.append("high_sequence_complexity")
    elif sequence_complexity < 0.82:
        future_floor += 0.03
        reasons.append("low_sequence_complexity_guard")
    elif sequence_complexity < 0.86:
        future_floor += 0.015
        reasons.append("moderate_sequence_complexity_guard")

    if depth_over_length >= 12.0:
        future_floor -= 0.015
        reasons.append("deep_msa")
    elif depth_over_length >= 8.0:
        future_floor -= 0.01
        reasons.append("usable_msa_depth")
    elif 0.0 < depth_over_length < 3.0:
        future_floor += 0.04
        reasons.append("low_msa_depth_guard")

    if target_coverage >= 0.95:
        future_floor -= 0.005
        reasons.append("full_target_coverage")
    elif 0.0 < target_coverage < 0.85:
        future_floor += 0.03
        reasons.append("low_target_coverage_guard")

    if future_ceiling < COUPLING_ADAPTIVE_LOW_SIGNAL_FUTURE_MIN:
        future_floor = max(
            future_floor,
            COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR_DEFAULT,
        )
        reasons.append("low_future_ceiling_noise_guard")

    future_floor = _rounded(
        _bounded_float(
            future_floor,
            minimum=COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR_MIN,
            maximum=COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR_MAX,
        )
    )

    physical_floor = COUPLING_ADAPTIVE_LOW_SIGNAL_PHYSICAL_SCORE_FLOOR_DEFAULT
    if phase_mode in {"diagonal", "square"}:
        physical_floor -= 0.035
    if sequence_complexity >= 0.88:
        physical_floor -= 0.015
    elif sequence_complexity < 0.82:
        physical_floor += 0.03
    if depth_over_length >= 8.0:
        physical_floor -= 0.02
    elif 0.0 < depth_over_length < 3.0:
        physical_floor += 0.04
    if target_coverage >= 0.95:
        physical_floor -= 0.005
    elif 0.0 < target_coverage < 0.85:
        physical_floor += 0.03
    if future_ceiling < COUPLING_ADAPTIVE_LOW_SIGNAL_FUTURE_MIN:
        physical_floor += 0.05

    physical_floor = _rounded(
        _bounded_float(
            physical_floor,
            minimum=COUPLING_ADAPTIVE_LOW_SIGNAL_PHYSICAL_SCORE_FLOOR_MIN,
            maximum=COUPLING_ADAPTIVE_LOW_SIGNAL_PHYSICAL_SCORE_FLOOR_MAX,
        )
    )

    adaptive_gate_enabled = (
        future_floor < COUPLING_INTERLOBE_FUTURE_PRESERVATION_FLOOR_DEFAULT
        and depth_over_length >= 8.0
        and target_coverage >= 0.85
        and sequence_complexity >= 0.82
        and future_ceiling >= COUPLING_ADAPTIVE_LOW_SIGNAL_FUTURE_MIN
    )
    low_signal_rescue_enabled = (
        adaptive_gate_enabled
        and physical_floor <= 0.38
        and future_ceiling >= COUPLING_ADAPTIVE_LOW_SIGNAL_FUTURE_MIN
    )
    return AdaptiveCouplingFloorProfile(
        row_id=event.row_id,
        source_accession=source_accession,
        phase_mode=phase_mode,
        sequence_complexity=sequence_complexity,
        coupling_depth_over_length=depth_over_length,
        target_coverage=target_coverage,
        row_future_preservation_ceiling=future_ceiling,
        future_preservation_floor=future_floor,
        physical_score_floor=physical_floor,
        adaptive_gate_enabled=adaptive_gate_enabled,
        low_signal_rescue_enabled=low_signal_rescue_enabled,
        signal_reason=";".join(reasons),
    )


def adaptive_coupling_floor_report(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> dict[str, object]:
    return _adaptive_coupling_floor_profile(event, context).to_dict()


def _top_confidence_pairs(
    constraints: Sequence[CouplingConstraint],
    *,
    fraction: float,
) -> frozenset[tuple[int, int]]:
    ranked = sorted(
        constraints,
        key=lambda constraint: (
            -constraint.confidence,
            constraint.i,
            constraint.j,
            constraint.constraint_id,
        ),
    )
    top_count = max(1, int(len(ranked) * fraction))
    return frozenset(constraint.pair() for constraint in ranked[:top_count])


def _segment_trace_distance(
    left: NucleusClosureEvent,
    right: NucleusClosureEvent,
) -> int:
    return abs(left.segment_a_start - right.segment_a_start) + abs(
        left.segment_b_start - right.segment_b_start
    )


def row_rank_confidence_consistency(
    constraints: Sequence[CouplingConstraint],
) -> float:
    if not constraints:
        return 0.0
    if any(constraint.rank_fraction <= 0.0 for constraint in constraints):
        return 0.0
    ranked_by_confidence = sorted(
        constraints,
        key=lambda constraint: (
            -constraint.confidence,
            constraint.i,
            constraint.j,
            constraint.constraint_id,
        ),
    )
    denominator = max(1, len(ranked_by_confidence))
    confidence_rank_fraction = {
        constraint.constraint_id: rank / denominator
        for rank, constraint in enumerate(ranked_by_confidence, start=1)
    }
    mean_rank_drift = mean(
        abs(
            confidence_rank_fraction[constraint.constraint_id]
            - constraint.rank_fraction
        )
        for constraint in constraints
    )
    return _rounded(1.0 - 2.0 * mean_rank_drift)


def row_score_confidence_calibration(
    constraints: Sequence[CouplingConstraint],
) -> float:
    if not constraints:
        return 0.0
    max_positive_apc = max(
        (max(0.0, constraint.apc_corrected_score) for constraint in constraints),
        default=0.0,
    )
    if max_positive_apc <= 0.0:
        return 0.0
    mean_calibration_drift = mean(
        abs(
            constraint.confidence
            - max(0.01, max(0.0, constraint.apc_corrected_score) / max_positive_apc)
        )
        for constraint in constraints
    )
    return _rounded(1.0 - 2.0 * mean_calibration_drift)


def row_rank_length_calibration(
    *,
    sequence_length: int,
    constraints: Sequence[CouplingConstraint],
) -> float:
    if not constraints or sequence_length <= 0:
        return 0.0
    if any(constraint.rank <= 0 for constraint in constraints):
        return 0.0
    mean_rank_length_drift = mean(
        abs(constraint.rank_fraction - round(constraint.rank / sequence_length, 6))
        for constraint in constraints
    )
    return _rounded(1.0 - 100.0 * mean_rank_length_drift)


def build_coupling_nucleus_context(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    coupling_dataset: CouplingDataset,
    physical_context: ActivePhysicalContext | None = None,
) -> CouplingNucleusContext:
    physical_context = physical_context or build_active_physical_context(rows)
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


def _is_interlobe_contact_signature(
    event: NucleusClosureEvent,
    assessment: 'CouplingClosureAssessment',
) -> bool:
    """
    Detect inter-lobe contacts from their intrinsic signal characteristics,
    not from hardcoded PDB ID lists.

    Inter-lobe contacts in multi-domain proteins (e.g., calmodulin N-lobe/C-lobe)
    consistently show: (1) large normalized span, (2) weak direct evolutionary
    conservation, and (3) low local cluster density. This is because the two
    lobes evolve semi-independently — their coupling signal comes from allosteric
    coevolution rather than direct sequence conservation.
    """
    return (
        event.normalized_span >= COUPLING_INTERLOBE_NORMALIZED_SPAN_MIN
        and assessment.direct_support_score < COUPLING_INTERLOBE_DIRECT_SUPPORT_MAX
        and event.contact_cluster_gain < COUPLING_INTERLOBE_CLUSTER_GAIN_MAX
    )


def _adaptive_future_preservation_threshold(
    row_assessments: list['CouplingClosureAssessment'],
    *,
    default: float,
    floor: float,
) -> float:
    """
    Find the natural boundary in future_preservation_score distribution
    using the maximum-gap principle.

    Same principle as _multiscale_critical_boundary_items: sort scores
    descending, find the largest gap between consecutive values, and
    use the value just below the gap as the adaptive threshold.

    Safety: never goes below `floor` (prevents over-relaxation on noisy data).
    """
    scores = sorted(
        [a.future_preservation_score for a in row_assessments],
        reverse=True,
    )
    if len(scores) < 2:
        return default
    gaps = [scores[i] - scores[i + 1] for i in range(len(scores) - 1)]
    max_gap_index = max(range(len(gaps)), key=lambda i: gaps[i])
    # Threshold is set at the score just below the largest gap
    natural_threshold = scores[max_gap_index + 1]
    return max(floor, min(default, natural_threshold))





def _has_strict_interlobe_row_signature(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> bool:
    for candidate in context.competitive_events:
        if candidate.row_id != event.row_id:
            continue
        assessment = context.assessment_by_event_id[candidate.event_id]
        if _is_interlobe_contact_signature(candidate, assessment):
            return True
    return False


def _is_multidomain_protein(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> bool:
    """
    Dynamic adaptive-mode detector.

    The old gate switched only on a strict inter-lobe signature. That catches
    1CLL, but it misses weak/segmented rows such as 4AKE and underuses MSA
    quality when the signal is real but shallow. The new decision still accepts
    the strict signature, then adds a non-oracle row-level profile based on
    phase mode, sequence complexity, coupling depth and coverage.
    """
    if _has_strict_interlobe_row_signature(event, context):
        return True
    return _adaptive_coupling_floor_profile(event, context).adaptive_gate_enabled


def _passes_low_signal_adaptive_rescue(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
    profile: AdaptiveCouplingFloorProfile,
) -> bool:
    if not profile.low_signal_rescue_enabled:
        return False
    assessment = context.assessment_by_event_id[event.event_id]
    return (
        assessment.future_preservation_score >= COUPLING_ADAPTIVE_LOW_SIGNAL_FUTURE_MIN
        and coupling_nucleus_score(event, context) >= profile.physical_score_floor
        and assessment.direct_support_score
        >= COUPLING_ADAPTIVE_LOW_SIGNAL_DIRECT_SUPPORT_MIN
        and event.contact_cluster_gain >= COUPLING_ADAPTIVE_LOW_SIGNAL_CLUSTER_MIN
        and assessment.blocked_future_pressure
        <= COUPLING_ADAPTIVE_LOW_SIGNAL_BLOCKED_FUTURE_CEILING
    )


def _passes_coupling_future_gate(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> bool:
    assessment = context.assessment_by_event_id[event.event_id]
    profile = _adaptive_coupling_floor_profile(event, context)

    # Standard rows: keep the original strict gates. This is the 1TIM guard.
    if not _is_multidomain_protein(event, context):
        return (
            assessment.direct_support_score >= COUPLING_DIRECT_SUPPORT_MIN
            and assessment.future_preservation_score >= COUPLING_FUTURE_PRESERVATION_MIN
            and assessment.blocked_future_pressure <= COUPLING_BLOCKED_FUTURE_MAX
        )

    # Adaptive rows:
    # 1. Skip direct_support for inter-lobe-style contacts.
    # 2. Use a row gap threshold, but with a floor derived from phase/sequence/MSA.
    # 3. Allow a narrow low-signal rescue only when physical/coupling support is high.
    row_assessments = [
        context.assessment_by_event_id[e.event_id]
        for e in context.competitive_events
        if e.row_id == event.row_id
    ]
    adaptive_fp = _adaptive_future_preservation_threshold(
        row_assessments,
        default=COUPLING_FUTURE_PRESERVATION_MIN,
        floor=profile.future_preservation_floor,
    )
    adaptive_future_pass = (
        assessment.future_preservation_score >= adaptive_fp
        and assessment.blocked_future_pressure <= COUPLING_INTERLOBE_BLOCKED_FUTURE_CEILING
    )
    return adaptive_future_pass or _passes_low_signal_adaptive_rescue(
        event,
        context,
        profile,
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
    if selector_name == "coupling_trace_loop_top_rank_gated":
        return select_coupling_trace_loop_events(
            context,
            top_confidence_fraction=TRACE_LOOP_TOP_RANK_FRACTION,
            min_new_top_confidence_pairs=TRACE_LOOP_TOP_RANK_MIN_NEW_PAIRS,
        )
    if selector_name == "coupling_trace_loop_core_expanded":
        return select_coupling_trace_loop_core_expanded_events(context)
    if selector_name == "coupling_trace_loop_cluster_gated_core_expanded":
        return select_coupling_trace_loop_core_expanded_events(
            context,
            min_contact_cluster_gain=TRACE_LOOP_CLUSTER_GATE_MIN,
        )
    if selector_name == "coupling_trace_loop_rank_consistent_cluster_gated":
        return select_coupling_trace_loop_rank_consistent_cluster_gated_events(
            context
        )
    if (
        selector_name
        == "coupling_trace_loop_persistent_rank_consistent_cluster_gated"
    ):
        return select_coupling_trace_loop_persistent_rank_consistent_cluster_gated_events(
            context
        )
    if selector_name == "coupling_trace_loop_score_margin_expanded":
        return select_coupling_trace_loop_score_margin_expanded_events(context)
    if selector_name == "coupling_trace_loop_boundary_continuity_expanded":
        return select_coupling_trace_loop_boundary_continuity_expanded_events(
            context
        )
    if selector_name == "coupling_trace_loop_edge_continuity_expanded":
        return select_coupling_trace_loop_edge_continuity_expanded_events(context)
    if selector_name == "coupling_trace_loop_pressure_release_expanded":
        return select_coupling_trace_loop_pressure_release_expanded_events(context)
    if selector_name == "coupling_trace_loop_registry_extension_expanded":
        return select_coupling_trace_loop_registry_extension_expanded_events(context)
    if selector_name == "coupling_trace_loop_terminal_bridge_expanded":
        return select_coupling_trace_loop_terminal_bridge_expanded_events(context)
    if selector_name == "coupling_trace_loop_boundary_field_replacement_probe":
        return select_coupling_trace_loop_boundary_field_replacement_probe_events(
            context
        )
    if selector_name == "coupling_trace_loop_macro_scale_future_preserved":
        return select_coupling_trace_loop_macro_scale_future_preserved_events(
            context
        )
    if selector_name == "coupling_trace_loop_self_deciding_frontier_expanded":
        return select_coupling_trace_loop_self_deciding_frontier_expanded_events(context)

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
    top_confidence_fraction: float | None = None,
    min_new_top_confidence_pairs: int = 0,
    min_rank_confidence_consistency: float | None = None,
    min_score_confidence_calibration: float | None = None,
    min_rank_length_calibration: float | None = None,
) -> tuple[NucleusClosureEvent, ...]:
    constraints_by_row = context.coupling_dataset.constraints_by_row_id()
    competitive_by_row = _events_by_row(context.competitive_events)
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        row_constraints = tuple(constraints_by_row.get(row.row_id, ()))
        if (
            min_rank_confidence_consistency is not None
            and row_rank_confidence_consistency(row_constraints)
            < min_rank_confidence_consistency
        ):
            continue
        if (
            min_score_confidence_calibration is not None
            and row_score_confidence_calibration(row_constraints)
            < min_score_confidence_calibration
        ):
            continue
        if (
            min_rank_length_calibration is not None
            and row_rank_length_calibration(
                sequence_length=row.sequence_length,
                constraints=row_constraints,
            )
            < min_rank_length_calibration
        ):
            continue
        confidence_by_pair = _constraint_confidence_by_pair(row_constraints)
        top_confidence_pairs = (
            _top_confidence_pairs(row_constraints, fraction=top_confidence_fraction)
            if top_confidence_fraction is not None and row_constraints
            else frozenset()
        )
        uncovered = set(confidence_by_pair)
        row_selected: list[NucleusClosureEvent] = []
        selected_ids: set[str] = set()
        candidate_pairs_by_event_id = {
            event.event_id: set(event.candidate_region_pairs())
            for event in competitive_by_row.get(row.row_id, ())
        }
        row_candidates = tuple(competitive_by_row.get(row.row_id, ()))
        while uncovered and len(row_selected) < SELECTED_EVENTS_PER_ROW:
            scored_candidates: list[tuple[float, NucleusClosureEvent]] = []
            for event in row_candidates:
                if event.event_id in selected_ids:
                    continue
                if any(
                    not compatible_future_event(selected_event, event)
                    for selected_event in row_selected
                ):
                    continue
                event_pairs = candidate_pairs_by_event_id[event.event_id]
                newly_covered = event_pairs & uncovered
                if not newly_covered:
                    continue
                if (
                    min_new_top_confidence_pairs
                    and len(newly_covered & top_confidence_pairs)
                    < min_new_top_confidence_pairs
                ):
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
            selected_ids.add(chosen.event_id)
            uncovered -= candidate_pairs_by_event_id[chosen.event_id]
        selected.extend(row_selected)
    return tuple(selected)


def _compatible_merge_by_row(
    context: CouplingNucleusContext,
    *,
    core_events: Sequence[NucleusClosureEvent],
    expansion_events: Sequence[NucleusClosureEvent],
    min_contact_cluster_gain: float | None = None,
) -> tuple[NucleusClosureEvent, ...]:
    core_by_row = _events_by_row(core_events)
    expansion_by_row = _events_by_row(expansion_events)
    merged: list[NucleusClosureEvent] = []
    for row in context.rows:
        row_selected = list(core_by_row.get(row.row_id, ()))
        if min_contact_cluster_gain is not None:
            row_selected = [
                event
                for event in row_selected
                if event.contact_cluster_gain >= min_contact_cluster_gain
            ]
        selected_ids = {event.event_id for event in row_selected}
        for event in expansion_by_row.get(row.row_id, ()):
            if len(row_selected) >= SELECTED_EVENTS_PER_ROW:
                break
            if event.event_id in selected_ids:
                continue
            if (
                min_contact_cluster_gain is not None
                and event.contact_cluster_gain < min_contact_cluster_gain
            ):
                continue
            if any(
                not compatible_future_event(selected_event, event)
                for selected_event in row_selected
            ):
                continue
            row_selected.append(event)
            selected_ids.add(event.event_id)
        merged.extend(row_selected)
    return tuple(merged)


def select_coupling_trace_loop_core_expanded_events(
    context: CouplingNucleusContext,
    *,
    min_contact_cluster_gain: float | None = None,
    min_rank_confidence_consistency: float | None = None,
    min_score_confidence_calibration: float | None = None,
    min_rank_length_calibration: float | None = None,
) -> tuple[NucleusClosureEvent, ...]:
    core_events = select_coupling_trace_loop_events(
        context,
        top_confidence_fraction=TRACE_LOOP_CORE_TOP_RANK_FRACTION,
        min_new_top_confidence_pairs=TRACE_LOOP_TOP_RANK_MIN_NEW_PAIRS,
        min_rank_confidence_consistency=min_rank_confidence_consistency,
        min_score_confidence_calibration=min_score_confidence_calibration,
        min_rank_length_calibration=min_rank_length_calibration,
    )
    expansion_events = select_coupling_trace_loop_events(
        context,
        top_confidence_fraction=TRACE_LOOP_EXPANSION_TOP_RANK_FRACTION,
        min_new_top_confidence_pairs=TRACE_LOOP_TOP_RANK_MIN_NEW_PAIRS,
        min_rank_confidence_consistency=min_rank_confidence_consistency,
        min_score_confidence_calibration=min_score_confidence_calibration,
        min_rank_length_calibration=min_rank_length_calibration,
    )
    return _compatible_merge_by_row(
        context,
        core_events=core_events,
        expansion_events=expansion_events,
        min_contact_cluster_gain=min_contact_cluster_gain,
    )


def _passes_rank_consistent_recovery_gate(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> bool:
    assessment = context.assessment_by_event_id[event.event_id]
    coupling_decoy_margin = context.coupling_decoy_margin_by_event_id[event.event_id]
    return (
        event.contact_cluster_gain
        >= TRACE_LOOP_RANK_CONSISTENT_RECOVERY_CLUSTER_GATE_MIN
        and assessment.direct_support_score
        >= TRACE_LOOP_RANK_CONSISTENT_RECOVERY_DIRECT_SUPPORT_MIN
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_RANK_CONSISTENT_RECOVERY_BLOCKED_FUTURE_MAX
        and (
            assessment.future_preservation_score
            >= TRACE_LOOP_RANK_CONSISTENT_RECOVERY_FUTURE_PRESERVATION_MIN
            or coupling_decoy_margin
            >= TRACE_LOOP_RANK_CONSISTENT_RECOVERY_DECOY_MARGIN_MIN
        )
    )


def trace_loop_persistence_evidence(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
    candidate_events: Sequence[NucleusClosureEvent],
) -> TraceLoopPersistenceEvidence:
    neighbors: list[tuple[int, NucleusClosureEvent, CouplingClosureAssessment]] = []
    for candidate in candidate_events:
        if candidate.event_id == event.event_id or candidate.row_id != event.row_id:
            continue
        if not compatible_future_event(event, candidate):
            continue
        trace_distance = _segment_trace_distance(event, candidate)
        if trace_distance > TRACE_LOOP_PERSISTENT_NEIGHBOR_WINDOW:
            continue
        assessment = context.assessment_by_event_id[candidate.event_id]
        if candidate.contact_cluster_gain < TRACE_LOOP_PERSISTENT_NEIGHBOR_CLUSTER_MIN:
            continue
        if (
            assessment.direct_support_score
            < TRACE_LOOP_PERSISTENT_NEIGHBOR_DIRECT_SUPPORT_MIN
        ):
            continue
        if (
            assessment.blocked_future_pressure
            > TRACE_LOOP_PERSISTENT_NEIGHBOR_BLOCKED_FUTURE_MAX
        ):
            continue
        neighbors.append((trace_distance, candidate, assessment))

    if not neighbors:
        return TraceLoopPersistenceEvidence(
            trace_loop_persistence_score=0.0,
            persistent_neighbor_count=0,
            mean_neighbor_direct_support=0.0,
            mean_neighbor_future_preservation=0.0,
            local_neighbor_fraction=0.0,
        )

    strongest = tuple(
        sorted(
            neighbors,
            key=lambda item: (
                item[0],
                -item[2].direct_support_score,
                -item[2].future_preservation_score,
                item[1].event_id,
            ),
        )[:TRACE_LOOP_PERSISTENT_NEIGHBOR_LIMIT]
    )
    count_score = min(
        1.0,
        len(strongest) / TRACE_LOOP_PERSISTENT_RECOVERY_NEIGHBOR_COUNT_MIN,
    )
    mean_direct = mean(item[2].direct_support_score for item in strongest)
    mean_future = mean(item[2].future_preservation_score for item in strongest)
    local_fraction = sum(
        1
        for item in strongest
        if item[0] <= TRACE_LOOP_PERSISTENT_LOCAL_NEIGHBOR_WINDOW
    ) / len(strongest)
    persistence_score = _rounded(
        0.30 * count_score
        + 0.28 * mean_direct
        + 0.22 * mean_future
        + 0.12 * event.closure_event_stability
        + 0.08 * local_fraction
    )
    return TraceLoopPersistenceEvidence(
        trace_loop_persistence_score=persistence_score,
        persistent_neighbor_count=len(strongest),
        mean_neighbor_direct_support=_rounded(mean_direct),
        mean_neighbor_future_preservation=_rounded(mean_future),
        local_neighbor_fraction=_rounded(local_fraction),
    )


def _passes_persistent_recovery_gate(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
    candidate_events: Sequence[NucleusClosureEvent],
) -> bool:
    assessment = context.assessment_by_event_id[event.event_id]
    persistence = trace_loop_persistence_evidence(
        event,
        context,
        candidate_events,
    )
    return (
        event.contact_cluster_gain
        >= TRACE_LOOP_PERSISTENT_RECOVERY_CLUSTER_GATE_MIN
        and assessment.direct_support_score
        >= TRACE_LOOP_PERSISTENT_RECOVERY_DIRECT_SUPPORT_MIN
        and assessment.future_preservation_score
        >= TRACE_LOOP_PERSISTENT_RECOVERY_FUTURE_PRESERVATION_MIN
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_PERSISTENT_RECOVERY_BLOCKED_FUTURE_MAX
        and persistence.trace_loop_persistence_score
        >= TRACE_LOOP_PERSISTENT_RECOVERY_SCORE_MIN
        and persistence.persistent_neighbor_count
        >= TRACE_LOOP_PERSISTENT_RECOVERY_NEIGHBOR_COUNT_MIN
    )


def _selector_score_decoy_margin(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> float:
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
    return _rounded(
        coupling_nucleus_score(event, context)
        - coupling_nucleus_score(decoy, context)
    )


def _passes_score_margin_expansion_gate(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> bool:
    assessment = context.assessment_by_event_id[event.event_id]
    evidence = _direct_constraint_trace_evidence(event, context)
    return (
        coupling_nucleus_score(event, context)
        >= TRACE_LOOP_SCORE_MARGIN_EXPANSION_SCORE_MIN
        and _selector_score_decoy_margin(event, context)
        >= TRACE_LOOP_SCORE_MARGIN_EXPANSION_DECOY_MARGIN_MIN
        and event.contact_cluster_gain
        >= TRACE_LOOP_SCORE_MARGIN_EXPANSION_CLUSTER_MIN
        and assessment.direct_support_score
        >= TRACE_LOOP_SCORE_MARGIN_EXPANSION_DIRECT_SUPPORT_MIN
        and evidence["direct_constraint_count"]
        >= TRACE_LOOP_SCORE_MARGIN_EXPANSION_DIRECT_CONSTRAINT_COUNT_MIN
        and assessment.future_preservation_score
        >= TRACE_LOOP_SCORE_MARGIN_EXPANSION_FUTURE_PRESERVATION_MIN
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_SCORE_MARGIN_EXPANSION_BLOCKED_FUTURE_MAX
    )


def _passes_boundary_continuity_rescue_gate(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> bool:
    assessment = context.assessment_by_event_id[event.event_id]
    score = coupling_nucleus_score(event, context)
    return (
        score >= TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_SCORE_MIN
        and score <= TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_SCORE_MAX
        and _selector_score_decoy_margin(event, context)
        >= TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_DECOY_MARGIN_MIN
        and event.contact_cluster_gain
        >= TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_CLUSTER_MIN
        and assessment.direct_support_score
        >= TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_DIRECT_SUPPORT_MIN
        and assessment.future_preservation_score
        >= TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_FUTURE_PRESERVATION_MIN
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_BLOCKED_FUTURE_MAX
        and event.secondary_structure_compatibility
        >= TRACE_LOOP_BOUNDARY_CONTINUITY_RESCUE_SECONDARY_STRUCTURE_MIN
    )


def _passes_edge_continuity_rescue_gate(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> bool:
    assessment = context.assessment_by_event_id[event.event_id]
    score = coupling_nucleus_score(event, context)
    return (
        score >= TRACE_LOOP_EDGE_CONTINUITY_RESCUE_SCORE_MIN
        and score <= TRACE_LOOP_EDGE_CONTINUITY_RESCUE_SCORE_MAX
        and _selector_score_decoy_margin(event, context)
        >= TRACE_LOOP_EDGE_CONTINUITY_RESCUE_DECOY_MARGIN_MIN
        and event.contact_cluster_gain >= TRACE_LOOP_EDGE_CONTINUITY_RESCUE_CLUSTER_MIN
        and assessment.direct_support_score
        >= TRACE_LOOP_EDGE_CONTINUITY_RESCUE_DIRECT_SUPPORT_MIN
        and assessment.future_preservation_score
        >= TRACE_LOOP_EDGE_CONTINUITY_RESCUE_FUTURE_PRESERVATION_MIN
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_EDGE_CONTINUITY_RESCUE_BLOCKED_FUTURE_MAX
        and event.secondary_structure_compatibility
        >= TRACE_LOOP_EDGE_CONTINUITY_RESCUE_SECONDARY_STRUCTURE_MIN
    )


def _direct_constraint_trace_evidence(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> dict[str, float | int]:
    constraints = context.coupling_dataset.constraints_by_row_id().get(
        event.row_id,
        (),
    )
    region_pairs = set(event.candidate_region_pairs())
    direct = tuple(
        constraint for constraint in constraints if constraint.pair() in region_pairs
    )
    return {
        "direct_constraint_count": len(direct),
        "direct_constraint_confidence_sum": _rounded(
            sum(constraint.confidence for constraint in direct)
        ),
        "direct_top_10pct_rank_count": sum(
            1 for constraint in direct if constraint.rank_fraction <= 0.10
        ),
    }




def _self_deciding_frontier_expansion_score(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> dict[str, object]:
    """Native-free score for bringing extra event regions into the contact frontier.

    This is deliberately different from contact collapse.  It asks whether an
    event region is a good *candidate space* because it covers external DCA roots
    and has enough internal physical/future support to deserve collapse.  Native
    labels are not read here; they remain audit-only after selection.
    """
    constraints = context.coupling_dataset.constraints_by_row_id().get(
        event.row_id,
        (),
    )
    region_pairs = set(event.candidate_region_pairs())
    direct = tuple(
        constraint for constraint in constraints if constraint.pair() in region_pairs
    )
    assessment = context.assessment_by_event_id[event.event_id]
    direct_count = len(direct)
    direct_confidence_sum = sum(constraint.confidence for constraint in direct)
    direct_top_rank_count = sum(1 for constraint in direct if constraint.rank_fraction <= 0.10)
    confidence_density = direct_confidence_sum / max(1, len(region_pairs))
    coverage_score = min(1.0, direct_confidence_sum / 2.0)
    count_score = min(1.0, direct_count / 4.0)
    top_rank_score = min(1.0, direct_top_rank_count / 2.0)
    score = _rounded(
        0.36 * coverage_score
        + 0.22 * count_score
        + 0.16 * coupling_nucleus_score(event, context)
        + 0.12 * assessment.future_preservation_score
        + 0.10 * event.contact_cluster_gain
        + 0.08 * top_rank_score
        - 0.10 * assessment.blocked_future_pressure
    )
    return {
        "event_id": event.event_id,
        "row_id": event.row_id,
        "source_accession": event.source_accession,
        "self_deciding_frontier_expansion_score": score,
        "direct_constraint_count": direct_count,
        "direct_constraint_confidence_sum": _rounded(direct_confidence_sum),
        "direct_top_10pct_rank_count": direct_top_rank_count,
        "direct_constraint_confidence_density": _rounded(confidence_density),
        "coupling_nucleus_score": coupling_nucleus_score(event, context),
        "future_preservation_score": assessment.future_preservation_score,
        "blocked_future_pressure": assessment.blocked_future_pressure,
        "contact_cluster_gain": event.contact_cluster_gain,
        "native_truth_used_before_frontier_expansion": False,
        "coordinate_truth_used_before_frontier_expansion": False,
    }


def _collapse_inputs_for_context(
    context: CouplingNucleusContext,
) -> tuple[
    Mapping[str, Sequence[object]],
    Mapping[str, Sequence[CouplingConstraint]],
]:
    # The feature rows are sequence/contact-law derived.  They do not inspect the
    # native contact map and are safe for pre-evaluation controller decisions.
    features_by_row = feature_rows_by_row_id(contact_law_feature_rows(context.rows))
    constraints_by_row = context.coupling_dataset.constraints_by_row_id()
    return features_by_row, constraints_by_row


def _selected_pair_degree_consistency(selected_pairs) -> float:
    """Native-free residue-degree regularity for a collapsed candidate set.

    This is not the old fixed ``2..5 contacts`` rule.  It measures whether the
    collapsed set is spread across its own residue footprint instead of being
    dominated by a single hub.  The scale comes from the selected set itself.
    """
    if not selected_pairs:
        return 0.0
    degree: Counter[int] = Counter()
    for pair in selected_pairs:
        degree[int(pair.i)] += 1
        degree[int(pair.j)] += 1
    if not degree:
        return 0.0
    values = tuple(degree.values())
    max_degree = max(values)
    mean_degree = mean(values)
    footprint_fraction = len(values) / max(1.0, 2.0 * len(selected_pairs))
    hub_balance = mean_degree / max(1.0, float(max_degree))
    return _rounded(0.55 * footprint_fraction + 0.45 * hub_balance)


def _self_collapse_confidence_components(
    *,
    selected_pairs,
    summary,
) -> dict[str, float | str]:
    """Native-free confidence that a candidate region collapses cleanly.

    The confidence is phase-specific but still self-deciding: no native map,
    accession name, coordinate truth, fixed confidence threshold, or row-specific
    override enters this score.  The downstream selector still cuts by the
    largest internal gap in the row's candidate-confidence distribution.
    """
    selected_count = int(summary.selected_pair_count)
    candidate_count = int(summary.candidate_region_pair_count)
    if selected_count <= 0 or candidate_count <= 0:
        return {
            "self_collapse_confidence": 0.0,
            "self_collapse_phase_specific_confidence": 0.0,
            "self_collapse_degree_consistency": 0.0,
            "self_collapse_confidence_mode": "closed",
        }
    mean_pair_score = mean([float(pair.collapse_score) for pair in selected_pairs]) if selected_pairs else 0.0
    mean_ridge = mean([float(pair.ridge_coherence_score) for pair in selected_pairs]) if selected_pairs else 0.0
    mean_density = mean([float(pair.coupling_density_score) for pair in selected_pairs]) if selected_pairs else 0.0
    mean_sequence = mean([float(pair.sequence_law_support_score) for pair in selected_pairs]) if selected_pairs else 0.0
    mean_boundary = mean([float(pair.boundary_coherence_score) for pair in selected_pairs]) if selected_pairs else 0.0
    mean_internal_support = mean(
        [
            0.34 * float(pair.ridge_coherence_score)
            + 0.28 * float(pair.coupling_density_score)
            + 0.22 * float(pair.sequence_law_support_score)
            + 0.16 * float(pair.boundary_coherence_score)
            for pair in selected_pairs
        ]
    ) if selected_pairs else 0.0
    gap_signal = float(summary.self_deciding_gap_clarity) / max(
        1.0,
        float(summary.self_deciding_gap_clarity) + 12.0,
    )
    direct_root_signal = min(
        1.0,
        float(summary.direct_coupling_count_in_region) / max(1.0, sqrt(candidate_count)),
    )
    # Good expansion candidates are neither closed nor the full 8x8 block.  A
    # wide tier is allowed, but a near-total region is not trusted as expansion.
    width_signal = _rounded(1.0 - selected_count / max(1.0, candidate_count))
    degree_consistency = _selected_pair_degree_consistency(selected_pairs)
    profile_signal = 0.0
    if str(summary.self_deciding_profile) == SELF_VERIFIED_EXPANSION_PROFILE:
        profile_signal += 0.55
    if str(summary.self_deciding_tier_mode) == SELF_VERIFIED_EXPANSION_TIER_MODE:
        profile_signal += 0.45
    base_confidence = _rounded(
        0.22 * mean_pair_score
        + 0.20 * mean_internal_support
        + 0.18 * gap_signal
        + 0.14 * direct_root_signal
        + 0.12 * width_signal
        + 0.10 * degree_consistency
        + 0.04 * profile_signal
    )

    profile = str(summary.self_deciding_profile)
    phase = str(summary.self_deciding_phase_mode)
    if profile == "direct_ridge_trace":
        # Ridge traces trust coherent external-root density, but they also need a
        # sane degree footprint so low-score expansion does not reopen an 8x8 map.
        phase_confidence = _rounded(
            0.20 * mean_pair_score
            + 0.18 * mean_ridge
            + 0.18 * mean_density
            + 0.16 * mean_sequence
            + 0.14 * degree_consistency
            + 0.08 * direct_root_signal
            + 0.06 * width_signal
        )
        mode = "phase_direct_ridge_trace"
    elif phase == "sequence_inferred_lattice_or_beta":
        phase_confidence = _rounded(
            0.20 * mean_pair_score
            + 0.22 * mean_boundary
            + 0.18 * degree_consistency
            + 0.16 * mean_sequence
            + 0.12 * gap_signal
            + 0.12 * direct_root_signal
        )
        mode = "phase_lattice_boundary_degree"
    elif phase == "sequence_inferred_alpha_strip":
        phase_confidence = _rounded(
            0.20 * mean_pair_score
            + 0.24 * mean_sequence
            + 0.20 * mean_boundary
            + 0.16 * degree_consistency
            + 0.10 * direct_root_signal
            + 0.10 * width_signal
        )
        mode = "phase_alpha_strip_compactness"
    else:
        phase_confidence = _rounded(
            0.22 * mean_pair_score
            + 0.20 * mean_internal_support
            + 0.20 * degree_consistency
            + 0.16 * gap_signal
            + 0.12 * direct_root_signal
            + 0.10 * width_signal
        )
        mode = "phase_mixed_internal_evidence"

    return {
        "self_collapse_confidence": _rounded(0.45 * base_confidence + 0.55 * phase_confidence),
        "self_collapse_phase_specific_confidence": phase_confidence,
        "self_collapse_degree_consistency": degree_consistency,
        "self_collapse_confidence_mode": mode,
    }


def _self_collapse_confidence(
    *,
    selected_pairs,
    summary,
) -> float:
    return float(
        _self_collapse_confidence_components(
            selected_pairs=selected_pairs,
            summary=summary,
        )["self_collapse_confidence"]
    )


def _self_verified_collapse_row(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
    *,
    features_by_row: Mapping[str, Sequence[object]],
    constraints_by_row: Mapping[str, Sequence[CouplingConstraint]],
    cache: dict[str, dict[str, object]] | None = None,
) -> dict[str, object]:
    if cache is not None and event.event_id in cache:
        return cache[event.event_id]
    selected_pairs, summary = collapse_event_region_contacts(
        event=event,
        row_features=features_by_row.get(event.row_id, ()),  # type: ignore[arg-type]
        row_constraints=constraints_by_row.get(event.row_id, ()),
        collapse_strategy=SELF_DECIDING_STRATEGY_NAME,
        min_pairs_per_event=0,
        max_pairs_per_event=0,
    )
    confidence_components = _self_collapse_confidence_components(
        selected_pairs=selected_pairs,
        summary=summary,
    )
    confidence = float(confidence_components["self_collapse_confidence"])
    output = {
        **confidence_components,
        "self_collapse_selected_pair_count": summary.selected_pair_count,
        "self_collapse_candidate_region_pair_count": summary.candidate_region_pair_count,
        "self_collapse_profile": summary.self_deciding_profile,
        "self_collapse_phase_mode": summary.self_deciding_phase_mode,
        "self_collapse_tier_mode": summary.self_deciding_tier_mode,
        "self_collapse_decision_reason": summary.collapse_decision_reason,
        "self_collapse_gap_clarity": summary.self_deciding_gap_clarity,
        "self_collapse_cutoff_signal": summary.self_deciding_cutoff_signal,
        "self_collapse_direct_coupling_count": summary.direct_coupling_count_in_region,
        "self_collapse_native_truth_used_before_selection": False,
        "self_collapse_coordinate_truth_used_before_selection": False,
    }
    if cache is not None:
        cache[event.event_id] = output
    return output


def _self_verified_expandable_profiles(
    seed_events: Sequence[NucleusClosureEvent],
    context: CouplingNucleusContext,
    *,
    features_by_row: Mapping[str, Sequence[object]],
    constraints_by_row: Mapping[str, Sequence[CouplingConstraint]],
    cache: dict[str, dict[str, object]] | None = None,
) -> set[str]:
    profiles: set[str] = set()
    for event in seed_events:
        collapse_row = _self_verified_collapse_row(
            event,
            context,
            features_by_row=features_by_row,
            constraints_by_row=constraints_by_row,
            cache=cache,
        )
        if (
            int(collapse_row["self_collapse_selected_pair_count"]) > 0
            and str(collapse_row["self_collapse_profile"]) == SELF_VERIFIED_EXPANSION_PROFILE
            and str(collapse_row["self_collapse_tier_mode"]) == SELF_VERIFIED_EXPANSION_TIER_MODE
        ):
            profiles.add(str(collapse_row["self_collapse_profile"]))
    return profiles


def _self_verified_frontier_expansion_acceptance_score(
    score_row: Mapping[str, object],
    collapse_row: Mapping[str, object],
) -> float:
    """Native-free multi-evidence ranking score for expansion candidates.

    This replaces the earlier product-style score.  It is still not a threshold:
    candidates are ranked by this score and accepted only at the row-local largest
    gap.  The separate terms let phase-specific collapse confidence rescue a
    clean low-score ridge-shadow region without lowering a global floor.
    """
    expansion_score = float(score_row["self_deciding_frontier_expansion_score"])
    phase_confidence = float(collapse_row.get("self_collapse_phase_specific_confidence", 0.0))
    degree_consistency = float(collapse_row.get("self_collapse_degree_consistency", 0.0))
    selected_count = float(collapse_row.get("self_collapse_selected_pair_count", 0.0))
    candidate_count = float(collapse_row.get("self_collapse_candidate_region_pair_count", 0.0))
    breadth_signal = sqrt(selected_count / candidate_count) if candidate_count > 0.0 else 0.0
    tier_signal = (
        1.0
        if str(collapse_row.get("self_collapse_tier_mode", ""))
        == SELF_VERIFIED_EXPANSION_TIER_MODE
        else 0.0
    )
    return _rounded(
        0.38 * expansion_score
        + 0.22 * phase_confidence
        + 0.18 * breadth_signal
        + 0.10 * degree_consistency
        + 0.12 * tier_signal
    )




def _median_positive(values: Sequence[float]) -> float:
    positive = sorted(float(value) for value in values if float(value) > 0.0)
    if not positive:
        return 0.0
    midpoint = len(positive) // 2
    if len(positive) % 2:
        return _rounded(positive[midpoint])
    return _rounded((positive[midpoint - 1] + positive[midpoint]) / 2.0)


def _self_verified_identity_baseline_rows(
    seed_events: Sequence[NucleusClosureEvent],
    context: CouplingNucleusContext,
    *,
    features_by_row: Mapping[str, Sequence[object]],
    constraints_by_row: Mapping[str, Sequence[CouplingConstraint]],
    cache: dict[str, dict[str, object]] | None = None,
) -> dict[str, dict[str, object]]:
    """Return row-local, identity-derived baselines for expansion scoring.

    The baseline is not a global confidence threshold.  It is derived from the
    already accepted seed frontier for the same row and profile.  Candidate
    expansion regions are normalized against this seed identity before the
    largest-gap cutoff is applied.  Native/contact labels are not read here.
    """
    seed_by_row = _events_by_row(seed_events)
    output: dict[str, dict[str, object]] = {}
    for row_id, events in seed_by_row.items():
        by_profile: dict[str, list[float]] = defaultdict(list)
        all_scores: list[float] = []
        mode_counts: Counter[str] = Counter()
        for event in events:
            score_row = _self_deciding_frontier_expansion_score(event, context)
            collapse_row = _self_verified_collapse_row(
                event,
                context,
                features_by_row=features_by_row,
                constraints_by_row=constraints_by_row,
                cache=cache,
            )
            if int(collapse_row["self_collapse_selected_pair_count"]) <= 0:
                continue
            acceptance = _self_verified_frontier_expansion_acceptance_score(
                score_row,
                collapse_row,
            )
            profile = str(collapse_row["self_collapse_profile"])
            mode = str(collapse_row["self_collapse_confidence_mode"])
            by_profile[profile].append(acceptance)
            all_scores.append(acceptance)
            mode_counts[mode] += 1
        profile_baseline = {
            profile: _median_positive(scores)
            for profile, scores in by_profile.items()
            if _median_positive(scores) > 0.0
        }
        baseline = _median_positive(all_scores)
        lower_envelope = min((score for score in all_scores if score > 0.0), default=0.0)
        dominant_mode = mode_counts.most_common(1)[0][0] if mode_counts else "closed"
        output[row_id] = {
            "self_verified_identity_baseline": baseline,
            "self_verified_identity_lower_envelope": _rounded(lower_envelope),
            "self_verified_identity_profile_baselines": profile_baseline,
            "self_verified_identity_dominant_mode": dominant_mode,
            "self_verified_identity_seed_count": len(events),
            "native_truth_used_before_frontier_expansion": False,
            "coordinate_truth_used_before_frontier_expansion": False,
        }
    return output


def _identity_normalized_expansion_score(
    *,
    row_id: str,
    profile: str,
    acceptance_score: float,
    baseline_rows: Mapping[str, Mapping[str, object]],
) -> tuple[float, float, str]:
    baseline_row = baseline_rows.get(row_id, {})
    profile_baselines = baseline_row.get("self_verified_identity_profile_baselines", {})
    profile_baseline = 0.0
    if isinstance(profile_baselines, Mapping):
        profile_baseline = float(profile_baselines.get(profile, 0.0))
    row_baseline = float(baseline_row.get("self_verified_identity_baseline", 0.0))
    baseline = profile_baseline or row_baseline
    if baseline <= 0.0:
        return _rounded(acceptance_score), 0.0, "raw_score_no_seed_identity_baseline"
    return _rounded(acceptance_score / baseline), _rounded(baseline), "identity_normalized_seed_baseline"


def _self_deciding_frontier_expansion_cutoff(
    scores: Sequence[float],
    *,
    seed_lower_envelope: float | None = None,
) -> tuple[int, str, float]:
    """Return a native-free, gap-only expansion cutoff.

    Frontier expansion must not depend on an absolute confidence threshold such
    as 0.55 or 0.60.  The candidate regions are ranked by their own
    self-collapse acceptance score and the acceptance boundary is placed at the
    largest internal score gap.  Native/contact labels are intentionally absent
    from this decision.
    """
    if not scores:
        return 0, "empty_expansion_distribution", 0.0
    ordered = sorted(scores, reverse=True)
    if len(ordered) == 1:
        return 1, "single_self_verified_external_root_region", _rounded(ordered[0])
    gaps = [ordered[index] - ordered[index + 1] for index in range(len(ordered) - 1)]
    max_gap_index = max(range(len(gaps)), key=lambda index: gaps[index])
    max_gap = gaps[max_gap_index]
    if max_gap <= 0.0:
        return 0, "flat_self_collapse_distribution_no_internal_gap", _rounded(ordered[0])
    cutoff = min(SELECTED_EVENTS_PER_ROW, max_gap_index + 1)
    natural_boundary = ordered[cutoff - 1]
    if seed_lower_envelope is not None and seed_lower_envelope > 0.0:
        # Multi-tier identity expansion: after the strongest natural gap, keep
        # any additional candidate whose normalized score still lives inside the
        # seed frontier's own lower envelope.  This is a row-local identity
        # baseline, not an absolute confidence threshold.
        envelope_cutoff = sum(1 for score in ordered if score >= seed_lower_envelope)
        if envelope_cutoff > cutoff:
            cutoff = min(SELECTED_EVENTS_PER_ROW, envelope_cutoff)
            natural_boundary = ordered[cutoff - 1]
            return cutoff, "identity_normalized_multitier_gap_frontier", _rounded(natural_boundary)
    return cutoff, "self_collapse_verified_internal_gap_frontier", _rounded(natural_boundary)


def select_coupling_trace_loop_self_deciding_frontier_expanded_events(
    context: CouplingNucleusContext,
    *,
    seed_events: Sequence[NucleusClosureEvent] = (),
) -> tuple[NucleusClosureEvent, ...]:
    """Select a self-verified external-root frontier expansion.

    This is not a raw low-floor expansion.  It starts from the already accepted
    frontier and opens extra regions only for rows whose seed frontier contains a
    broad, low-score ridge trace.  Each candidate is collapsed first with the
    native-free self-deciding collapse layer, then ranked by its own expansion
    evidence and collapse confidence distribution.  Native labels are not read.
    """
    features_by_row, constraints_by_row = _collapse_inputs_for_context(context)
    collapse_cache: dict[str, dict[str, object]] = {}
    seed_by_row = _events_by_row(seed_events)
    identity_baseline_rows = _self_verified_identity_baseline_rows(
        seed_events,
        context,
        features_by_row=features_by_row,
        constraints_by_row=constraints_by_row,
        cache=collapse_cache,
    )
    selected: list[NucleusClosureEvent] = []
    competitive_by_row = _events_by_row(context.competitive_events)
    for row in context.rows:
        row_selected: list[NucleusClosureEvent] = list(seed_by_row.get(row.row_id, ()))
        selected_ids = {event.event_id for event in row_selected}
        expandable_profiles = _self_verified_expandable_profiles(
            row_selected,
            context,
            features_by_row=features_by_row,
            constraints_by_row=constraints_by_row,
            cache=collapse_cache,
        )
        if not expandable_profiles:
            selected.extend(row_selected)
            continue
        candidates: list[tuple[float, NucleusClosureEvent, dict[str, object]]] = []
        for event in competitive_by_row.get(row.row_id, ()):
            if event.event_id in selected_ids:
                continue
            score_row = _self_deciding_frontier_expansion_score(event, context)
            if int(score_row["direct_constraint_count"]) <= 0:
                continue
            collapse_row = _self_verified_collapse_row(
                event,
                context,
                features_by_row=features_by_row,
                constraints_by_row=constraints_by_row,
                cache=collapse_cache,
            )
            if int(collapse_row["self_collapse_selected_pair_count"]) <= 0:
                continue
            if str(collapse_row["self_collapse_profile"]) not in expandable_profiles:
                continue
            acceptance_score = _self_verified_frontier_expansion_acceptance_score(
                score_row,
                collapse_row,
            )
            normalized_score, identity_baseline, normalization_reason = (
                _identity_normalized_expansion_score(
                    row_id=event.row_id,
                    profile=str(collapse_row["self_collapse_profile"]),
                    acceptance_score=acceptance_score,
                    baseline_rows=identity_baseline_rows,
                )
            )
            merged_row = {
                **score_row,
                **collapse_row,
                "self_verified_frontier_expansion_acceptance_score": acceptance_score,
                "self_verified_frontier_expansion_identity_normalized_score": normalized_score,
                "self_verified_frontier_expansion_identity_baseline": identity_baseline,
                "self_verified_frontier_expansion_normalization_reason": normalization_reason,
                "self_verified_frontier_expansion_gate_reason": "matching_seed_low_score_ridge_trace_self_collapse_survived",
                "native_truth_used_before_frontier_expansion": False,
                "coordinate_truth_used_before_frontier_expansion": False,
            }
            candidates.append((normalized_score, event, merged_row))
        candidates.sort(
            key=lambda item: (
                -item[0],
                -float(item[2]["self_verified_frontier_expansion_acceptance_score"]),
                -float(item[2]["self_deciding_frontier_expansion_score"]),
                -float(item[2]["self_collapse_confidence"]),
                item[1].segment_a_start,
                item[1].segment_b_start,
                item[1].event_id,
            )
        )
        baseline_row = identity_baseline_rows.get(row.row_id, {})
        identity_baseline = float(baseline_row.get("self_verified_identity_baseline", 0.0))
        lower_envelope = float(baseline_row.get("self_verified_identity_lower_envelope", 0.0))
        normalized_lower_envelope = (
            _rounded(lower_envelope / identity_baseline)
            if identity_baseline > 0.0 and lower_envelope > 0.0
            else None
        )
        cutoff, _reason, _natural_boundary = _self_deciding_frontier_expansion_cutoff(
            [score for score, _event, _score_row in candidates],
            seed_lower_envelope=normalized_lower_envelope,
        )
        for _score, event, _score_row in candidates[:cutoff]:
            if len(row_selected) >= SELECTED_EVENTS_PER_ROW:
                break
            if event.event_id in selected_ids:
                continue
            row_selected.append(event)
            selected_ids.add(event.event_id)
        selected.extend(row_selected)
    return tuple(selected)


def _closed_self_verified_collapse_row() -> dict[str, object]:
    return {
        "self_collapse_confidence": 0.0,
        "self_collapse_phase_specific_confidence": 0.0,
        "self_collapse_degree_consistency": 0.0,
        "self_collapse_confidence_mode": "closed",
        "self_collapse_selected_pair_count": 0,
        "self_collapse_candidate_region_pair_count": 0,
        "self_collapse_profile": "not_evaluated_without_seed_expandable_profile",
        "self_collapse_phase_mode": "",
        "self_collapse_tier_mode": "closed",
        "self_collapse_decision_reason": "row_has_no_self_verified_expandable_seed",
        "self_collapse_gap_clarity": 0.0,
        "self_collapse_cutoff_signal": 0.0,
        "self_collapse_direct_coupling_count": 0,
        "self_collapse_native_truth_used_before_selection": False,
        "self_collapse_coordinate_truth_used_before_selection": False,
    }


def self_deciding_frontier_expansion_rows(
    context: CouplingNucleusContext,
    *,
    seed_events: Sequence[NucleusClosureEvent] = (),
) -> list[dict[str, object]]:
    features_by_row, constraints_by_row = _collapse_inputs_for_context(context)
    collapse_cache: dict[str, dict[str, object]] = {}
    seed_by_row = _events_by_row(seed_events)
    identity_baseline_rows = _self_verified_identity_baseline_rows(
        seed_events,
        context,
        features_by_row=features_by_row,
        constraints_by_row=constraints_by_row,
        cache=collapse_cache,
    )
    seed_ids = {event.event_id for event in seed_events}
    selected_ids = {
        event.event_id
        for event in select_coupling_trace_loop_self_deciding_frontier_expanded_events(
            context,
            seed_events=seed_events,
        )
    }
    expandable_profiles_by_row = {
        row_id: _self_verified_expandable_profiles(
            events,
            context,
            features_by_row=features_by_row,
            constraints_by_row=constraints_by_row,
            cache=collapse_cache,
        )
        for row_id, events in seed_by_row.items()
    }
    rows: list[dict[str, object]] = []
    for event in context.competitive_events:
        score_row = _self_deciding_frontier_expansion_score(event, context)
        if int(score_row["direct_constraint_count"]) <= 0 and event.event_id not in selected_ids:
            continue
        expandable_profiles = expandable_profiles_by_row.get(event.row_id, set())
        if event.event_id in seed_ids or expandable_profiles:
            collapse_row = _self_verified_collapse_row(
                event,
                context,
                features_by_row=features_by_row,
                constraints_by_row=constraints_by_row,
                cache=collapse_cache,
            )
        else:
            collapse_row = _closed_self_verified_collapse_row()
        candidate_profile_allowed = (
            str(collapse_row["self_collapse_profile"]) in expandable_profiles
            and int(collapse_row["self_collapse_selected_pair_count"]) > 0
        )
        acceptance_score = (
            _self_verified_frontier_expansion_acceptance_score(score_row, collapse_row)
            if candidate_profile_allowed
            else 0.0
        )
        normalized_score, identity_baseline, normalization_reason = (
            _identity_normalized_expansion_score(
                row_id=event.row_id,
                profile=str(collapse_row["self_collapse_profile"]),
                acceptance_score=acceptance_score,
                baseline_rows=identity_baseline_rows,
            )
            if candidate_profile_allowed
            else (0.0, 0.0, "not_candidate_for_identity_normalization")
        )
        baseline_row = identity_baseline_rows.get(event.row_id, {})
        score_row.update(
            collapse_row
        )
        score_row.update(
            {
                "seed_event": event.event_id in seed_ids,
                "self_verified_candidate_profile_allowed": candidate_profile_allowed,
                "self_verified_frontier_expansion_acceptance_score": acceptance_score,
                "self_verified_frontier_expansion_identity_normalized_score": normalized_score,
                "self_verified_frontier_expansion_identity_baseline": identity_baseline,
                "self_verified_frontier_expansion_identity_lower_envelope": baseline_row.get(
                    "self_verified_identity_lower_envelope",
                    0.0,
                ),
                "self_verified_frontier_expansion_identity_dominant_mode": baseline_row.get(
                    "self_verified_identity_dominant_mode",
                    "closed",
                ),
                "self_verified_frontier_expansion_normalization_reason": normalization_reason,
                "self_deciding_frontier_expansion_selected": event.event_id in selected_ids,
                "self_verified_frontier_expansion_selected": event.event_id in selected_ids,
                "self_verified_frontier_expansion_gate_reason": (
                    "seed_event"
                    if event.event_id in seed_ids
                    else "matching_seed_low_score_ridge_trace_self_collapse_survived"
                    if candidate_profile_allowed
                    else "closed_or_profile_not_seed_verified"
                ),
                "native_truth_used_before_frontier_expansion": False,
                "coordinate_truth_used_before_frontier_expansion": False,
            }
        )
        rows.append(score_row)

    candidate_rows_by_row: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        if row["seed_event"] is True:
            continue
        if row["self_verified_candidate_profile_allowed"] is not True:
            continue
        candidate_rows_by_row[str(row["row_id"])].append(row)

    cutoff_metadata_by_row: dict[str, tuple[int, str, float]] = {}
    for row_id, candidate_rows in candidate_rows_by_row.items():
        ordered_rows = sorted(
            candidate_rows,
            key=lambda row: (
                -float(row["self_verified_frontier_expansion_identity_normalized_score"]),
                -float(row["self_verified_frontier_expansion_acceptance_score"]),
                -float(row["self_deciding_frontier_expansion_score"]),
                -float(row["self_collapse_confidence"]),
                str(row["event_id"]),
            ),
        )
        baseline_row = identity_baseline_rows.get(row_id, {})
        identity_baseline = float(baseline_row.get("self_verified_identity_baseline", 0.0))
        lower_envelope = float(baseline_row.get("self_verified_identity_lower_envelope", 0.0))
        normalized_lower_envelope = (
            _rounded(lower_envelope / identity_baseline)
            if identity_baseline > 0.0 and lower_envelope > 0.0
            else None
        )
        cutoff, reason, natural_boundary = _self_deciding_frontier_expansion_cutoff(
            [
                float(row["self_verified_frontier_expansion_identity_normalized_score"])
                for row in ordered_rows
            ],
            seed_lower_envelope=normalized_lower_envelope,
        )
        cutoff_metadata_by_row[row_id] = (cutoff, reason, natural_boundary)
        for rank, row in enumerate(ordered_rows, start=1):
            row["self_verified_frontier_expansion_gap_rank"] = rank
            row["self_verified_frontier_expansion_cutoff_count"] = cutoff
            row["self_verified_frontier_expansion_cutoff_reason"] = reason
            row["self_verified_frontier_expansion_natural_boundary"] = natural_boundary
            row["self_verified_frontier_expansion_gap_selected"] = rank <= cutoff

    for row in rows:
        if "self_verified_frontier_expansion_cutoff_reason" in row:
            continue
        cutoff, reason, natural_boundary = cutoff_metadata_by_row.get(
            str(row["row_id"]),
            (0, "not_candidate_for_gap_cutoff", 0.0),
        )
        row["self_verified_frontier_expansion_gap_rank"] = 0
        row["self_verified_frontier_expansion_cutoff_count"] = cutoff
        row["self_verified_frontier_expansion_cutoff_reason"] = reason
        row["self_verified_frontier_expansion_natural_boundary"] = natural_boundary
        row["self_verified_frontier_expansion_gap_selected"] = False

    return sorted(
        rows,
        key=lambda row: (
            str(row["row_id"]),
            -float(row["self_verified_frontier_expansion_acceptance_score"]),
            -float(row["self_deciding_frontier_expansion_score"]),
            str(row["event_id"]),
        ),
    )


def _passes_pressure_release_rescue_gate(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> bool:
    assessment = context.assessment_by_event_id[event.event_id]
    score = coupling_nucleus_score(event, context)
    evidence = _direct_constraint_trace_evidence(event, context)
    return (
        score >= TRACE_LOOP_PRESSURE_RELEASE_RESCUE_SCORE_MIN
        and score <= TRACE_LOOP_PRESSURE_RELEASE_RESCUE_SCORE_MAX
        and _selector_score_decoy_margin(event, context)
        >= TRACE_LOOP_PRESSURE_RELEASE_RESCUE_DECOY_MARGIN_MIN
        and event.contact_cluster_gain
        >= TRACE_LOOP_PRESSURE_RELEASE_RESCUE_CLUSTER_MIN
        and assessment.direct_support_score
        >= TRACE_LOOP_PRESSURE_RELEASE_RESCUE_DIRECT_SUPPORT_MIN
        and assessment.future_preservation_score
        >= TRACE_LOOP_PRESSURE_RELEASE_RESCUE_FUTURE_PRESERVATION_MIN
        and assessment.blocked_future_pressure
        >= TRACE_LOOP_PRESSURE_RELEASE_RESCUE_BLOCKED_FUTURE_MIN
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_PRESSURE_RELEASE_RESCUE_BLOCKED_FUTURE_MAX
        and event.secondary_structure_compatibility
        >= TRACE_LOOP_PRESSURE_RELEASE_RESCUE_SECONDARY_STRUCTURE_MIN
        and evidence["direct_constraint_count"]
        >= TRACE_LOOP_PRESSURE_RELEASE_RESCUE_DIRECT_CONSTRAINT_COUNT_MIN
        and evidence["direct_constraint_confidence_sum"]
        >= TRACE_LOOP_PRESSURE_RELEASE_RESCUE_DIRECT_CONFIDENCE_SUM_MIN
        and evidence["direct_top_10pct_rank_count"]
        >= TRACE_LOOP_PRESSURE_RELEASE_RESCUE_DIRECT_TOP_RANK_COUNT_MIN
    )


def _registry_extension_anchor(
    event: NucleusClosureEvent,
    selected_events: Sequence[NucleusClosureEvent],
    context: CouplingNucleusContext,
) -> NucleusClosureEvent | None:
    anchors = tuple(
        selected
        for selected in selected_events
        if (
            selected.segment_b_start == event.segment_b_start
            and abs(selected.segment_a_start - event.segment_a_start)
            <= TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_ANCHOR_WINDOW
        )
        or (
            selected.segment_a_start == event.segment_a_start
            and abs(selected.segment_b_start - event.segment_b_start)
            <= TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_ANCHOR_WINDOW
        )
    )
    if not anchors:
        return None
    return max(
        anchors,
        key=lambda selected: (
            coupling_nucleus_score(selected, context),
            -_segment_trace_distance(selected, event),
            selected.event_id,
        ),
    )


def _passes_registry_extension_rescue_gate(
    event: NucleusClosureEvent,
    selected_events: Sequence[NucleusClosureEvent],
    context: CouplingNucleusContext,
) -> bool:
    assessment = context.assessment_by_event_id[event.event_id]
    score = coupling_nucleus_score(event, context)
    anchor = _registry_extension_anchor(event, selected_events, context)
    if anchor is None:
        return False
    anchor_assessment = context.assessment_by_event_id[anchor.event_id]
    return (
        score >= TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_SCORE_MIN
        and score <= TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_SCORE_MAX
        and _selector_score_decoy_margin(event, context)
        >= TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_DECOY_MARGIN_MIN
        and event.contact_cluster_gain
        >= TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_CLUSTER_MIN
        and assessment.direct_support_score
        >= TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_DIRECT_SUPPORT_MIN
        and assessment.future_preservation_score
        >= TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_FUTURE_PRESERVATION_MIN
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_BLOCKED_FUTURE_MAX
        and event.secondary_structure_compatibility
        >= TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_SECONDARY_STRUCTURE_MIN
        and coupling_nucleus_score(anchor, context)
        >= TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_ANCHOR_SCORE_MIN
        and anchor_assessment.direct_support_score
        >= TRACE_LOOP_REGISTRY_EXTENSION_RESCUE_ANCHOR_DIRECT_SUPPORT_MIN
    )


def _passes_direct_margin_tail_rescue_gate(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> bool:
    assessment = context.assessment_by_event_id[event.event_id]
    score = coupling_nucleus_score(event, context)
    evidence = _direct_constraint_trace_evidence(event, context)
    return (
        score >= TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_SCORE_MIN
        and score <= TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_SCORE_MAX
        and _selector_score_decoy_margin(event, context)
        >= TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_DECOY_MARGIN_MIN
        and event.contact_cluster_gain
        >= TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_CLUSTER_MIN
        and assessment.direct_support_score
        >= TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_DIRECT_SUPPORT_MIN
        and assessment.future_preservation_score
        >= TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_FUTURE_PRESERVATION_MIN
        and assessment.future_preservation_score
        <= TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_FUTURE_PRESERVATION_MAX
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_BLOCKED_FUTURE_MAX
        and context.coupling_decoy_margin_by_event_id[event.event_id]
        >= TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_COUPLING_MARGIN_MIN
        and event.secondary_structure_compatibility
        >= TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_SECONDARY_STRUCTURE_MIN
        and evidence["direct_constraint_count"]
        >= TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_DIRECT_CONSTRAINT_COUNT_MIN
        and evidence["direct_constraint_confidence_sum"]
        >= TRACE_LOOP_DIRECT_MARGIN_TAIL_RESCUE_DIRECT_CONFIDENCE_SUM_MIN
    )


def _passes_terminal_bridge_rescue_gate(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> bool:
    assessment = context.assessment_by_event_id[event.event_id]
    score = coupling_nucleus_score(event, context)
    evidence = _direct_constraint_trace_evidence(event, context)
    return (
        score >= TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_SCORE_MIN
        and score <= TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_SCORE_MAX
        and _selector_score_decoy_margin(event, context)
        >= TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_DECOY_MARGIN_MIN
        and event.contact_cluster_gain >= TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_CLUSTER_MIN
        and assessment.direct_support_score
        >= TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_DIRECT_SUPPORT_MIN
        and assessment.future_preservation_score
        >= TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_FUTURE_PRESERVATION_MIN
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_BLOCKED_FUTURE_MAX
        and context.coupling_decoy_margin_by_event_id[event.event_id]
        >= TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_COUPLING_MARGIN_MIN
        and event.secondary_structure_compatibility
        >= TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_SECONDARY_STRUCTURE_MIN
        and evidence["direct_constraint_count"]
        >= TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_DIRECT_CONSTRAINT_COUNT_MIN
        and evidence["direct_constraint_confidence_sum"]
        >= TRACE_LOOP_TERMINAL_BRIDGE_RESCUE_DIRECT_CONFIDENCE_SUM_MIN
    )


def _passes_high_pressure_tail_rescue_gate(
    event: NucleusClosureEvent,
    context: CouplingNucleusContext,
) -> bool:
    assessment = context.assessment_by_event_id[event.event_id]
    score = coupling_nucleus_score(event, context)
    evidence = _direct_constraint_trace_evidence(event, context)
    return (
        score >= TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_SCORE_MIN
        and score <= TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_SCORE_MAX
        and _selector_score_decoy_margin(event, context)
        >= TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_DECOY_MARGIN_MIN
        and event.contact_cluster_gain
        >= TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_CLUSTER_MIN
        and assessment.direct_support_score
        >= TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_DIRECT_SUPPORT_MIN
        and assessment.future_preservation_score
        >= TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_FUTURE_PRESERVATION_MIN
        and assessment.blocked_future_pressure
        >= TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_BLOCKED_FUTURE_MIN
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_BLOCKED_FUTURE_MAX
        and context.coupling_decoy_margin_by_event_id[event.event_id]
        >= TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_COUPLING_MARGIN_MIN
        and event.secondary_structure_compatibility
        >= TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_SECONDARY_STRUCTURE_MIN
        and evidence["direct_constraint_count"]
        >= TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_DIRECT_CONSTRAINT_COUNT_MIN
        and evidence["direct_constraint_confidence_sum"]
        >= TRACE_LOOP_HIGH_PRESSURE_TAIL_RESCUE_DIRECT_CONFIDENCE_SUM_MIN
    )


def select_coupling_trace_loop_rank_consistent_cluster_gated_events(
    context: CouplingNucleusContext,
) -> tuple[NucleusClosureEvent, ...]:
    recovery_candidates = select_coupling_trace_loop_core_expanded_events(
        context,
        min_contact_cluster_gain=TRACE_LOOP_RANK_CONSISTENT_RECOVERY_CLUSTER_GATE_MIN,
        min_rank_confidence_consistency=TRACE_LOOP_RANK_CONFIDENCE_CONSISTENCY_MIN,
        min_score_confidence_calibration=TRACE_LOOP_SCORE_CONFIDENCE_CALIBRATION_MIN,
        min_rank_length_calibration=TRACE_LOOP_RANK_LENGTH_CALIBRATION_MIN,
    )
    high_confidence_core = tuple(
        event
        for event in recovery_candidates
        if event.contact_cluster_gain >= TRACE_LOOP_RANK_CONSISTENT_CLUSTER_GATE_MIN
    )
    core_by_row = _events_by_row(high_confidence_core)
    recovery_by_row = _events_by_row(recovery_candidates)
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        row_selected = list(core_by_row.get(row.row_id, ()))
        selected_ids = {event.event_id for event in row_selected}
        row_recovery = sorted(
            (
                event
                for event in recovery_by_row.get(row.row_id, ())
                if event.event_id not in selected_ids
            ),
            key=lambda event: (
                -context.assessment_by_event_id[
                    event.event_id
                ].coupling_selectivity_score,
                -context.assessment_by_event_id[event.event_id].direct_support_score,
                -context.coupling_decoy_margin_by_event_id[event.event_id],
                event.segment_a_start,
                event.segment_b_start,
                event.event_id,
            ),
        )
        for event in row_recovery:
            if len(row_selected) >= SELECTED_EVENTS_PER_ROW:
                break
            if not _passes_rank_consistent_recovery_gate(event, context):
                continue
            if any(
                not compatible_future_event(selected_event, event)
                for selected_event in row_selected
            ):
                continue
            row_selected.append(event)
            selected_ids.add(event.event_id)
        selected.extend(row_selected)
    return tuple(selected)


def select_coupling_trace_loop_score_margin_expanded_events(
    context: CouplingNucleusContext,
) -> tuple[NucleusClosureEvent, ...]:
    persistent_events = (
        select_coupling_trace_loop_persistent_rank_consistent_cluster_gated_events(
            context
        )
    )
    persistent_by_row = _events_by_row(persistent_events)
    trace_by_row = _events_by_row(select_coupling_trace_loop_events(context))
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        row_selected = list(persistent_by_row.get(row.row_id, ()))
        selected_ids = {event.event_id for event in row_selected}
        expansion_candidates = sorted(
            (
                event
                for event in trace_by_row.get(row.row_id, ())
                if event.event_id not in selected_ids
                and _passes_score_margin_expansion_gate(event, context)
            ),
            key=lambda event: (
                -coupling_nucleus_score(event, context),
                -_selector_score_decoy_margin(event, context),
                -context.assessment_by_event_id[event.event_id].direct_support_score,
                event.segment_a_start,
                event.segment_b_start,
                event.event_id,
            ),
        )
        for event in expansion_candidates:
            if len(row_selected) >= SELECTED_EVENTS_PER_ROW:
                break
            if any(
                not compatible_future_event(selected_event, event)
                for selected_event in row_selected
            ):
                continue
            row_selected.append(event)
            selected_ids.add(event.event_id)
        selected.extend(row_selected)
    return tuple(selected)


def select_coupling_trace_loop_boundary_continuity_expanded_events(
    context: CouplingNucleusContext,
) -> tuple[NucleusClosureEvent, ...]:
    score_margin_events = select_coupling_trace_loop_score_margin_expanded_events(
        context
    )
    score_margin_by_row = _events_by_row(score_margin_events)
    trace_by_row = _events_by_row(select_coupling_trace_loop_events(context))
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        row_selected = list(score_margin_by_row.get(row.row_id, ()))
        selected_ids = {event.event_id for event in row_selected}
        rescue_candidates = sorted(
            (
                event
                for event in trace_by_row.get(row.row_id, ())
                if event.event_id not in selected_ids
                and _passes_boundary_continuity_rescue_gate(event, context)
            ),
            key=lambda event: (
                -coupling_nucleus_score(event, context),
                -_selector_score_decoy_margin(event, context),
                -context.assessment_by_event_id[event.event_id].direct_support_score,
                -event.secondary_structure_compatibility,
                event.segment_a_start,
                event.segment_b_start,
                event.event_id,
            ),
        )
        for event in rescue_candidates:
            if len(row_selected) >= SELECTED_EVENTS_PER_ROW:
                break
            if any(
                not compatible_future_event(selected_event, event)
                for selected_event in row_selected
            ):
                continue
            row_selected.append(event)
            selected_ids.add(event.event_id)
        selected.extend(row_selected)
    return tuple(selected)


def select_coupling_trace_loop_edge_continuity_expanded_events(
    context: CouplingNucleusContext,
) -> tuple[NucleusClosureEvent, ...]:
    boundary_events = select_coupling_trace_loop_boundary_continuity_expanded_events(
        context
    )
    boundary_by_row = _events_by_row(boundary_events)
    trace_by_row = _events_by_row(select_coupling_trace_loop_events(context))
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        row_selected = list(boundary_by_row.get(row.row_id, ()))
        selected_ids = {event.event_id for event in row_selected}
        rescue_candidates = sorted(
            (
                event
                for event in trace_by_row.get(row.row_id, ())
                if event.event_id not in selected_ids
                and _passes_edge_continuity_rescue_gate(event, context)
            ),
            key=lambda event: (
                -coupling_nucleus_score(event, context),
                -_selector_score_decoy_margin(event, context),
                -context.assessment_by_event_id[event.event_id].direct_support_score,
                -event.contact_cluster_gain,
                -event.secondary_structure_compatibility,
                event.segment_a_start,
                event.segment_b_start,
                event.event_id,
            ),
        )
        for event in rescue_candidates:
            if len(row_selected) >= SELECTED_EVENTS_PER_ROW:
                break
            if any(
                not compatible_future_event(selected_event, event)
                for selected_event in row_selected
            ):
                continue
            row_selected.append(event)
            selected_ids.add(event.event_id)
        selected.extend(row_selected)
    return tuple(selected)


def select_coupling_trace_loop_pressure_release_expanded_events(
    context: CouplingNucleusContext,
) -> tuple[NucleusClosureEvent, ...]:
    edge_events = select_coupling_trace_loop_edge_continuity_expanded_events(context)
    edge_by_row = _events_by_row(edge_events)
    trace_by_row = _events_by_row(select_coupling_trace_loop_events(context))
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        row_selected = list(edge_by_row.get(row.row_id, ()))
        selected_ids = {event.event_id for event in row_selected}
        rescue_candidates = sorted(
            (
                event
                for event in trace_by_row.get(row.row_id, ())
                if event.event_id not in selected_ids
                and _passes_pressure_release_rescue_gate(event, context)
            ),
            key=lambda event: (
                -coupling_nucleus_score(event, context),
                -_selector_score_decoy_margin(event, context),
                -context.assessment_by_event_id[event.event_id].direct_support_score,
                -_direct_constraint_trace_evidence(
                    event,
                    context,
                )["direct_constraint_confidence_sum"],
                event.segment_a_start,
                event.segment_b_start,
                event.event_id,
            ),
        )
        for event in rescue_candidates:
            if len(row_selected) >= SELECTED_EVENTS_PER_ROW:
                break
            if any(
                not compatible_future_event(selected_event, event)
                for selected_event in row_selected
            ):
                continue
            row_selected.append(event)
            selected_ids.add(event.event_id)
        selected.extend(row_selected)
    return tuple(selected)


def select_coupling_trace_loop_registry_extension_expanded_events(
    context: CouplingNucleusContext,
) -> tuple[NucleusClosureEvent, ...]:
    pressure_events = select_coupling_trace_loop_pressure_release_expanded_events(
        context
    )
    pressure_by_row = _events_by_row(pressure_events)
    trace_by_row = _events_by_row(select_coupling_trace_loop_events(context))
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        row_selected = list(pressure_by_row.get(row.row_id, ()))
        selected_ids = {event.event_id for event in row_selected}
        rescue_candidates = sorted(
            (
                event
                for event in trace_by_row.get(row.row_id, ())
                if event.event_id not in selected_ids
                and (
                    _passes_registry_extension_rescue_gate(
                        event,
                        row_selected,
                        context,
                    )
                    or _passes_direct_margin_tail_rescue_gate(event, context)
                )
            ),
            key=lambda event: (
                0
                if _passes_registry_extension_rescue_gate(
                    event,
                    row_selected,
                    context,
                )
                else 1,
                -coupling_nucleus_score(event, context),
                -event.contact_cluster_gain,
                -context.assessment_by_event_id[event.event_id].direct_support_score,
                event.segment_a_start,
                event.segment_b_start,
                event.event_id,
            ),
        )
        for event in rescue_candidates:
            if len(row_selected) >= SELECTED_EVENTS_PER_ROW:
                break
            if any(
                not compatible_future_event(selected_event, event)
                for selected_event in row_selected
            ):
                continue
            row_selected.append(event)
            selected_ids.add(event.event_id)
        selected.extend(row_selected)
    return tuple(selected)


def select_coupling_trace_loop_terminal_bridge_expanded_events(
    context: CouplingNucleusContext,
) -> tuple[NucleusClosureEvent, ...]:
    registry_events = select_coupling_trace_loop_registry_extension_expanded_events(
        context
    )
    registry_by_row = _events_by_row(registry_events)
    trace_by_row = _events_by_row(select_coupling_trace_loop_events(context))
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        row_selected = list(registry_by_row.get(row.row_id, ()))
        selected_ids = {event.event_id for event in row_selected}
        rescue_candidates = sorted(
            (
                event
                for event in trace_by_row.get(row.row_id, ())
                if event.event_id not in selected_ids
                and (
                    _passes_terminal_bridge_rescue_gate(event, context)
                    or _passes_high_pressure_tail_rescue_gate(event, context)
                )
            ),
            key=lambda event: (
                0 if _passes_terminal_bridge_rescue_gate(event, context) else 1,
                -coupling_nucleus_score(event, context),
                -context.assessment_by_event_id[event.event_id].direct_support_score,
                event.segment_a_start,
                event.segment_b_start,
                event.event_id,
            ),
        )
        for event in rescue_candidates:
            if len(row_selected) >= SELECTED_EVENTS_PER_ROW:
                break
            if any(
                not compatible_future_event(selected_event, event)
                for selected_event in row_selected
            ):
                continue
            row_selected.append(event)
            selected_ids.add(event.event_id)
        selected.extend(row_selected)
    return tuple(selected)


def _passes_boundary_field_replacement_gate(
    event: NucleusClosureEvent,
    blockers: Sequence[NucleusClosureEvent],
    context: CouplingNucleusContext,
) -> bool:
    if not blockers:
        return False
    assessment = context.assessment_by_event_id[event.event_id]
    score = coupling_nucleus_score(event, context)
    evidence = _direct_constraint_trace_evidence(event, context)
    blocker_evidence = tuple(
        _direct_constraint_trace_evidence(blocker, context) for blocker in blockers
    )
    max_blocker_confidence_sum = max(
        float(item["direct_constraint_confidence_sum"])
        for item in blocker_evidence
    )
    max_blocker_direct_count = max(
        int(item["direct_constraint_count"]) for item in blocker_evidence
    )
    max_blocker_direct_support = max(
        context.assessment_by_event_id[blocker.event_id].direct_support_score
        for blocker in blockers
    )
    max_blocker_future = max(
        context.assessment_by_event_id[blocker.event_id].future_preservation_score
        for blocker in blockers
    )
    max_blocker_cluster = max(blocker.contact_cluster_gain for blocker in blockers)
    max_blocker_secondary_structure = max(
        blocker.secondary_structure_compatibility for blocker in blockers
    )
    denser_direct_trace = (
        float(evidence["direct_constraint_confidence_sum"])
        - max_blocker_confidence_sum
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_DIRECT_CONFIDENCE_DELTA_MIN
        or int(evidence["direct_constraint_count"]) > max_blocker_direct_count
    )
    return (
        int(evidence["direct_constraint_count"])
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_DIRECT_CONSTRAINT_COUNT_MIN
        and float(evidence["direct_constraint_confidence_sum"])
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_DIRECT_CONFIDENCE_SUM_MIN
        and denser_direct_trace
        and assessment.future_preservation_score
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_PRESERVATION_MIN
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_BLOCKED_FUTURE_MAX
        and event.secondary_structure_compatibility
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_SECONDARY_STRUCTURE_MIN
        and event.contact_cluster_gain
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_CLUSTER_MIN
        and score >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_SCORE_MIN
        and context.coupling_decoy_margin_by_event_id[event.event_id]
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_COUPLING_MARGIN_MIN
        and assessment.direct_support_score
        >= max_blocker_direct_support
        - TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_DIRECT_SUPPORT_DROP_MAX
        and assessment.future_preservation_score
        >= max_blocker_future
        - TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_PRESERVATION_DROP_MAX
        and event.contact_cluster_gain
        >= max_blocker_cluster - TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_CLUSTER_DROP_MAX
        and event.secondary_structure_compatibility
        >= max_blocker_secondary_structure
        - TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_SECONDARY_STRUCTURE_DROP_MAX
    )


def _boundary_field_replacement_span_delta(
    event: NucleusClosureEvent,
    blockers: Sequence[NucleusClosureEvent],
) -> int:
    return max(
        abs(event.sequence_span - blocker.sequence_span) for blocker in blockers
    )


def _passes_boundary_field_future_direct_replacement_gate(
    event: NucleusClosureEvent,
    blockers: Sequence[NucleusClosureEvent],
    context: CouplingNucleusContext,
) -> bool:
    if len(blockers) != 1:
        return False
    assessment = context.assessment_by_event_id[event.event_id]
    state = context.physical_context.state_by_event_id[event.event_id]
    evidence = _direct_constraint_trace_evidence(event, context)
    return (
        _boundary_field_replacement_span_delta(event, blockers)
        <= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_SPAN_CONTINUITY_DELTA_MAX
        and event.normalized_span
        <= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_NORMALIZED_SPAN_MAX
        and coupling_nucleus_score(event, context)
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_SCORE_MIN
        and state.physical_state_score
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_PHYSICAL_MIN
        and state.physical_state_score
        <= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_PHYSICAL_MAX
        and event.contact_cluster_gain
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_CLUSTER_MIN
        and assessment.direct_support_score
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_SUPPORT_MIN
        and assessment.future_preservation_score
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_PRESERVATION_MIN
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_BLOCKED_MAX
        and int(evidence["direct_constraint_count"])
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_COUNT_MIN
        and float(evidence["direct_constraint_confidence_sum"])
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_FUTURE_DIRECT_CONFIDENCE_MIN
    )


def _passes_boundary_field_physical_release_replacement_gate(
    event: NucleusClosureEvent,
    blockers: Sequence[NucleusClosureEvent],
    context: CouplingNucleusContext,
) -> bool:
    if len(blockers) != 1:
        return False
    assessment = context.assessment_by_event_id[event.event_id]
    state = context.physical_context.state_by_event_id[event.event_id]
    evidence = _direct_constraint_trace_evidence(event, context)
    return (
        _boundary_field_replacement_span_delta(event, blockers)
        <= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_SPAN_CONTINUITY_DELTA_MAX
        and coupling_nucleus_score(event, context)
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_SCORE_MIN
        and state.physical_state_score
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_STATE_MIN
        and state.burial_gain
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_BURIAL_MIN
        and event.contact_cluster_gain
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_CLUSTER_MIN
        and event.secondary_structure_compatibility
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_SECONDARY_STRUCTURE_MIN
        and assessment.future_preservation_score
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_FUTURE_MIN
        and assessment.blocked_future_pressure
        <= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_BLOCKED_MAX
        and int(evidence["direct_constraint_count"])
        >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_PHYSICAL_RELEASE_DIRECT_COUNT_MIN
    )


def _boundary_field_replacement_mode(
    event: NucleusClosureEvent,
    blockers: Sequence[NucleusClosureEvent],
    context: CouplingNucleusContext,
) -> str | None:
    if _passes_boundary_field_future_direct_replacement_gate(
        event,
        blockers,
        context,
    ):
        return "future_direct"
    if _passes_boundary_field_physical_release_replacement_gate(
        event,
        blockers,
        context,
    ):
        return "physical_release"
    return None


def _boundary_field_replacement_sort_key(
    mode: str,
    event: NucleusClosureEvent,
    blockers: Sequence[NucleusClosureEvent],
    context: CouplingNucleusContext,
) -> tuple[object, ...]:
    assessment = context.assessment_by_event_id[event.event_id]
    state = context.physical_context.state_by_event_id[event.event_id]
    evidence = _direct_constraint_trace_evidence(event, context)
    span_delta = _boundary_field_replacement_span_delta(event, blockers)
    if mode == "future_direct":
        return (
            0,
            span_delta,
            -float(evidence["direct_constraint_confidence_sum"]),
            -assessment.future_preservation_score,
            -coupling_nucleus_score(event, context),
            event.segment_a_start,
            event.segment_b_start,
            event.event_id,
        )
    physical_integrity = (
        state.physical_state_score
        + state.burial_gain
        + event.contact_cluster_gain
        - state.unsatisfied_polar_penalty
    )
    return (
        1,
        span_delta,
        -physical_integrity,
        -assessment.future_preservation_score,
        -coupling_nucleus_score(event, context),
        event.segment_a_start,
        event.segment_b_start,
        event.event_id,
    )


def select_coupling_trace_loop_boundary_field_replacement_probe_events(
    context: CouplingNucleusContext,
) -> tuple[NucleusClosureEvent, ...]:
    terminal_events = select_coupling_trace_loop_terminal_bridge_expanded_events(
        context
    )
    terminal_by_row = _events_by_row(terminal_events)
    competitive_by_row = _events_by_row(context.competitive_events)
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        terminal_row_selected = list(terminal_by_row.get(row.row_id, ()))
        terminal_selected_ids = {
            event.event_id for event in terminal_row_selected
        }
        replacement_candidates: list[
            tuple[str, NucleusClosureEvent, tuple[NucleusClosureEvent, ...]]
        ] = []
        for event in competitive_by_row.get(row.row_id, ()):
            if event.event_id in terminal_selected_ids:
                continue
            blockers = [
                selected_event
                for selected_event in terminal_row_selected
                if not compatible_future_event(selected_event, event)
            ]
            if not blockers:
                continue
            blocker_ids = {blocker.event_id for blocker in blockers}
            if not blocker_ids <= terminal_selected_ids:
                continue
            remaining = [
                selected_event
                for selected_event in terminal_row_selected
                if selected_event.event_id
                not in blocker_ids
            ]
            if any(
                not compatible_future_event(selected_event, event)
                for selected_event in remaining
            ):
                continue
            mode = _boundary_field_replacement_mode(
                event,
                tuple(blockers),
                context,
            )
            if mode is not None:
                replacement_candidates.append((mode, event, tuple(blockers)))

        chosen_replacements: list[
            tuple[str, NucleusClosureEvent, tuple[NucleusClosureEvent, ...]]
        ] = []
        replaced_terminal_ids: set[str] = set()
        replacement_candidates = sorted(
            replacement_candidates,
            key=lambda item: _boundary_field_replacement_sort_key(
                item[0],
                item[1],
                item[2],
                context,
            ),
        )
        for _mode, event, blockers in replacement_candidates:
            if (
                len(chosen_replacements)
                >= TRACE_LOOP_BOUNDARY_FIELD_REPLACEMENT_LIMIT_PER_ROW
            ):
                break
            blocker_ids = {blocker.event_id for blocker in blockers}
            if replaced_terminal_ids & blocker_ids:
                continue
            if any(
                not compatible_future_event(replacement_event, event)
                for _, replacement_event, _ in chosen_replacements
            ):
                continue
            chosen_replacements.append((mode, event, blockers))
            replaced_terminal_ids.update(blocker_ids)
        row_selected = [
            event
            for event in terminal_row_selected
            if event.event_id not in replaced_terminal_ids
        ]
        row_selected.extend(event for _, event, _ in chosen_replacements)
        selected.extend(row_selected)
    return tuple(selected)


def select_coupling_trace_loop_macro_scale_future_preserved_events(
    context: CouplingNucleusContext,
) -> tuple[NucleusClosureEvent, ...]:
    selected = select_coupling_trace_loop_boundary_field_replacement_probe_events(
        context
    )
    return tuple(
        event
        for event in selected
        if context.assessment_by_event_id[
            event.event_id
        ].future_preservation_score
        >= TRACE_LOOP_MACRO_SCALE_FUTURE_PRESERVATION_MIN
    )


def select_coupling_trace_loop_persistent_rank_consistent_cluster_gated_events(
    context: CouplingNucleusContext,
) -> tuple[NucleusClosureEvent, ...]:
    recovery_candidates = select_coupling_trace_loop_core_expanded_events(
        context,
        min_contact_cluster_gain=TRACE_LOOP_RANK_CONSISTENT_RECOVERY_CLUSTER_GATE_MIN,
        min_rank_confidence_consistency=TRACE_LOOP_RANK_CONFIDENCE_CONSISTENCY_MIN,
        min_score_confidence_calibration=TRACE_LOOP_SCORE_CONFIDENCE_CALIBRATION_MIN,
        min_rank_length_calibration=TRACE_LOOP_RANK_LENGTH_CALIBRATION_MIN,
    )
    high_confidence_core = tuple(
        event
        for event in recovery_candidates
        if event.contact_cluster_gain >= TRACE_LOOP_RANK_CONSISTENT_CLUSTER_GATE_MIN
    )
    core_by_row = _events_by_row(high_confidence_core)
    recovery_by_row = _events_by_row(recovery_candidates)
    selected: list[NucleusClosureEvent] = []
    for row in context.rows:
        row_selected = list(core_by_row.get(row.row_id, ()))
        selected_ids = {event.event_id for event in row_selected}
        row_recovery_pool = tuple(recovery_by_row.get(row.row_id, ()))
        row_recovery = sorted(
            (
                event
                for event in row_recovery_pool
                if event.event_id not in selected_ids
            ),
            key=lambda event: (
                -context.assessment_by_event_id[
                    event.event_id
                ].coupling_selectivity_score,
                -context.assessment_by_event_id[event.event_id].direct_support_score,
                -context.coupling_decoy_margin_by_event_id[event.event_id],
                event.segment_a_start,
                event.segment_b_start,
                event.event_id,
            ),
        )
        for event in row_recovery:
            if len(row_selected) >= SELECTED_EVENTS_PER_ROW:
                break
            if not (
                _passes_rank_consistent_recovery_gate(event, context)
                or _passes_persistent_recovery_gate(
                    event,
                    context,
                    row_recovery_pool,
                )
            ):
                continue
            if any(
                not compatible_future_event(selected_event, event)
                for selected_event in row_selected
            ):
                continue
            row_selected.append(event)
            selected_ids.add(event.event_id)
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
        possible_region_pair_count = sum(
            len(event.candidate_region_pairs()) for event in row_events
        )
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
    selected_scores = [
        comparison.real_coupling_selectivity_score for comparison in comparisons
    ]
    decoy_scores = [
        comparison.decoy_coupling_selectivity_score for comparison in comparisons
    ]
    selected_selectivity_mean = _rounded(
        mean(selected_scores) if selected_scores else 0.0
    )
    decoy_selectivity_mean = _rounded(
        mean(decoy_scores) if decoy_scores else 0.0
    )
    selectivity_margin_mean = _rounded(
        selected_selectivity_mean - decoy_selectivity_mean
    )
    selected_nucleus_scores = [
        coupling_nucleus_score(
            context.event_by_id[comparison.real_event_id],
            context,
        )
        for comparison in comparisons
    ]
    decoy_nucleus_scores = [
        coupling_nucleus_score(
            context.event_by_id[comparison.decoy_event_id],
            context,
        )
        for comparison in comparisons
    ]
    selected_nucleus_mean = _rounded(
        mean(selected_nucleus_scores) if selected_nucleus_scores else 0.0
    )
    decoy_nucleus_mean = _rounded(
        mean(decoy_nucleus_scores) if decoy_nucleus_scores else 0.0
    )
    nucleus_score_enrichment = _rounded(
        selected_nucleus_mean / decoy_nucleus_mean if decoy_nucleus_mean else 0.0
    )
    nucleus_score_beats = _rounded(
        sum(
            1
            for selected_score, decoy_score in zip(
                selected_nucleus_scores,
                decoy_nucleus_scores,
            )
            if selected_score > decoy_score
        )
        / len(selected_nucleus_scores)
        if selected_nucleus_scores
        else 0.0
    )
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
        mean_selected_coupling_selectivity_score=selected_selectivity_mean,
        mean_decoy_coupling_selectivity_score=decoy_selectivity_mean,
        mean_coupling_decoy_selectivity_margin=selectivity_margin_mean,
        mean_coupling_nucleus_score=selected_nucleus_mean,
        mean_decoy_coupling_nucleus_score=decoy_nucleus_mean,
        mean_coupling_nucleus_decoy_margin=_rounded(
            selected_nucleus_mean - decoy_nucleus_mean
        ),
        real_vs_decoy_coupling_nucleus_enrichment_ratio=(
            nucleus_score_enrichment
        ),
        real_beats_decoy_coupling_nucleus_score_rate=nucleus_score_beats,
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
            adaptive_profile = _adaptive_coupling_floor_profile(event, context)
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
                    "adaptive_phase_mode": adaptive_profile.phase_mode,
                    "adaptive_sequence_complexity": (
                        adaptive_profile.sequence_complexity
                    ),
                    "adaptive_coupling_depth_over_length": (
                        adaptive_profile.coupling_depth_over_length
                    ),
                    "adaptive_target_coverage": adaptive_profile.target_coverage,
                    "adaptive_future_preservation_floor": (
                        adaptive_profile.future_preservation_floor
                    ),
                    "adaptive_physical_score_floor": (
                        adaptive_profile.physical_score_floor
                    ),
                    "adaptive_gate_enabled": (
                        adaptive_profile.adaptive_gate_enabled
                    ),
                    "adaptive_low_signal_rescue_enabled": (
                        adaptive_profile.low_signal_rescue_enabled
                    ),
                    "adaptive_signal_reason": adaptive_profile.signal_reason,
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


def contact_collapse_results_for_selected_events(
    context: CouplingNucleusContext,
    selected_events: Sequence[NucleusClosureEvent],
    *,
    collapse_strategy: str = PRIMARY_CONTACT_COLLAPSE_STRATEGY,
) -> tuple[RowCollapseResult, ...]:
    """Collapse selected 8x8 event regions into residue-pair contact maps.

    This is part of the selector pipeline, not an oracle filter.  Native labels are
    attached only inside RowCollapseEvaluation after the sequence/coupling-only
    pair subset has already been selected.
    """
    row_by_id = {row.row_id: row for row in context.rows}
    events_by_row = _events_by_row(selected_events)
    constraints_by_row = context.coupling_dataset.constraints_by_row_id()
    features_by_row = feature_rows_by_row_id(contact_law_feature_rows(context.rows))
    results: list[RowCollapseResult] = []
    for row_id in sorted(events_by_row):
        row = row_by_id.get(row_id)
        if row is None:
            continue
        row_events = tuple(events_by_row[row_id])
        if not row_events:
            continue
        if collapse_strategy == SELF_DECIDING_STRATEGY_NAME:
            # The self-deciding collapse strategy derives its own per-event cutoff
            # from the score distribution, long-range candidate space, sequence-
            # inferred phase shape, direct-coupling roots, and gap clarity. These
            # placeholders are intentionally zero so no fixed event budget leaks
            # into the decision path.
            min_pairs_per_event = 0
            max_pairs_per_event = 0
        elif collapse_strategy == "frontier_internal_gap_balanced":
            min_pairs_per_event = 1
            max_pairs_per_event = DEFAULT_BALANCED_PAIRS_PER_EVENT
        else:
            min_pairs_per_event = DEFAULT_BALANCED_PAIRS_PER_EVENT
            max_pairs_per_event = DEFAULT_BALANCED_PAIRS_PER_EVENT
        results.append(
            collapse_row_event_regions(
                row=row,
                events=row_events,
                row_features=features_by_row.get(row_id, ()),
                row_constraints=constraints_by_row.get(row_id, ()),
                collapse_strategy=collapse_strategy,
                min_pairs_per_event=min_pairs_per_event,
                max_pairs_per_event=max_pairs_per_event,
            )
        )
    return tuple(results)


def _row_collapse_report_row(result: RowCollapseResult) -> dict[str, object]:
    collapsed_pairs = tuple(pair.pair() for pair in result.collapsed_pairs)
    return {
        **result.evaluation.to_dict(),
        "collapsed_contact_map_hash": contact_map_hash(collapsed_pairs),
        "uncollapsed_region_contact_map_hash": contact_map_hash(result.uncollapsed_pairs),
    }


def contact_collapse_row_rows(
    results: Sequence[RowCollapseResult],
) -> list[dict[str, object]]:
    return [_row_collapse_report_row(result) for result in results]


def contact_collapse_pair_rows(
    results: Sequence[RowCollapseResult],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        rows.extend(pair.to_dict() for pair in result.collapsed_pairs)
    return rows


def contact_collapse_event_rows(
    results: Sequence[RowCollapseResult],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        rows.extend(summary.to_dict() for summary in result.event_summaries)
    return rows


def contact_collapse_summary(
    results: Sequence[RowCollapseResult],
) -> dict[str, object]:
    row_reports = contact_collapse_row_rows(results)
    one_cll = next(
        (row for row in row_reports if row.get("source_accession") == "1CLL:A"),
        None,
    )
    return {
        "contact_collapse_integrated": True,
        "contact_collapse_kind": EVENT_REGION_CONTACT_COLLAPSE_KIND,
        "contact_collapse_boundary": EVENT_REGION_CONTACT_COLLAPSE_BOUNDARY,
        "contact_collapse_source_selector": PRIMARY_CONTACT_COLLAPSE_SELECTOR_NAME,
        "contact_collapse_strategy": PRIMARY_CONTACT_COLLAPSE_STRATEGY,
        "contact_collapse_row_count": len(row_reports),
        "contact_collapse_collapsed_pair_count": sum(
            int(row["collapsed_pair_count"]) for row in row_reports
        ),
        "contact_collapse_mean_contact_precision": _rounded(
            mean([float(row["collapsed_contact_precision"]) for row in row_reports])
            if row_reports
            else 0.0
        ),
        "contact_collapse_mean_long_range_precision": _rounded(
            mean([float(row["collapsed_long_range_precision"]) for row in row_reports])
            if row_reports
            else 0.0
        ),
        "contact_collapse_mean_long_range_recall": _rounded(
            mean([float(row["collapsed_long_range_recall"]) for row in row_reports])
            if row_reports
            else 0.0
        ),
        "contact_collapse_mean_long_range_f1": _rounded(
            mean([float(row["collapsed_long_range_f1"]) for row in row_reports])
            if row_reports
            else 0.0
        ),
        "contact_collapse_native_truth_used_before_selection": False,
        "contact_collapse_coordinate_truth_used_before_selection": False,
        "contact_collapse_native_truth_attached_after_selection_for_evaluation": True,
        "contact_collapse_1cll_self_deciding": one_cll or {},
        "contact_collapse_1cll_balanced": one_cll or {},
        "contact_collapse_1cll_internal_gap": one_cll or {},
        "contact_collapse_rows": row_reports,
    }


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
    contact_collapse_results: Sequence[RowCollapseResult] = (),
) -> dict[str, object]:
    selector_lookup = {row.selector_name: row for row in selector_rows}
    coupling_target_survives = any(row.survives_targets for row in selector_rows)
    claim_mode_failures = coupling_claim_mode_validation_failures(
        context.coupling_dataset
    )
    claim_mode_validation_passed = not claim_mode_failures
    claim_allowed = coupling_target_survives and claim_mode_validation_passed
    contact_collapse_payload = contact_collapse_summary(contact_collapse_results)
    return {
        **contact_collapse_payload,
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
        "contact_collapse_integrated": report.get("contact_collapse_integrated", False),
        "contact_collapse_strategy": report.get("contact_collapse_strategy", ""),
        "contact_collapse_mean_contact_precision": report.get(
            "contact_collapse_mean_contact_precision",
            0.0,
        ),
        "contact_collapse_mean_long_range_f1": report.get(
            "contact_collapse_mean_long_range_f1",
            0.0,
        ),
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
    contact_collapse_results = contact_collapse_results_for_selected_events(
        context,
        selections[PRIMARY_CONTACT_COLLAPSE_SELECTOR_NAME],
    )
    report = build_coupling_nucleus_selector_report(
        context=context,
        selector_rows=selector_rows,
        source_benchmark_file=benchmark_file,
        coupling_file=coupling_file,
        contact_collapse_results=contact_collapse_results,
    )
    outputs = write_coupling_nucleus_selector_outputs(
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
    write_csv_rows(
        contact_collapse_row_rows(contact_collapse_results),
        selected_events_path.with_name("coupling_nucleus_selector_contact_collapse_rows.csv"),
    )
    write_csv_rows(
        contact_collapse_pair_rows(contact_collapse_results),
        selected_events_path.with_name("coupling_nucleus_selector_collapsed_contacts.csv"),
    )
    write_csv_rows(
        contact_collapse_event_rows(contact_collapse_results),
        selected_events_path.with_name("coupling_nucleus_selector_contact_collapse_events.csv"),
    )
    return outputs


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
