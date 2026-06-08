from __future__ import annotations

"""Dense native-free geometry-field annealing challenge.

This layer tests the next hypothesis after the native-free geometry-field V0
result: can sparse external DCA anchors be amplified into a dense, accurate
contact field without AlphaFold/ESMFold, templates, a trained structure model,
or native/coordinate truth before selection?

The implementation is intentionally bounded and deterministic.  It has three
lanes:

* pure_sequence_coarse_annealing: sequence priors + coarse chain relaxation only.
* external_dca_compact_annealing: DCA anchors pull a coarse C-alpha trace into a
  compact geometry field.
* external_dca_multifield_annealing: union of compact annealing and a long-range
  anchor-temperature lane.  This explicitly tests whether we can recover the
  missing long-range contacts by increasing exploration pressure.

Native coordinates/native contacts are used only after selection for audit and
matched controls.  The claim gate stays closed unless every row clears strict
precision/recall/long-range gates.
"""

import hashlib
from collections import defaultdict
from dataclasses import asdict, dataclass
from math import cos, exp, pi, sin, sqrt
from statistics import mean
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_five_axis_physics import matched_control_pairs
from pharmacotopology.folding_global_factor_graph_ensemble import (
    GlobalContactFactor,
    GlobalFactorGraphSolution,
    build_contact_factors,
    solve_top_k_factor_graph,
)
from pharmacotopology.folding_native_contact_eval import (
    ContactMetricPacket,
    ContactPair,
    contact_map_hash,
    evaluate_contact_prediction,
    normalized_contact_pairs,
)
from pharmacotopology.folding_native_free_geometry_field import (
    EXTERNAL_DCA_GEOMETRY_FIELD_MODE,
    PURE_SEQUENCE_SYMBOLIC_GEOMETRY_MODE,
    _consensus_pairs,
    _top_pair_pool,
    build_native_free_geometry_scores,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    MIN_SEQUENCE_SEPARATION,
    RealCoordinateVisualRow,
)
from pharmacotopology.folding_sequence_physical_priors import predict_lightweight_secondary_structure


DENSE_GEOMETRY_ANNEALING_KIND = "dense_native_free_geometry_annealing_v0"
DENSE_GEOMETRY_ANNEALING_DECISION_KIND = "dense_geometry_annealing_contact_decision_v0"
PURE_SEQUENCE_COARSE_ANNEALING_MODE = "pure_sequence_coarse_annealing"
EXTERNAL_DCA_COMPACT_ANNEALING_MODE = "external_dca_compact_annealing"
EXTERNAL_DCA_MULTIFIELD_ANNEALING_MODE = "external_dca_multifield_annealing"
COMPACT_CHANNEL = "compact_distance_geometry_channel"
LONG_RANGE_CHANNEL = "long_range_anchor_temperature_channel"
MULTIFIELD_CHANNEL = "compact_plus_long_range_union_channel"
DENSE_GEOMETRY_ANNEALING_RULE = (
    "safe_dca_anchors_or_sequence_only;bounded_coarse_chain_relaxation;"
    "secondary_structure_local_springs;steric_repulsion;compactness_pull;"
    "factor_graph_contact_patch_selection;optional_long_range_temperature_union;"
    "native_audit_after_selection_only"
)
DENSE_GEOMETRY_ANNEALING_CLAIM_RULE = (
    "claim_requires_no_coordinate_truth_no_native_truth_no_structure_model_no_learned_geometry;"
    "all_rows_precision_ge_0_70_and_recall_ge_0_70_and_long_range_recall_ge_0_70;"
    "matched_control_f1_and_long_range_margin_ge_0_15;"
    "pure_sequence_mode_required_for_universal_physical_law_claim"
)


@dataclass(frozen=True)
class DenseGeometryContactDecision:
    kind: str
    row_id: str
    source_accession: str
    source_mode: str
    selection_channel: str
    i: int
    j: int
    sequence_separation: int
    base_score: float
    compact_geometry_score: float
    long_range_temperature_score: float
    final_score: float
    selected: bool
    selected_by_compact_channel: bool
    selected_by_long_range_channel: bool
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False
    msa_dca_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def pair(self) -> ContactPair:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class DenseGeometryRowReport:
    row_id: str
    source_accession: str
    source_mode: str
    sequence_hash: str
    sequence_length: int
    safe_anchor_count: int
    relaxation_step_count: int
    compact_channel_selected_count: int
    long_range_channel_selected_count: int
    selected_contact_count: int
    selected_long_range_contact_count: int
    selected_contact_map_hash: str
    compact_channel_metric_after_native_audit: ContactMetricPacket
    long_range_channel_metric_after_native_audit: ContactMetricPacket
    metric_after_native_audit: ContactMetricPacket
    matched_control_count: int
    best_control_f1_after_audit: float
    best_control_long_range_recall_after_audit: float
    f1_margin_vs_best_control: float
    long_range_recall_margin_vs_best_control: float
    row_dense_geometry_claim_allowed: bool
    row_universal_physical_law_claim_allowed: bool
    row_claim_rejection_reason: str
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    native_truth_attached_after_selection_for_evaluation: bool = True
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False
    msa_dca_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["compact_channel_metric_after_native_audit"] = self.compact_channel_metric_after_native_audit.to_dict()
        payload["long_range_channel_metric_after_native_audit"] = self.long_range_channel_metric_after_native_audit.to_dict()
        payload["metric_after_native_audit"] = self.metric_after_native_audit.to_dict()
        return payload


@dataclass(frozen=True)
class DenseGeometryAnnealingPacket:
    kind: str
    source_mode: str
    row_count: int
    decision_rule: str
    claim_rule: str
    pure_sequence_mode_included: bool
    external_dca_anchor_field_included: bool
    coarse_annealing_included: bool
    long_range_temperature_channel_included: bool
    global_factor_graph_included: bool
    mean_native_contact_precision_after_audit: float
    mean_native_contact_recall_after_audit: float
    mean_long_range_contact_recall_after_audit: float
    mean_contact_map_f1_after_audit: float
    mean_compact_channel_f1_after_audit: float
    mean_compact_channel_long_range_recall_after_audit: float
    mean_long_range_channel_f1_after_audit: float
    mean_long_range_channel_long_range_recall_after_audit: float
    min_row_precision_after_audit: float
    min_row_recall_after_audit: float
    min_row_long_range_recall_after_audit: float
    mean_f1_margin_vs_best_control: float
    mean_long_range_recall_margin_vs_best_control: float
    dense_geometry_annealing_claim_allowed: bool
    universal_physical_law_claim_allowed: bool
    folding_problem_solved: bool
    claim_rejection_reason: str
    rows: tuple[DenseGeometryRowReport, ...]
    factors: tuple[GlobalContactFactor, ...]
    solutions: tuple[GlobalFactorGraphSolution, ...]
    decisions: tuple[DenseGeometryContactDecision, ...]
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    structure_model_used_before_selection: bool = False
    learned_geometry_prior_used_before_selection: bool = False
    msa_dca_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "source_mode": self.source_mode,
            "row_count": self.row_count,
            "decision_rule": self.decision_rule,
            "claim_rule": self.claim_rule,
            "pure_sequence_mode_included": self.pure_sequence_mode_included,
            "external_dca_anchor_field_included": self.external_dca_anchor_field_included,
            "coarse_annealing_included": self.coarse_annealing_included,
            "long_range_temperature_channel_included": self.long_range_temperature_channel_included,
            "global_factor_graph_included": self.global_factor_graph_included,
            "mean_native_contact_precision_after_audit": self.mean_native_contact_precision_after_audit,
            "mean_native_contact_recall_after_audit": self.mean_native_contact_recall_after_audit,
            "mean_long_range_contact_recall_after_audit": self.mean_long_range_contact_recall_after_audit,
            "mean_contact_map_f1_after_audit": self.mean_contact_map_f1_after_audit,
            "mean_compact_channel_f1_after_audit": self.mean_compact_channel_f1_after_audit,
            "mean_compact_channel_long_range_recall_after_audit": self.mean_compact_channel_long_range_recall_after_audit,
            "mean_long_range_channel_f1_after_audit": self.mean_long_range_channel_f1_after_audit,
            "mean_long_range_channel_long_range_recall_after_audit": self.mean_long_range_channel_long_range_recall_after_audit,
            "min_row_precision_after_audit": self.min_row_precision_after_audit,
            "min_row_recall_after_audit": self.min_row_recall_after_audit,
            "min_row_long_range_recall_after_audit": self.min_row_long_range_recall_after_audit,
            "mean_f1_margin_vs_best_control": self.mean_f1_margin_vs_best_control,
            "mean_long_range_recall_margin_vs_best_control": self.mean_long_range_recall_margin_vs_best_control,
            "dense_geometry_annealing_claim_allowed": self.dense_geometry_annealing_claim_allowed,
            "universal_physical_law_claim_allowed": self.universal_physical_law_claim_allowed,
            "folding_problem_solved": self.folding_problem_solved,
            "claim_rejection_reason": self.claim_rejection_reason,
            "rows": [row.to_dict() for row in self.rows],
            "factors": [factor.to_dict() for factor in self.factors],
            "solutions": [solution.to_dict() for solution in self.solutions],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "coordinate_truth_used_before_selection": self.coordinate_truth_used_before_selection,
            "native_truth_used_before_selection": self.native_truth_used_before_selection,
            "structure_model_used_before_selection": self.structure_model_used_before_selection,
            "learned_geometry_prior_used_before_selection": self.learned_geometry_prior_used_before_selection,
            "msa_dca_used_before_selection": self.msa_dca_used_before_selection,
            "raw_sequence_exposed": self.raw_sequence_exposed,
        }


def _score(value: float) -> float:
    return round(float(value), 6)


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, float(value)))


def _rounded(value: float) -> float:
    return round(_clamp(value), 6)


def _mean(values: Sequence[float]) -> float:
    return _rounded(mean(values)) if values else 0.0


def _pair_hash(pairs: Iterable[ContactPair]) -> str:
    return contact_map_hash(normalized_contact_pairs(pairs))


def _safe_constraints_for_row(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    source_mode: str,
) -> tuple[tuple[CouplingConstraint, ...], bool, bool, bool, bool]:
    if source_mode == PURE_SEQUENCE_COARSE_ANNEALING_MODE:
        return (), False, False, False, False
    safe: list[CouplingConstraint] = []
    coordinate_taint = False
    native_taint = False
    structure_model = False
    for constraint in constraints:
        if constraint.source_accession != row.source_accession and constraint.row_id != row.row_id:
            continue
        if constraint.coordinate_truth_used_to_build_constraint:
            coordinate_taint = True
            continue
        if constraint.native_truth_used_before_coupling_selection:
            native_taint = True
            continue
        if constraint.structure_model_used:
            structure_model = True
            continue
        if constraint.sequence_separation < MIN_SEQUENCE_SEPARATION:
            continue
        safe.append(constraint)
    safe.sort(key=lambda item: (-(float(item.confidence)), item.rank if item.rank else 999999, item.i, item.j))
    return tuple(safe), coordinate_taint, native_taint, structure_model, bool(safe)


def _initial_coarse_trace(row: RealCoordinateVisualRow) -> list[list[float]]:
    ss = predict_lightweight_secondary_structure(row.sequence)
    coords: list[list[float]] = []
    x = y = z = 0.0
    angle = 0.0
    for idx, label in enumerate(ss):
        if label == "H":
            angle += 100.0 * pi / 180.0
            x += 1.55
            y += 2.10 * cos(angle)
            z += 2.10 * sin(angle)
        elif label == "E":
            x += 3.30
            y += 1.10 * (-1.0 if idx % 2 else 1.0)
            z += 0.22 * sin(idx * 0.6)
        else:
            angle += 47.0 * pi / 180.0
            x += 2.35
            y += 1.45 * cos(angle)
            z += 1.45 * sin(angle)
        coords.append([x, y, z])
    return coords


def _relax_coarse_trace(
    *,
    row: RealCoordinateVisualRow,
    anchors: Sequence[CouplingConstraint],
    steps: int,
    max_anchors: int = 72,
) -> list[list[float]]:
    coords = _initial_coarse_trace(row)
    n = row.sequence_length
    ss = predict_lightweight_secondary_structure(row.sequence)
    edges: list[tuple[int, int, float, float]] = []
    for left in range(1, n):
        edges.append((left - 1, left, 3.80, 0.55))
    for index, label in enumerate(ss, start=1):
        if label == "H":
            for offset, target, weight in ((3, 5.20, 0.10), (4, 6.20, 0.10)):
                if index + offset <= n:
                    edges.append((index - 1, index + offset - 1, target, weight))
        elif label == "E" and index + 2 <= n:
            edges.append((index - 1, index + 1, 6.80, 0.04))
    for constraint in anchors[:max_anchors]:
        confidence = float(constraint.confidence)
        left = constraint.i
        right = constraint.j
        target = 7.40 if abs(right - left) >= 8 else 6.20
        edges.append((left - 1, right - 1, target, 0.75 * confidence))
        for delta in (-1, 1):
            if 1 <= left + delta <= n:
                edges.append((left + delta - 1, right - 1, target + 0.30, 0.22 * confidence))
            if 1 <= right + delta <= n:
                edges.append((left - 1, right + delta - 1, target + 0.30, 0.22 * confidence))

    step_count = max(8, int(steps))
    for step in range(step_count):
        learning_rate = 0.030 * (1.0 - step / (step_count * 1.2))
        for left, right, target, weight in edges:
            lx, ly, lz = coords[left]
            rx, ry, rz = coords[right]
            dx = rx - lx
            dy = ry - ly
            dz = rz - lz
            distance = sqrt(dx * dx + dy * dy + dz * dz) + 1e-6
            force = learning_rate * weight * (distance - target) / distance
            mx = dx * force * 0.5
            my = dy * force * 0.5
            mz = dz * force * 0.5
            coords[left][0] += mx
            coords[left][1] += my
            coords[left][2] += mz
            coords[right][0] -= mx
            coords[right][1] -= my
            coords[right][2] -= mz

        center_x = sum(point[0] for point in coords) / n
        center_y = sum(point[1] for point in coords) / n
        center_z = sum(point[2] for point in coords) / n
        radius = 2.2 * (n ** 0.38)
        for point in coords:
            point[0] -= center_x * 0.020
            point[1] -= center_y * 0.020
            point[2] -= center_z * 0.020
            distance = sqrt(point[0] * point[0] + point[1] * point[1] + point[2] * point[2]) + 1e-6
            if distance > radius:
                force = learning_rate * 0.035 * (distance - radius) / distance
                point[0] -= point[0] * force
                point[1] -= point[1] * force
                point[2] -= point[2] * force

        if step % 4 == 0:
            for left in range(0, n, 3):
                for right in range(left + 6, n, 5):
                    dx = coords[right][0] - coords[left][0]
                    dy = coords[right][1] - coords[left][1]
                    dz = coords[right][2] - coords[left][2]
                    distance = sqrt(dx * dx + dy * dy + dz * dz) + 1e-6
                    if distance >= 3.20:
                        continue
                    force = learning_rate * 0.12 * (3.20 - distance) / distance
                    mx = dx * force * 0.5
                    my = dy * force * 0.5
                    mz = dz * force * 0.5
                    coords[left][0] -= mx
                    coords[left][1] -= my
                    coords[left][2] -= mz
                    coords[right][0] += mx
                    coords[right][1] += my
                    coords[right][2] += mz
    return coords


def _compact_geometry_scores(row: RealCoordinateVisualRow, coords: Sequence[Sequence[float]]) -> dict[ContactPair, float]:
    scores: dict[ContactPair, float] = {}
    for left in range(1, row.sequence_length + 1):
        lx, ly, lz = coords[left - 1]
        for right in range(left + MIN_SEQUENCE_SEPARATION, row.sequence_length + 1):
            rx, ry, rz = coords[right - 1]
            distance = sqrt((rx - lx) ** 2 + (ry - ly) ** 2 + (rz - lz) ** 2)
            if distance <= 7.8:
                value = 1.0 - 0.025 * max(0.0, distance - 5.0)
            else:
                value = exp(-(distance - 7.8) / 4.2)
            scores[(left, right)] = _rounded(value)
    return scores


def _base_scores_for_mode(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    source_mode: str,
) -> tuple[dict[ContactPair, float], bool]:
    if source_mode == PURE_SEQUENCE_COARSE_ANNEALING_MODE:
        scores, _components, _anchor_count, _coordinate_taint, _native_taint, _structure_model, _dca_used = build_native_free_geometry_scores(
            row=row,
            constraints=(),
            source_mode=PURE_SEQUENCE_SYMBOLIC_GEOMETRY_MODE,
        )
        return dict(scores), False
    scores, _components, _anchor_count, _coordinate_taint, _native_taint, _structure_model, dca_used = build_native_free_geometry_scores(
        row=row,
        constraints=constraints,
        source_mode=EXTERNAL_DCA_GEOMETRY_FIELD_MODE,
    )
    return dict(scores), dca_used


def _compose_compact_scores(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    source_mode: str,
    relaxation_steps: int,
) -> tuple[dict[ContactPair, float], dict[ContactPair, float], bool]:
    safe_constraints, _coordinate_taint, _native_taint, _structure_model, dca_used = _safe_constraints_for_row(
        row=row,
        constraints=constraints,
        source_mode=source_mode,
    )
    base_scores, base_dca_used = _base_scores_for_mode(row=row, constraints=constraints, source_mode=source_mode)
    coords = _relax_coarse_trace(row=row, anchors=safe_constraints, steps=relaxation_steps)
    compact_scores = _compact_geometry_scores(row, coords)
    composed: dict[ContactPair, float] = {}
    for pair, base in base_scores.items():
        compact = compact_scores.get(pair, 0.0)
        if source_mode == PURE_SEQUENCE_COARSE_ANNEALING_MODE:
            value = 0.64 * compact + 0.36 * base
        else:
            value = 0.52 * compact + 0.48 * base
        composed[pair] = _rounded(value)
    return composed, compact_scores, dca_used or base_dca_used


def _compose_long_range_temperature_scores(base_scores: Mapping[ContactPair, float]) -> dict[ContactPair, float]:
    scores: dict[ContactPair, float] = {}
    for pair, value in base_scores.items():
        separation = pair[1] - pair[0]
        long_range_boost = 0.16 if separation >= 24 else 0.0
        local_penalty = 0.10 if separation < 8 else 0.0
        scores[pair] = _rounded(float(value) + long_range_boost - local_penalty)
    return scores


def _solve_channel(
    *,
    row: RealCoordinateVisualRow,
    scores: Mapping[ContactPair, float],
    channel: str,
    pair_pool_multiplier: float,
    max_selected_contact_multiplier: float,
    ensemble_size: int,
    max_residue_degree: int,
) -> tuple[tuple[ContactPair, ...], tuple[GlobalContactFactor, ...], tuple[GlobalFactorGraphSolution, ...]]:
    pair_pool = _top_pair_pool(scores, sequence_length=row.sequence_length, multiplier=pair_pool_multiplier)
    factors = build_contact_factors(
        row=row,
        source_mode=channel,
        scored_pairs=scores,
        pair_pool=pair_pool,
        factor_source=f"{channel}:dense_geometry_annealing_v0",
        neighbourhood_radius=2,
        min_relative_score=0.0,
        max_factor_count=max(64, int(round(row.sequence_length * 1.65))),
        learned_geometry_prior_used_before_selection=False,
    )
    solutions, solution_pairs = solve_top_k_factor_graph(
        row=row,
        factors=factors,
        source_mode=channel,
        ensemble_size=ensemble_size,
        max_selected_contacts=max(64, int(round(row.sequence_length * max_selected_contact_multiplier))),
        max_residue_degree=max_residue_degree,
        diversity_temperature=0.020,
    )
    consensus = _consensus_pairs(solution_pairs, threshold=0.50)
    return consensus, factors, solutions


def _matched_controls_for_report(
    *,
    row: RealCoordinateVisualRow,
    selected_pairs: Sequence[ContactPair],
    candidate_pairs: Sequence[ContactPair],
) -> tuple[float, float, int]:
    controls: list[ContactMetricPacket] = []
    count = max(2, min(6, int(sqrt(max(1, len(selected_pairs))))))
    for control_index in range(1, count + 1):
        control_pairs = matched_control_pairs(
            row=row,
            selected_pairs=selected_pairs,
            candidate_pairs=candidate_pairs,
            control_index=control_index,
        )
        if not control_pairs:
            continue
        controls.append(evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=control_pairs))
    return (
        max((metric.contact_map_f1 for metric in controls), default=0.0),
        max((metric.long_range_contact_recall for metric in controls), default=0.0),
        len(controls),
    )


def _row_claim_decision(
    *,
    source_mode: str,
    metric: ContactMetricPacket,
    f1_margin: float,
    long_range_margin: float,
    coordinate_taint: bool,
    native_taint: bool,
    structure_model: bool,
    anchor_count: int,
) -> tuple[bool, bool, str]:
    if coordinate_taint or native_taint:
        return False, False, "row_claim_rejected_coordinate_or_native_tainted_anchor"
    if structure_model:
        return False, False, "row_claim_rejected_structure_model_used_before_selection"
    if source_mode != PURE_SEQUENCE_COARSE_ANNEALING_MODE and anchor_count <= 0:
        return False, False, "row_claim_rejected_no_safe_external_dca_anchors"
    claim_allowed = (
        metric.native_contact_precision >= 0.70
        and metric.native_contact_recall >= 0.70
        and metric.long_range_contact_recall >= 0.70
        and f1_margin >= 0.15
        and long_range_margin >= 0.15
    )
    universal_allowed = claim_allowed and source_mode == PURE_SEQUENCE_COARSE_ANNEALING_MODE
    if claim_allowed and source_mode != PURE_SEQUENCE_COARSE_ANNEALING_MODE:
        return True, False, "dense_geometry_row_survived_external_dca_gate_not_universal_physics"
    if universal_allowed:
        return True, True, "pure_sequence_dense_geometry_row_survived_universal_physics_gate"
    if metric.native_contact_precision < 0.70:
        return False, False, "row_claim_rejected_precision_below_0_70"
    if metric.native_contact_recall < 0.70:
        return False, False, "row_claim_rejected_recall_below_0_70"
    if metric.long_range_contact_recall < 0.70:
        return False, False, "row_claim_rejected_long_range_recall_below_0_70"
    return False, False, "row_claim_rejected_matched_control_margin_below_gate"


def _decisions_for_row(
    *,
    row: RealCoordinateVisualRow,
    source_mode: str,
    base_scores: Mapping[ContactPair, float],
    compact_scores: Mapping[ContactPair, float],
    long_scores: Mapping[ContactPair, float],
    selected_pairs: Sequence[ContactPair],
    compact_pairs: Sequence[ContactPair],
    long_pairs: Sequence[ContactPair],
    dca_used: bool,
    max_report_decisions: int = 900,
) -> tuple[DenseGeometryContactDecision, ...]:
    selected = set(selected_pairs)
    compact = set(compact_pairs)
    long_range = set(long_pairs)
    candidate_pairs = selected | compact | long_range
    ordered = sorted(
        candidate_pairs,
        key=lambda pair: (
            pair not in selected,
            pair not in compact,
            pair not in long_range,
            -max(float(base_scores.get(pair, 0.0)), float(compact_scores.get(pair, 0.0)), float(long_scores.get(pair, 0.0))),
            pair[0],
            pair[1],
        ),
    )[:max_report_decisions]
    decisions: list[DenseGeometryContactDecision] = []
    for pair in ordered:
        compact_value = float(compact_scores.get(pair, 0.0))
        long_value = float(long_scores.get(pair, 0.0))
        base_value = float(base_scores.get(pair, 0.0))
        final_value = max(compact_value, long_value, base_value)
        if pair in compact and pair in long_range:
            channel = MULTIFIELD_CHANNEL
        elif pair in compact:
            channel = COMPACT_CHANNEL
        else:
            channel = LONG_RANGE_CHANNEL
        decisions.append(
            DenseGeometryContactDecision(
                kind=DENSE_GEOMETRY_ANNEALING_DECISION_KIND,
                row_id=row.row_id,
                source_accession=row.source_accession,
                source_mode=source_mode,
                selection_channel=channel,
                i=pair[0],
                j=pair[1],
                sequence_separation=pair[1] - pair[0],
                base_score=_rounded(base_value),
                compact_geometry_score=_rounded(compact_value),
                long_range_temperature_score=_rounded(long_value),
                final_score=_rounded(final_value),
                selected=pair in selected,
                selected_by_compact_channel=pair in compact,
                selected_by_long_range_channel=pair in long_range,
                msa_dca_used_before_selection=dca_used,
            )
        )
    return tuple(decisions)


def run_dense_geometry_annealing_row(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint] = (),
    source_mode: str = EXTERNAL_DCA_MULTIFIELD_ANNEALING_MODE,
    relaxation_steps: int = 32,
    pair_pool_multiplier: float = 4.5,
    max_selected_contact_multiplier: float = 2.0,
    ensemble_size: int = 3,
    max_residue_degree: int = 12,
) -> tuple[
    DenseGeometryRowReport,
    tuple[GlobalContactFactor, ...],
    tuple[GlobalFactorGraphSolution, ...],
    tuple[DenseGeometryContactDecision, ...],
]:
    safe_constraints, coordinate_taint, native_taint, structure_model, safe_dca_used = _safe_constraints_for_row(
        row=row,
        constraints=constraints,
        source_mode=source_mode,
    )
    base_scores, base_dca_used = _base_scores_for_mode(row=row, constraints=constraints, source_mode=source_mode)
    compact_scores, compact_geometry_scores, compact_dca_used = _compose_compact_scores(
        row=row,
        constraints=constraints,
        source_mode=source_mode,
        relaxation_steps=relaxation_steps,
    )
    long_scores = _compose_long_range_temperature_scores(base_scores)

    compact_pairs, compact_factors, compact_solutions = _solve_channel(
        row=row,
        scores=compact_scores,
        channel=COMPACT_CHANNEL,
        pair_pool_multiplier=pair_pool_multiplier,
        max_selected_contact_multiplier=max_selected_contact_multiplier,
        ensemble_size=ensemble_size,
        max_residue_degree=max_residue_degree,
    )
    if source_mode == PURE_SEQUENCE_COARSE_ANNEALING_MODE:
        long_pairs: tuple[ContactPair, ...] = ()
        long_factors: tuple[GlobalContactFactor, ...] = ()
        long_solutions: tuple[GlobalFactorGraphSolution, ...] = ()
        selected_pairs = compact_pairs
    else:
        long_pairs, long_factors, long_solutions = _solve_channel(
            row=row,
            scores=long_scores,
            channel=LONG_RANGE_CHANNEL,
            pair_pool_multiplier=pair_pool_multiplier,
            max_selected_contact_multiplier=max_selected_contact_multiplier,
            ensemble_size=ensemble_size,
            max_residue_degree=max_residue_degree,
        )
        if source_mode == EXTERNAL_DCA_COMPACT_ANNEALING_MODE:
            selected_pairs = compact_pairs
        else:
            selected_pairs = normalized_contact_pairs(set(compact_pairs) | set(long_pairs))

    compact_metric = evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=compact_pairs)
    long_metric = evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=long_pairs)
    metric = evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=selected_pairs)
    candidate_pairs = _top_pair_pool(base_scores, sequence_length=row.sequence_length, multiplier=pair_pool_multiplier)
    best_control_f1, best_control_lr, control_count = _matched_controls_for_report(
        row=row,
        selected_pairs=selected_pairs,
        candidate_pairs=candidate_pairs,
    )
    f1_margin = _score(metric.contact_map_f1 - best_control_f1)
    lr_margin = _score(metric.long_range_contact_recall - best_control_lr)
    dca_used = safe_dca_used or base_dca_used or compact_dca_used
    row_claim, row_universal, rejection = _row_claim_decision(
        source_mode=source_mode,
        metric=metric,
        f1_margin=f1_margin,
        long_range_margin=lr_margin,
        coordinate_taint=coordinate_taint,
        native_taint=native_taint,
        structure_model=structure_model,
        anchor_count=len(safe_constraints),
    )
    decisions = _decisions_for_row(
        row=row,
        source_mode=source_mode,
        base_scores=base_scores,
        compact_scores=compact_geometry_scores,
        long_scores=long_scores,
        selected_pairs=selected_pairs,
        compact_pairs=compact_pairs,
        long_pairs=long_pairs,
        dca_used=dca_used,
    )
    report = DenseGeometryRowReport(
        row_id=row.row_id,
        source_accession=row.source_accession,
        source_mode=source_mode,
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        safe_anchor_count=len(safe_constraints),
        relaxation_step_count=relaxation_steps,
        compact_channel_selected_count=len(compact_pairs),
        long_range_channel_selected_count=len(long_pairs),
        selected_contact_count=len(selected_pairs),
        selected_long_range_contact_count=sum(1 for pair in selected_pairs if pair[1] - pair[0] >= 24),
        selected_contact_map_hash=_pair_hash(selected_pairs),
        compact_channel_metric_after_native_audit=compact_metric,
        long_range_channel_metric_after_native_audit=long_metric,
        metric_after_native_audit=metric,
        matched_control_count=control_count,
        best_control_f1_after_audit=best_control_f1,
        best_control_long_range_recall_after_audit=best_control_lr,
        f1_margin_vs_best_control=f1_margin,
        long_range_recall_margin_vs_best_control=lr_margin,
        row_dense_geometry_claim_allowed=row_claim,
        row_universal_physical_law_claim_allowed=row_universal,
        row_claim_rejection_reason=rejection,
        coordinate_truth_used_before_selection=coordinate_taint,
        native_truth_used_before_selection=native_taint,
        structure_model_used_before_selection=structure_model,
        msa_dca_used_before_selection=dca_used,
    )
    return report, compact_factors + long_factors, compact_solutions + long_solutions, decisions


def run_dense_geometry_annealing_packet(
    rows: Sequence[RealCoordinateVisualRow],
    *,
    constraints: Sequence[CouplingConstraint] = (),
    source_mode: str = EXTERNAL_DCA_MULTIFIELD_ANNEALING_MODE,
    relaxation_steps: int = 32,
    pair_pool_multiplier: float = 4.5,
    max_selected_contact_multiplier: float = 2.0,
    ensemble_size: int = 3,
) -> DenseGeometryAnnealingPacket:
    row_reports: list[DenseGeometryRowReport] = []
    all_factors: list[GlobalContactFactor] = []
    all_solutions: list[GlobalFactorGraphSolution] = []
    all_decisions: list[DenseGeometryContactDecision] = []
    for row in rows:
        report, factors, solutions, decisions = run_dense_geometry_annealing_row(
            row=row,
            constraints=constraints,
            source_mode=source_mode,
            relaxation_steps=relaxation_steps,
            pair_pool_multiplier=pair_pool_multiplier,
            max_selected_contact_multiplier=max_selected_contact_multiplier,
            ensemble_size=ensemble_size,
        )
        row_reports.append(report)
        all_factors.extend(factors[:260])
        all_solutions.extend(solutions)
        selected = [decision for decision in decisions if decision.selected]
        rejected = [decision for decision in decisions if not decision.selected][:80]
        all_decisions.extend(selected + rejected)

    precision_values = [row.metric_after_native_audit.native_contact_precision for row in row_reports]
    recall_values = [row.metric_after_native_audit.native_contact_recall for row in row_reports]
    long_range_values = [row.metric_after_native_audit.long_range_contact_recall for row in row_reports]
    f1_values = [row.metric_after_native_audit.contact_map_f1 for row in row_reports]
    compact_f1_values = [row.compact_channel_metric_after_native_audit.contact_map_f1 for row in row_reports]
    compact_lr_values = [row.compact_channel_metric_after_native_audit.long_range_contact_recall for row in row_reports]
    long_f1_values = [row.long_range_channel_metric_after_native_audit.contact_map_f1 for row in row_reports]
    long_lr_values = [row.long_range_channel_metric_after_native_audit.long_range_contact_recall for row in row_reports]
    f1_margins = [row.f1_margin_vs_best_control for row in row_reports]
    lr_margins = [row.long_range_recall_margin_vs_best_control for row in row_reports]

    coordinate_taint = any(row.coordinate_truth_used_before_selection for row in row_reports)
    native_taint = any(row.native_truth_used_before_selection for row in row_reports)
    structure_model = any(row.structure_model_used_before_selection for row in row_reports)
    dca_used = any(row.msa_dca_used_before_selection for row in row_reports)
    all_rows_claim = bool(row_reports) and all(row.row_dense_geometry_claim_allowed for row in row_reports)
    all_rows_universal = bool(row_reports) and all(row.row_universal_physical_law_claim_allowed for row in row_reports)
    mean_gate = _mean(precision_values) >= 0.70 and _mean(recall_values) >= 0.70 and _mean(long_range_values) >= 0.70
    dense_claim = all_rows_claim and mean_gate and not coordinate_taint and not native_taint and not structure_model
    universal_claim = all_rows_universal and dense_claim and source_mode == PURE_SEQUENCE_COARSE_ANNEALING_MODE and not dca_used

    if dense_claim and source_mode != PURE_SEQUENCE_COARSE_ANNEALING_MODE:
        rejection = "external_dca_dense_geometry_annealing_survived_gate_not_universal_physics"
    elif universal_claim:
        rejection = "pure_sequence_dense_geometry_annealing_survived_universal_physical_law_gate"
    else:
        failed = [row.source_accession for row in row_reports if not row.row_dense_geometry_claim_allowed]
        rejection = "dense_geometry_annealing_claim_rejected_for_rows:" + ",".join(failed[:12])

    return DenseGeometryAnnealingPacket(
        kind=DENSE_GEOMETRY_ANNEALING_KIND,
        source_mode=source_mode,
        row_count=len(row_reports),
        decision_rule=DENSE_GEOMETRY_ANNEALING_RULE,
        claim_rule=DENSE_GEOMETRY_ANNEALING_CLAIM_RULE,
        pure_sequence_mode_included=source_mode == PURE_SEQUENCE_COARSE_ANNEALING_MODE,
        external_dca_anchor_field_included=source_mode != PURE_SEQUENCE_COARSE_ANNEALING_MODE,
        coarse_annealing_included=True,
        long_range_temperature_channel_included=source_mode == EXTERNAL_DCA_MULTIFIELD_ANNEALING_MODE,
        global_factor_graph_included=True,
        mean_native_contact_precision_after_audit=_mean(precision_values),
        mean_native_contact_recall_after_audit=_mean(recall_values),
        mean_long_range_contact_recall_after_audit=_mean(long_range_values),
        mean_contact_map_f1_after_audit=_mean(f1_values),
        mean_compact_channel_f1_after_audit=_mean(compact_f1_values),
        mean_compact_channel_long_range_recall_after_audit=_mean(compact_lr_values),
        mean_long_range_channel_f1_after_audit=_mean(long_f1_values),
        mean_long_range_channel_long_range_recall_after_audit=_mean(long_lr_values),
        min_row_precision_after_audit=_rounded(min(precision_values)) if precision_values else 0.0,
        min_row_recall_after_audit=_rounded(min(recall_values)) if recall_values else 0.0,
        min_row_long_range_recall_after_audit=_rounded(min(long_range_values)) if long_range_values else 0.0,
        mean_f1_margin_vs_best_control=_mean(f1_margins),
        mean_long_range_recall_margin_vs_best_control=_mean(lr_margins),
        dense_geometry_annealing_claim_allowed=dense_claim,
        universal_physical_law_claim_allowed=universal_claim,
        folding_problem_solved=dense_claim or universal_claim,
        claim_rejection_reason=rejection,
        rows=tuple(row_reports),
        factors=tuple(all_factors),
        solutions=tuple(all_solutions),
        decisions=tuple(all_decisions),
        coordinate_truth_used_before_selection=coordinate_taint,
        native_truth_used_before_selection=native_taint,
        structure_model_used_before_selection=structure_model,
        learned_geometry_prior_used_before_selection=False,
        msa_dca_used_before_selection=dca_used,
        raw_sequence_exposed=False,
    )
