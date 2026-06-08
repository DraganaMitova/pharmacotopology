from __future__ import annotations

"""Native-free global geometry-field generator plus factor graph.

This layer tests the next honest question after the ESMFold/factor-graph result:
can a usable global geometry prior be generated without AlphaFold/ESMFold,
templates, or any trained structure model?

The implementation deliberately separates three modes:

* pure_sequence_symbolic_geometry: sequence chemistry + symbolic topology only.
* external_dca_geometry_field: safe external evolutionary couplings + sequence
  symbolic topology. This is native-free and template-free, but it is not a
  universal physical-law claim because the coupling channel comes from an MSA.
* oracle_anchor_geometry_field_control: coordinate-native oracle anchors. This
  mode is an upper-bound/debug control and can never open a claim gate.

Native/native-coordinate truth is attached only after selection for audit. The
factor graph itself receives only generated pair scores and all-or-none contact
patch factors.
"""

import hashlib
from collections import defaultdict
from dataclasses import asdict, dataclass
from math import exp, sqrt
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
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    MIN_SEQUENCE_SEPARATION,
    RealCoordinateVisualRow,
)
from pharmacotopology.folding_sequence_physical_priors import (
    NEGATIVE_CHARGED,
    POSITIVE_CHARGED,
    build_sequence_physical_prior_scores,
    predict_lightweight_secondary_structure,
)
from pharmacotopology.folding_topology import HYDROPHOBIC_AMINO_ACIDS


NATIVE_FREE_GEOMETRY_FIELD_KIND = "native_free_global_geometry_field_factor_graph_v0"
NATIVE_FREE_GEOMETRY_DECISION_KIND = "native_free_geometry_field_contact_decision_v0"
PURE_SEQUENCE_SYMBOLIC_GEOMETRY_MODE = "pure_sequence_symbolic_geometry"
EXTERNAL_DCA_GEOMETRY_FIELD_MODE = "external_dca_geometry_field"
ORACLE_ANCHOR_GEOMETRY_FIELD_CONTROL_MODE = "oracle_anchor_geometry_field_control"
GEOMETRY_FIELD_FACTOR_GRAPH_MODE = "native_free_geometry_field_factor_graph"
GEOMETRY_FIELD_RULE = (
    "sequence_symbolic_topology+optional_safe_dca_anchor_field;"
    "contact_map_grid_diffusion_without_mds;global_all_or_none_factor_graph;"
    "top_k_consensus;native_audit_after_selection_only"
)
GEOMETRY_FIELD_CLAIM_RULE = (
    "claim_requires_no_coordinate_truth_no_native_truth_no_structure_model_no_learned_geometry;"
    "all_rows_precision_ge_0_70_and_recall_ge_0_70_and_long_range_recall_ge_0_70;"
    "matched_control_f1_and_long_range_margin_ge_0_15;"
    "pure_sequence_only_required_for_universal_physical_law_claim"
)


@dataclass(frozen=True)
class GeometryFieldPairScore:
    row_id: str
    source_accession: str
    i: int
    j: int
    sequence_separation: int
    final_score: float
    physical_prior_score: float
    symbolic_topology_score: float
    separation_geometry_score: float
    anchor_field_score: float
    anchor_support_count: int
    selected: bool
    selection_probability: float
    factor_vote_count: int
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
class GeometryFieldRowReport:
    row_id: str
    source_accession: str
    source_mode: str
    sequence_hash: str
    sequence_length: int
    candidate_pair_count: int
    safe_anchor_count: int
    selected_pair_pool_count: int
    all_or_none_factor_count: int
    ensemble_solution_count: int
    consensus_probability_threshold: float
    selected_contact_count: int
    selected_long_range_contact_count: int
    selected_contact_map_hash: str
    metric_after_native_audit: ContactMetricPacket
    best_solution_metric_after_native_audit: ContactMetricPacket
    matched_control_count: int
    best_control_f1_after_audit: float
    best_control_long_range_recall_after_audit: float
    f1_margin_vs_best_control: float
    long_range_recall_margin_vs_best_control: float
    row_geometry_field_claim_allowed: bool
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
        payload["metric_after_native_audit"] = self.metric_after_native_audit.to_dict()
        payload["best_solution_metric_after_native_audit"] = self.best_solution_metric_after_native_audit.to_dict()
        return payload


@dataclass(frozen=True)
class NativeFreeGeometryFieldPacket:
    kind: str
    source_mode: str
    row_count: int
    decision_rule: str
    claim_rule: str
    pure_sequence_symbolic_geometry_included: bool
    external_dca_anchor_field_included: bool
    global_factor_graph_included: bool
    top_k_ensemble_included: bool
    mean_native_contact_precision_after_audit: float
    mean_native_contact_recall_after_audit: float
    mean_long_range_contact_recall_after_audit: float
    mean_contact_map_f1_after_audit: float
    min_row_precision_after_audit: float
    min_row_recall_after_audit: float
    min_row_long_range_recall_after_audit: float
    mean_f1_margin_vs_best_control: float
    mean_long_range_recall_margin_vs_best_control: float
    native_free_geometry_field_claim_allowed: bool
    universal_physical_law_claim_allowed: bool
    folding_problem_solved: bool
    claim_rejection_reason: str
    rows: tuple[GeometryFieldRowReport, ...]
    factors: tuple[GlobalContactFactor, ...]
    solutions: tuple[GlobalFactorGraphSolution, ...]
    decisions: tuple[GeometryFieldPairScore, ...]
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
            "pure_sequence_symbolic_geometry_included": self.pure_sequence_symbolic_geometry_included,
            "external_dca_anchor_field_included": self.external_dca_anchor_field_included,
            "global_factor_graph_included": self.global_factor_graph_included,
            "top_k_ensemble_included": self.top_k_ensemble_included,
            "mean_native_contact_precision_after_audit": self.mean_native_contact_precision_after_audit,
            "mean_native_contact_recall_after_audit": self.mean_native_contact_recall_after_audit,
            "mean_long_range_contact_recall_after_audit": self.mean_long_range_contact_recall_after_audit,
            "mean_contact_map_f1_after_audit": self.mean_contact_map_f1_after_audit,
            "min_row_precision_after_audit": self.min_row_precision_after_audit,
            "min_row_recall_after_audit": self.min_row_recall_after_audit,
            "min_row_long_range_recall_after_audit": self.min_row_long_range_recall_after_audit,
            "mean_f1_margin_vs_best_control": self.mean_f1_margin_vs_best_control,
            "mean_long_range_recall_margin_vs_best_control": self.mean_long_range_recall_margin_vs_best_control,
            "native_free_geometry_field_claim_allowed": self.native_free_geometry_field_claim_allowed,
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


def _candidate_pairs(row: RealCoordinateVisualRow) -> tuple[ContactPair, ...]:
    return tuple(
        (left, right)
        for left in range(1, row.sequence_length + 1)
        for right in range(left + MIN_SEQUENCE_SEPARATION, row.sequence_length + 1)
    )


def _safe_constraints_for_mode(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    source_mode: str,
) -> tuple[tuple[CouplingConstraint, ...], bool, bool, bool, bool]:
    if source_mode == PURE_SEQUENCE_SYMBOLIC_GEOMETRY_MODE:
        return (), False, False, False, False

    row_constraints = tuple(
        constraint
        for constraint in constraints
        if constraint.source_accession == row.source_accession or constraint.row_id == row.row_id
    )
    if source_mode == ORACLE_ANCHOR_GEOMETRY_FIELD_CONTROL_MODE:
        coordinate_taint = any(constraint.coordinate_truth_used_to_build_constraint for constraint in row_constraints)
        native_taint = any(constraint.native_truth_used_before_coupling_selection for constraint in row_constraints)
        # The legacy oracle-control file stores native taint primarily at the
        # dataset boundary and coordinate taint on each row. In this control mode
        # coordinate-derived anchors are treated as native/coordinate tainted as
        # a whole, so the claim gate can never accidentally open.
        native_taint = native_taint or coordinate_taint
        structure_model = any(constraint.structure_model_used for constraint in row_constraints)
        return row_constraints, coordinate_taint, native_taint, structure_model, False

    safe: list[CouplingConstraint] = []
    coordinate_taint = False
    native_taint = False
    structure_model = False
    for constraint in row_constraints:
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
    return tuple(safe), coordinate_taint, native_taint, structure_model, bool(safe)


def _constraint_sort_key(constraint: CouplingConstraint) -> tuple[float, int, int, int]:
    rank = constraint.rank if constraint.rank else 999999
    return (-float(constraint.confidence), rank, constraint.i, constraint.j)


def _anchor_field_scores(
    candidate_pairs: Sequence[ContactPair],
    anchors: Sequence[tuple[ContactPair, float]],
    *,
    max_anchors: int = 96,
    radius: float = 8.0,
) -> tuple[dict[ContactPair, float], dict[ContactPair, int]]:
    kept = tuple(anchors[:max_anchors])
    if not kept:
        return ({pair: 0.0 for pair in candidate_pairs}, {pair: 0 for pair in candidate_pairs})
    field: dict[ContactPair, float] = {}
    support_counts: dict[ContactPair, int] = {}
    for pair in candidate_pairs:
        i, j = pair
        best = 0.0
        support = 0
        for (ai, aj), confidence in kept:
            direct = abs(i - ai) + abs(j - aj)
            registry = abs((i - j) - (ai - aj))
            value = float(confidence) * max(
                exp(-direct / radius),
                0.48 * exp(-registry / (radius * 1.4)),
            )
            if value >= 0.34:
                support += 1
            if value > best:
                best = value
        field[pair] = _rounded(best)
        support_counts[pair] = support
    return field, support_counts


def _opposite_charge(left: str, right: str) -> bool:
    return (left in POSITIVE_CHARGED and right in NEGATIVE_CHARGED) or (
        right in POSITIVE_CHARGED and left in NEGATIVE_CHARGED
    )


def _symbolic_topology_score(row: RealCoordinateVisualRow, pair: ContactPair, secondary_structure: Sequence[str]) -> float:
    i, j = pair
    separation = j - i
    left = row.sequence[i - 1]
    right = row.sequence[j - 1]
    left_ss = secondary_structure[i - 1]
    right_ss = secondary_structure[j - 1]
    score = 0.0
    if left_ss == "H" and right_ss == "H" and separation in (3, 4, 5):
        score = max(score, 1.0)
    elif left_ss == "E" and right_ss == "E" and separation >= 8:
        score = max(score, 0.60)
    if left in HYDROPHOBIC_AMINO_ACIDS and right in HYDROPHOBIC_AMINO_ACIDS and separation >= 8:
        score = max(score, 0.64)
    if _opposite_charge(left, right) and separation >= 8:
        score = max(score, 0.58)
    return _rounded(score)


def _separation_geometry_score(row: RealCoordinateVisualRow, pair: ContactPair) -> float:
    # Native-free broad fold-scale prior: contacts often populate loop closures
    # around a middle normalized separation, but this is deliberately weak.
    normalized = (pair[1] - pair[0]) / max(1, row.sequence_length)
    return _rounded(1.0 / (1.0 + abs(normalized - 0.28) * 5.2))


def build_native_free_geometry_scores(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint] = (),
    source_mode: str = EXTERNAL_DCA_GEOMETRY_FIELD_MODE,
) -> tuple[
    dict[ContactPair, float],
    dict[ContactPair, dict[str, float | int]],
    int,
    bool,
    bool,
    bool,
    bool,
]:
    candidate_pairs = _candidate_pairs(row)
    safe_constraints, coordinate_taint, native_taint, structure_model, dca_used = _safe_constraints_for_mode(
        row=row,
        constraints=constraints,
        source_mode=source_mode,
    )
    anchors = tuple(
        (constraint.pair(), float(constraint.confidence))
        for constraint in sorted(safe_constraints, key=_constraint_sort_key)
        if constraint.sequence_separation >= MIN_SEQUENCE_SEPARATION
    )
    anchor_scores, anchor_counts = _anchor_field_scores(candidate_pairs, anchors)
    priors = build_sequence_physical_prior_scores(row=row, candidate_pairs=candidate_pairs, current_pairs=())
    secondary_structure = predict_lightweight_secondary_structure(row.sequence)

    final_scores: dict[ContactPair, float] = {}
    components: dict[ContactPair, dict[str, float | int]] = {}
    for pair in candidate_pairs:
        physical = float(priors[pair].physical_prior_score)
        symbolic = _symbolic_topology_score(row, pair, secondary_structure)
        separation = _separation_geometry_score(row, pair)
        anchor = float(anchor_scores[pair])
        if anchors:
            final = 0.47 * anchor + 0.25 * physical + 0.18 * symbolic + 0.10 * separation
        else:
            final = 0.55 * physical + 0.25 * symbolic + 0.20 * separation
        final_scores[pair] = _score(final)
        components[pair] = {
            "physical_prior_score": _rounded(physical),
            "symbolic_topology_score": _rounded(symbolic),
            "separation_geometry_score": _rounded(separation),
            "anchor_field_score": _rounded(anchor),
            "anchor_support_count": int(anchor_counts[pair]),
        }
    return final_scores, components, len(anchors), coordinate_taint, native_taint, structure_model, dca_used


def _top_pair_pool(
    scores: Mapping[ContactPair, float],
    *,
    sequence_length: int,
    multiplier: float,
) -> tuple[ContactPair, ...]:
    limit = max(80, int(round(sequence_length * multiplier)))
    return normalized_contact_pairs(
        pair
        for pair, _value in sorted(scores.items(), key=lambda item: (-item[1], item[0][0], item[0][1]))[:limit]
    )


def _consensus_pairs(solution_pairs: Sequence[Sequence[ContactPair]], *, threshold: float) -> tuple[ContactPair, ...]:
    counts: dict[ContactPair, int] = defaultdict(int)
    for pairs in solution_pairs:
        for pair in set(pairs):
            counts[pair] += 1
    denom = max(1, len(solution_pairs))
    return normalized_contact_pairs(pair for pair, count in counts.items() if count / denom >= threshold)


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


def _contact_decisions(
    *,
    row: RealCoordinateVisualRow,
    source_mode: str,
    scores: Mapping[ContactPair, float],
    components: Mapping[ContactPair, Mapping[str, float | int]],
    factors: Sequence[GlobalContactFactor],
    solution_pairs: Sequence[Sequence[ContactPair]],
    consensus_pairs: Sequence[ContactPair],
    coordinate_taint: bool,
    native_taint: bool,
    structure_model: bool,
    dca_used: bool,
    max_report_decisions: int = 1200,
) -> tuple[GeometryFieldPairScore, ...]:
    factor_votes: dict[ContactPair, int] = defaultdict(int)
    for factor in factors:
        for pair in factor.contacts:
            factor_votes[pair] += 1
    solution_counts: dict[ContactPair, int] = defaultdict(int)
    for pairs in solution_pairs:
        for pair in set(pairs):
            solution_counts[pair] += 1
    denom = max(1, len(solution_pairs))
    selected = set(consensus_pairs)
    candidate_report_pairs = set(selected) | set(factor_votes)
    ordered = sorted(
        candidate_report_pairs,
        key=lambda pair: (
            pair not in selected,
            -solution_counts.get(pair, 0) / denom,
            -float(scores.get(pair, 0.0)),
            pair[0],
            pair[1],
        ),
    )[:max_report_decisions]
    decisions: list[GeometryFieldPairScore] = []
    for pair in ordered:
        row_components = components.get(pair, {})
        decisions.append(
            GeometryFieldPairScore(
                row_id=row.row_id,
                source_accession=row.source_accession,
                i=pair[0],
                j=pair[1],
                sequence_separation=pair[1] - pair[0],
                final_score=_score(float(scores.get(pair, 0.0))),
                physical_prior_score=_rounded(float(row_components.get("physical_prior_score", 0.0))),
                symbolic_topology_score=_rounded(float(row_components.get("symbolic_topology_score", 0.0))),
                separation_geometry_score=_rounded(float(row_components.get("separation_geometry_score", 0.0))),
                anchor_field_score=_rounded(float(row_components.get("anchor_field_score", 0.0))),
                anchor_support_count=int(row_components.get("anchor_support_count", 0)),
                selected=pair in selected,
                selection_probability=_rounded(solution_counts.get(pair, 0) / denom),
                factor_vote_count=factor_votes.get(pair, 0),
                coordinate_truth_used_before_selection=coordinate_taint,
                native_truth_used_before_selection=native_taint,
                structure_model_used_before_selection=structure_model,
                msa_dca_used_before_selection=dca_used,
            )
        )
    return tuple(decisions)


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
    if source_mode == ORACLE_ANCHOR_GEOMETRY_FIELD_CONTROL_MODE or coordinate_taint or native_taint:
        return False, False, "oracle_or_coordinate_native_tainted_anchor_control_not_claimable"
    if structure_model:
        return False, False, "row_claim_rejected_structure_model_used_before_selection"
    if source_mode == EXTERNAL_DCA_GEOMETRY_FIELD_MODE and anchor_count <= 0:
        return False, False, "row_claim_rejected_no_safe_external_dca_anchors"
    claim_allowed = (
        metric.native_contact_precision >= 0.70
        and metric.native_contact_recall >= 0.70
        and metric.long_range_contact_recall >= 0.70
        and f1_margin >= 0.15
        and long_range_margin >= 0.15
    )
    universal_allowed = claim_allowed and source_mode == PURE_SEQUENCE_SYMBOLIC_GEOMETRY_MODE
    if claim_allowed and source_mode == EXTERNAL_DCA_GEOMETRY_FIELD_MODE:
        return True, False, "native_free_geometry_field_row_survived_external_dca_gate_not_universal_physics"
    if universal_allowed:
        return True, True, "pure_sequence_geometry_field_row_survived_universal_physics_gate"
    if metric.native_contact_precision < 0.70:
        return False, False, "row_claim_rejected_precision_below_0_70"
    if metric.native_contact_recall < 0.70:
        return False, False, "row_claim_rejected_recall_below_0_70"
    if metric.long_range_contact_recall < 0.70:
        return False, False, "row_claim_rejected_long_range_recall_below_0_70"
    return False, False, "row_claim_rejected_matched_control_margin_below_gate"


def run_native_free_geometry_field_row(
    *,
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint] = (),
    source_mode: str = EXTERNAL_DCA_GEOMETRY_FIELD_MODE,
    pair_pool_multiplier: float = 6.0,
    max_selected_contact_multiplier: float = 2.5,
    ensemble_size: int = 6,
    consensus_probability_threshold: float = 0.50,
    max_residue_degree: int = 12,
) -> tuple[
    GeometryFieldRowReport,
    tuple[GlobalContactFactor, ...],
    tuple[GlobalFactorGraphSolution, ...],
    tuple[GeometryFieldPairScore, ...],
]:
    scores, components, anchor_count, coordinate_taint, native_taint, structure_model, dca_used = build_native_free_geometry_scores(
        row=row,
        constraints=constraints,
        source_mode=source_mode,
    )
    pair_pool = _top_pair_pool(scores, sequence_length=row.sequence_length, multiplier=pair_pool_multiplier)
    factors = build_contact_factors(
        row=row,
        source_mode=GEOMETRY_FIELD_FACTOR_GRAPH_MODE,
        scored_pairs=scores,
        pair_pool=pair_pool,
        factor_source=f"{source_mode}_generated_contact_map_geometry_field",
        neighbourhood_radius=2,
        min_relative_score=0.0,
        max_factor_count=max(80, int(round(row.sequence_length * 2.5))),
        learned_geometry_prior_used_before_selection=False,
    )
    solutions, solution_pairs = solve_top_k_factor_graph(
        row=row,
        factors=factors,
        source_mode=GEOMETRY_FIELD_FACTOR_GRAPH_MODE,
        ensemble_size=ensemble_size,
        max_selected_contacts=max(64, int(round(row.sequence_length * max_selected_contact_multiplier))),
        max_residue_degree=max_residue_degree,
        diversity_temperature=0.020,
    )
    consensus = _consensus_pairs(solution_pairs, threshold=consensus_probability_threshold)
    metric = evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=consensus)
    best_solution_metric = max(
        (solution.metric_after_native_audit for solution in solutions),
        key=lambda item: (item.contact_map_f1, item.native_contact_recall, item.native_contact_precision),
        default=evaluate_contact_prediction(native_pairs=row.native_contact_pairs(), predicted_pairs=()),
    )
    best_control_f1, best_control_lr, control_count = _matched_controls_for_report(
        row=row,
        selected_pairs=consensus,
        candidate_pairs=pair_pool,
    )
    f1_margin = _score(metric.contact_map_f1 - best_control_f1)
    lr_margin = _score(metric.long_range_contact_recall - best_control_lr)
    row_claim, row_universal, rejection = _row_claim_decision(
        source_mode=source_mode,
        metric=metric,
        f1_margin=f1_margin,
        long_range_margin=lr_margin,
        coordinate_taint=coordinate_taint,
        native_taint=native_taint,
        structure_model=structure_model,
        anchor_count=anchor_count,
    )
    decisions = _contact_decisions(
        row=row,
        source_mode=source_mode,
        scores=scores,
        components=components,
        factors=factors,
        solution_pairs=solution_pairs,
        consensus_pairs=consensus,
        coordinate_taint=coordinate_taint,
        native_taint=native_taint,
        structure_model=structure_model,
        dca_used=dca_used,
    )
    report = GeometryFieldRowReport(
        row_id=row.row_id,
        source_accession=row.source_accession,
        source_mode=source_mode,
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        candidate_pair_count=len(scores),
        safe_anchor_count=anchor_count,
        selected_pair_pool_count=len(pair_pool),
        all_or_none_factor_count=len(factors),
        ensemble_solution_count=len(solutions),
        consensus_probability_threshold=consensus_probability_threshold,
        selected_contact_count=len(consensus),
        selected_long_range_contact_count=sum(1 for pair in consensus if pair[1] - pair[0] >= 24),
        selected_contact_map_hash=_pair_hash(consensus),
        metric_after_native_audit=metric,
        best_solution_metric_after_native_audit=best_solution_metric,
        matched_control_count=control_count,
        best_control_f1_after_audit=best_control_f1,
        best_control_long_range_recall_after_audit=best_control_lr,
        f1_margin_vs_best_control=f1_margin,
        long_range_recall_margin_vs_best_control=lr_margin,
        row_geometry_field_claim_allowed=row_claim,
        row_universal_physical_law_claim_allowed=row_universal,
        row_claim_rejection_reason=rejection,
        coordinate_truth_used_before_selection=coordinate_taint,
        native_truth_used_before_selection=native_taint,
        structure_model_used_before_selection=structure_model,
        msa_dca_used_before_selection=dca_used,
    )
    return report, factors, solutions, decisions


def run_native_free_geometry_field_packet(
    rows: Sequence[RealCoordinateVisualRow],
    *,
    constraints: Sequence[CouplingConstraint] = (),
    source_mode: str = EXTERNAL_DCA_GEOMETRY_FIELD_MODE,
    pair_pool_multiplier: float = 6.0,
    max_selected_contact_multiplier: float = 2.5,
    ensemble_size: int = 6,
) -> NativeFreeGeometryFieldPacket:
    row_reports: list[GeometryFieldRowReport] = []
    all_factors: list[GlobalContactFactor] = []
    all_solutions: list[GlobalFactorGraphSolution] = []
    all_decisions: list[GeometryFieldPairScore] = []
    for row in rows:
        report, factors, solutions, decisions = run_native_free_geometry_field_row(
            row=row,
            constraints=constraints,
            source_mode=source_mode,
            pair_pool_multiplier=pair_pool_multiplier,
            max_selected_contact_multiplier=max_selected_contact_multiplier,
            ensemble_size=ensemble_size,
        )
        row_reports.append(report)
        all_factors.extend(factors[:240])
        all_solutions.extend(solutions)
        selected = [decision for decision in decisions if decision.selected]
        rejected = [decision for decision in decisions if not decision.selected][:80]
        all_decisions.extend(selected + rejected)

    precision_values = [row.metric_after_native_audit.native_contact_precision for row in row_reports]
    recall_values = [row.metric_after_native_audit.native_contact_recall for row in row_reports]
    long_range_values = [row.metric_after_native_audit.long_range_contact_recall for row in row_reports]
    f1_values = [row.metric_after_native_audit.contact_map_f1 for row in row_reports]
    f1_margins = [row.f1_margin_vs_best_control for row in row_reports]
    lr_margins = [row.long_range_recall_margin_vs_best_control for row in row_reports]

    coordinate_taint = any(row.coordinate_truth_used_before_selection for row in row_reports)
    native_taint = any(row.native_truth_used_before_selection for row in row_reports)
    structure_model = any(row.structure_model_used_before_selection for row in row_reports)
    dca_used = any(row.msa_dca_used_before_selection for row in row_reports)
    all_rows_claim = bool(row_reports) and all(row.row_geometry_field_claim_allowed for row in row_reports)
    all_rows_universal = bool(row_reports) and all(row.row_universal_physical_law_claim_allowed for row in row_reports)
    mean_gate = _mean(precision_values) >= 0.70 and _mean(recall_values) >= 0.70 and _mean(long_range_values) >= 0.70
    geometry_claim = all_rows_claim and mean_gate and not coordinate_taint and not native_taint and not structure_model
    universal_claim = all_rows_universal and geometry_claim and source_mode == PURE_SEQUENCE_SYMBOLIC_GEOMETRY_MODE and not dca_used

    if source_mode == ORACLE_ANCHOR_GEOMETRY_FIELD_CONTROL_MODE or coordinate_taint or native_taint:
        rejection = "oracle_anchor_geometry_field_control_is_coordinate_native_tainted_not_claimable"
        geometry_claim = False
        universal_claim = False
    elif geometry_claim and source_mode == EXTERNAL_DCA_GEOMETRY_FIELD_MODE:
        rejection = "external_dca_geometry_field_claim_survived_native_free_gate_not_universal_physics"
    elif universal_claim:
        rejection = "pure_sequence_symbolic_geometry_field_survived_universal_physical_law_gate"
    else:
        failed = [row.source_accession for row in row_reports if not row.row_geometry_field_claim_allowed]
        rejection = "native_free_geometry_field_claim_rejected_for_rows:" + ",".join(failed[:12])

    return NativeFreeGeometryFieldPacket(
        kind=NATIVE_FREE_GEOMETRY_FIELD_KIND,
        source_mode=source_mode,
        row_count=len(row_reports),
        decision_rule=GEOMETRY_FIELD_RULE,
        claim_rule=GEOMETRY_FIELD_CLAIM_RULE,
        pure_sequence_symbolic_geometry_included=source_mode == PURE_SEQUENCE_SYMBOLIC_GEOMETRY_MODE,
        external_dca_anchor_field_included=source_mode == EXTERNAL_DCA_GEOMETRY_FIELD_MODE,
        global_factor_graph_included=True,
        top_k_ensemble_included=True,
        mean_native_contact_precision_after_audit=_mean(precision_values),
        mean_native_contact_recall_after_audit=_mean(recall_values),
        mean_long_range_contact_recall_after_audit=_mean(long_range_values),
        mean_contact_map_f1_after_audit=_mean(f1_values),
        min_row_precision_after_audit=_rounded(min(precision_values)) if precision_values else 0.0,
        min_row_recall_after_audit=_rounded(min(recall_values)) if recall_values else 0.0,
        min_row_long_range_recall_after_audit=_rounded(min(long_range_values)) if long_range_values else 0.0,
        mean_f1_margin_vs_best_control=_mean(f1_margins),
        mean_long_range_recall_margin_vs_best_control=_mean(lr_margins),
        native_free_geometry_field_claim_allowed=geometry_claim,
        universal_physical_law_claim_allowed=universal_claim,
        folding_problem_solved=geometry_claim or universal_claim,
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
