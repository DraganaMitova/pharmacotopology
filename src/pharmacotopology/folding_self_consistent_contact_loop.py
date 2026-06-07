from __future__ import annotations

import hashlib
from collections import defaultdict
from dataclasses import asdict, dataclass
from math import sqrt
from statistics import mean
from typing import Mapping, Sequence

from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_independent_contact_evidence import (
    IndependentContactEvidencePair,
)
from pharmacotopology.folding_native_contact_eval import (
    ContactMetricPacket,
    ContactPair,
    contact_map_hash,
    evaluate_contact_prediction,
    normalized_contact_pairs,
)
from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    MIN_SEQUENCE_SEPARATION,
    RealCoordinateVisualRow,
)
from pharmacotopology.folding_sequence_physical_priors import (
    SEQUENCE_PHYSICAL_PRIOR_KIND,
    SequenceContactPhysicalPrior,
    build_sequence_physical_prior_scores,
)


SELF_CONSISTENT_CONTACT_LOOP_KIND = "self_consistent_contact_loop_v0"
SELF_CONSISTENT_INTERNAL_CONTACT_SOURCE_KIND = "self_consistent_internal_contact_source_v0"
SELF_CONSISTENT_CONTROL_KIND = "self_consistent_matched_negative_control_v0"
SELF_CONSISTENT_SCORING_RULE = (
    "weighted_coupling_event_graph_plus_sequence_energy_secondary_structure_degree;"
    "row_local_largest_gap;adaptive_seed_derived_minimum;no_static_confidence_threshold"
)


@dataclass(frozen=True)
class GapBoundary:
    selected_count: int
    boundary_score: float
    reason: str
    largest_gap: float
    positive_score_count: int
    seed_identity_envelope_count: int

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class SelfConsistentIteration:
    loop_name: str
    iteration_index: int
    candidate_pair_count: int
    seed_pair_count: int
    selected_pair_count: int
    previous_selected_pair_count: int
    changed_pair_count: int
    previous_jaccard: float
    selected_contact_map_hash: str
    boundary: GapBoundary
    mean_selected_score: float
    mean_selected_contact_energy_score: float = 0.0
    mean_selected_secondary_structure_score: float = 0.0
    mean_selected_degree_consistency_score: float = 0.0
    mean_selected_physical_prior_score: float = 0.0
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["boundary"] = self.boundary.to_dict()
        return data


@dataclass(frozen=True)
class SelfConsistentControlResult:
    control_id: str
    seed_pair_count: int
    final_pair_count: int
    final_contact_map_hash: str
    self_consistency_strength: float
    final_stability: float
    seed_retention: float
    coupling_retention: float
    graph_closure_coherence: float
    boundary_gap: float
    mean_final_contact_energy_score: float = 0.0
    mean_final_secondary_structure_score: float = 0.0
    mean_final_degree_consistency_score: float = 0.0
    mean_final_physical_prior_score: float = 0.0
    control_kind: str = SELF_CONSISTENT_CONTROL_KIND
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class SelfConsistentContactLoopReport:
    kind: str
    row_id: str
    source_accession: str
    event_source: str
    candidate_event_count: int
    candidate_pair_count: int
    external_coupling_pair_count: int
    seed_pair_count: int
    final_pair_count: int
    final_long_range_pair_count: int
    final_contact_map_hash: str
    iteration_count: int
    decision_rule: str
    sequence_physical_prior_kind: str
    self_consistency_status: str
    self_consistent_internal_claim_allowed: bool
    external_independent_claim_allowed: bool
    global_folding_claim_allowed: bool
    folding_problem_solved: bool
    claim_rejection_reason: str
    negative_control_count: int
    real_rank_among_controls: int
    real_self_consistency_strength: float
    best_control_self_consistency_strength: float
    control_margin: float
    final_stability: float
    seed_retention: float
    coupling_retention: float
    graph_closure_coherence: float
    mean_final_contact_energy_score: float
    mean_final_secondary_structure_score: float
    mean_final_degree_consistency_score: float
    mean_final_physical_prior_score: float
    seed_contact_precision_after_native_audit: float
    seed_contact_recall_after_native_audit: float
    seed_long_range_recall_after_native_audit: float
    contact_precision_delta_vs_seed_after_native_audit: float
    long_range_recall_delta_vs_seed_after_native_audit: float
    contact_precision_after_native_audit: float
    contact_recall_after_native_audit: float
    long_range_precision_after_native_audit: float
    long_range_recall_after_native_audit: float
    native_contact_count_after_audit: int
    native_long_range_contact_count_after_audit: int
    true_positive_contacts_after_audit: int
    true_positive_long_range_contacts_after_audit: int
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    native_truth_attached_after_selection_for_evaluation: bool = True
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class SelfConsistentContactLoopPacket:
    report: SelfConsistentContactLoopReport
    iterations: tuple[SelfConsistentIteration, ...]
    controls: tuple[SelfConsistentControlResult, ...]
    selected_pairs: tuple[ContactPair, ...]
    self_evidence: tuple[IndependentContactEvidencePair, ...]
    metric: ContactMetricPacket


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, float(value))), 6)


def _score(value: float) -> float:
    return round(float(value), 6)


def _mean(values: Sequence[float]) -> float:
    return _score(mean(values)) if values else 0.0


def _jaccard(left: set[ContactPair], right: set[ContactPair]) -> float:
    if not left and not right:
        return 1.0
    return _rounded(len(left & right) / len(left | right))


def _candidate_pairs_from_events(
    row: RealCoordinateVisualRow,
    events: Sequence[NucleusClosureEvent],
) -> tuple[
    set[ContactPair],
    dict[ContactPair, float],
    dict[ContactPair, int],
]:
    candidate_pairs: set[ContactPair] = set()
    event_scores: dict[ContactPair, float] = {}
    event_counts: dict[ContactPair, int] = defaultdict(int)
    for event in events:
        if event.row_id != row.row_id:
            continue
        for raw_pair in event.candidate_region_pairs():
            if raw_pair[0] < 1 or raw_pair[1] > row.sequence_length:
                continue
            if raw_pair[1] - raw_pair[0] < MIN_SEQUENCE_SEPARATION:
                continue
            pair = (raw_pair[0], raw_pair[1])
            candidate_pairs.add(pair)
            event_scores[pair] = max(event_scores.get(pair, 0.0), float(event.nucleus_score))
            event_counts[pair] += 1
    return candidate_pairs, event_scores, dict(event_counts)


def _coupling_scores_by_candidate_pair(
    row: RealCoordinateVisualRow,
    constraints: Sequence[CouplingConstraint],
    candidate_pairs: set[ContactPair],
) -> dict[ContactPair, float]:
    scores: dict[ContactPair, float] = {}
    for constraint in constraints:
        if constraint.row_id != row.row_id:
            continue
        pair = constraint.pair()
        if pair not in candidate_pairs:
            continue
        if constraint.coordinate_truth_used_to_build_constraint:
            continue
        if constraint.native_truth_used_before_coupling_selection:
            continue
        if constraint.structure_model_used:
            continue
        scores[pair] = max(scores.get(pair, 0.0), float(constraint.confidence))
    return scores


def _rank_percentiles(values: Mapping[ContactPair, float]) -> dict[ContactPair, float]:
    positive_values = sorted({float(value) for value in values.values() if float(value) > 0.0})
    if not positive_values:
        return {pair: 0.0 for pair in values}
    denominator = max(1, len(positive_values))
    rank_by_value = {
        value: (index + 1) / denominator
        for index, value in enumerate(positive_values)
    }
    return {
        pair: _rounded(rank_by_value[float(value)]) if float(value) > 0.0 else 0.0
        for pair, value in values.items()
    }


def _gap_boundary(
    scores: Mapping[ContactPair, float],
    *,
    seed_size: int,
    minimum_selected_count: int = 0,
    enforce_seed_identity_envelope: bool = False,
) -> GapBoundary:
    ordered = sorted(scores.items(), key=lambda item: (-item[1], item[0][0], item[0][1]))
    positive_ordered = tuple(item for item in ordered if item[1] > 0.0)
    if not positive_ordered:
        return GapBoundary(
            selected_count=0,
            boundary_score=0.0,
            reason="empty_positive_signal_distribution",
            largest_gap=0.0,
            positive_score_count=0,
            seed_identity_envelope_count=seed_size,
        )
    if len(positive_ordered) == 1:
        return GapBoundary(
            selected_count=1,
            boundary_score=_score(positive_ordered[0][1]),
            reason="single_positive_signal_pair",
            largest_gap=0.0,
            positive_score_count=1,
            seed_identity_envelope_count=seed_size,
        )

    values = [score for _pair, score in positive_ordered]
    gaps = [values[index] - values[index + 1] for index in range(len(values) - 1)]
    minimum_count = max(0, min(int(minimum_selected_count), len(positive_ordered)))
    allowed_gap_indices = [
        index
        for index in range(len(gaps))
        if index + 1 >= minimum_count
    ]
    if not allowed_gap_indices:
        allowed_gap_indices = list(range(len(gaps)))
    if not allowed_gap_indices:
        cutoff = len(positive_ordered)
        largest_gap = 0.0
        reason = "all_positive_scores_tied"
    else:
        largest_gap_index = max(allowed_gap_indices, key=lambda index: gaps[index])
        largest_gap = gaps[largest_gap_index]
        if largest_gap <= 0.0:
            cutoff = len(positive_ordered)
            reason = "all_positive_scores_tied"
        else:
            cutoff = largest_gap_index + 1
            reason = "row_local_largest_gap"
            if minimum_count > 0:
                reason = f"{reason}_after_adaptive_seed_derived_minimum"

    if enforce_seed_identity_envelope and seed_size > 0 and cutoff < seed_size:
        cutoff = min(seed_size, len(positive_ordered))
        reason = f"{reason}_with_seed_identity_envelope"

    return GapBoundary(
        selected_count=cutoff,
        boundary_score=_score(positive_ordered[cutoff - 1][1]),
        reason=reason,
        largest_gap=_score(largest_gap),
        positive_score_count=len(positive_ordered),
        seed_identity_envelope_count=seed_size,
    )


def _adjacency(pairs: set[ContactPair]) -> dict[int, set[int]]:
    adjacency: dict[int, set[int]] = defaultdict(set)
    for left, right in pairs:
        adjacency[left].add(right)
        adjacency[right].add(left)
    return adjacency


def _degree_consistency_score_from_degrees(left_degree: int, right_degree: int) -> float:
    def one(degree: int) -> float:
        if 2 <= degree <= 5:
            return 1.0
        if degree == 1:
            return 0.70
        if degree == 6:
            return 0.82
        if degree == 7:
            return 0.62
        if degree == 8:
            return 0.42
        if degree == 9:
            return 0.25
        if degree >= 10:
            return 0.08
        return 0.45

    return _rounded(mean((one(left_degree), one(right_degree))))


def _summarize_prior_score_maps(
    *,
    selected_pairs: set[ContactPair],
    contact_energy_scores: Mapping[ContactPair, float],
    secondary_structure_scores: Mapping[ContactPair, float],
    degree_consistency_scores: Mapping[ContactPair, float],
) -> dict[str, float]:
    if not selected_pairs:
        return {
            "mean_contact_energy_score": 0.0,
            "mean_secondary_structure_score": 0.0,
            "mean_degree_consistency_score": 0.0,
            "mean_physical_prior_score": 0.0,
        }
    physical_scores = [
        0.45 * contact_energy_scores[pair]
        + 0.30 * secondary_structure_scores[pair]
        + 0.25 * degree_consistency_scores[pair]
        for pair in selected_pairs
    ]
    return {
        "mean_contact_energy_score": _rounded(mean(contact_energy_scores[pair] for pair in selected_pairs)),
        "mean_secondary_structure_score": _rounded(mean(secondary_structure_scores[pair] for pair in selected_pairs)),
        "mean_degree_consistency_score": _rounded(mean(degree_consistency_scores[pair] for pair in selected_pairs)),
        "mean_physical_prior_score": _rounded(mean(physical_scores)),
    }


def _score_candidate_pairs(
    *,
    candidate_pairs: set[ContactPair],
    current_pairs: set[ContactPair],
    coupling_scores: Mapping[ContactPair, float],
    event_scores: Mapping[ContactPair, float],
    event_counts: Mapping[ContactPair, int],
    static_physical_priors: Mapping[ContactPair, SequenceContactPhysicalPrior],
) -> tuple[dict[ContactPair, float], dict[str, dict[ContactPair, float]]]:
    adjacency = _adjacency(current_pairs)
    shared_neighbor_raw = {
        pair: float(len(adjacency[pair[0]] & adjacency[pair[1]]))
        for pair in candidate_pairs
    }
    prior_presence_raw = {
        pair: 1.0 if pair in current_pairs else 0.0
        for pair in candidate_pairs
    }
    direct_coupling_raw = {
        pair: float(coupling_scores.get(pair, 0.0))
        for pair in candidate_pairs
    }
    event_score_raw = {
        pair: float(event_scores.get(pair, 0.0))
        for pair in candidate_pairs
    }
    event_count_raw = {
        pair: float(event_counts.get(pair, 0))
        for pair in candidate_pairs
    }
    contact_energy_raw = {
        pair: float(static_physical_priors[pair].contact_energy_score)
        for pair in candidate_pairs
    }
    secondary_structure_raw = {
        pair: float(static_physical_priors[pair].secondary_structure_score)
        for pair in candidate_pairs
    }
    degree_consistency_raw = {}
    for pair in candidate_pairs:
        left, right = pair
        left_degree = len(adjacency[left]) + (0 if right in adjacency[left] else 1)
        right_degree = len(adjacency[right]) + (0 if left in adjacency[right] else 1)
        degree_consistency_raw[pair] = _degree_consistency_score_from_degrees(left_degree, right_degree)

    rank_layers = {
        "direct_coupling": _rank_percentiles(direct_coupling_raw),
        "prior_presence": _rank_percentiles(prior_presence_raw),
        "event_score": _rank_percentiles(event_score_raw),
        "event_count": _rank_percentiles(event_count_raw),
        "shared_neighbor": _rank_percentiles(shared_neighbor_raw),
        "contact_energy": _rank_percentiles(contact_energy_raw),
        "secondary_structure": _rank_percentiles(secondary_structure_raw),
        "degree_consistency": _rank_percentiles(degree_consistency_raw),
    }
    weights = {
        "direct_coupling": 0.42,
        "prior_presence": 0.20,
        "event_score": 0.08,
        "event_count": 0.04,
        "shared_neighbor": 0.06,
        "contact_energy": 0.12,
        "secondary_structure": 0.04,
        "degree_consistency": 0.04,
    }
    scores = {
        pair: _score(
            sum(weights[layer_name] * rank_layers[layer_name][pair] for layer_name in weights)
        )
        for pair in candidate_pairs
    }
    prior_score_maps = {
        "contact_energy": contact_energy_raw,
        "secondary_structure": secondary_structure_raw,
        "degree_consistency": degree_consistency_raw,
    }
    return scores, prior_score_maps


def _select_by_boundary(
    scores: Mapping[ContactPair, float],
    *,
    boundary: GapBoundary,
) -> set[ContactPair]:
    if boundary.selected_count <= 0:
        return set()
    ordered = sorted(scores.items(), key=lambda item: (-item[1], item[0][0], item[0][1]))
    positive_ordered = [item for item in ordered if item[1] > 0.0]
    return {pair for pair, _score_value in positive_ordered[: boundary.selected_count]}


def _loop_strength_components(
    *,
    final_pairs: set[ContactPair],
    previous_pairs: set[ContactPair],
    seed_pairs: set[ContactPair],
    coupling_scores: Mapping[ContactPair, float],
    boundary: GapBoundary,
) -> dict[str, float]:
    adjacency = _adjacency(final_pairs)
    closure_values = [
        len(adjacency[left] & adjacency[right])
        for left, right in final_pairs
    ]
    max_closure = max(closure_values) if closure_values else 0
    graph_closure = (
        _rounded(mean(value / max_closure for value in closure_values))
        if max_closure > 0
        else 0.0
    )
    seed_retention = _rounded(len(final_pairs & seed_pairs) / len(seed_pairs)) if seed_pairs else 0.0
    coupling_pairs = set(coupling_scores)
    coupling_retention = (
        _rounded(len(final_pairs & coupling_pairs) / len(coupling_pairs))
        if coupling_pairs
        else 0.0
    )
    final_stability = _jaccard(final_pairs, previous_pairs)
    boundary_gap = _rounded(boundary.largest_gap)
    strength = _rounded(
        mean(
            (
                final_stability,
                seed_retention,
                coupling_retention,
                graph_closure,
                boundary_gap,
            )
        )
    )
    return {
        "self_consistency_strength": strength,
        "final_stability": final_stability,
        "seed_retention": seed_retention,
        "coupling_retention": coupling_retention,
        "graph_closure_coherence": graph_closure,
        "boundary_gap": boundary_gap,
    }


def _run_self_loop(
    *,
    loop_name: str,
    row: RealCoordinateVisualRow,
    candidate_pairs: set[ContactPair],
    seed_pairs: set[ContactPair],
    coupling_scores: Mapping[ContactPair, float],
    event_scores: Mapping[ContactPair, float],
    event_counts: Mapping[ContactPair, int],
    static_physical_priors: Mapping[ContactPair, SequenceContactPhysicalPrior],
) -> tuple[set[ContactPair], tuple[SelfConsistentIteration, ...], dict[str, float]]:
    current_pairs = set(seed_pairs)
    previous_pairs = set(seed_pairs)
    iterations: list[SelfConsistentIteration] = []
    seen_hashes: set[str] = set()
    iteration_budget = max(1, int(sqrt(max(1, row.sequence_length))))
    last_boundary = GapBoundary(0, 0.0, "not_started", 0.0, 0, len(seed_pairs))

    for iteration_index in range(1, iteration_budget + 1):
        scores, prior_score_maps = _score_candidate_pairs(
            candidate_pairs=candidate_pairs,
            current_pairs=current_pairs,
            coupling_scores=coupling_scores,
            event_scores=event_scores,
            event_counts=event_counts,
            static_physical_priors=static_physical_priors,
        )
        boundary = _gap_boundary(
            scores,
            seed_size=len(seed_pairs),
            minimum_selected_count=max(1, len(seed_pairs) // 4),
            enforce_seed_identity_envelope=False,
        )
        selected_pairs = _select_by_boundary(scores, boundary=boundary)
        prior_summary = _summarize_prior_score_maps(
            selected_pairs=selected_pairs,
            contact_energy_scores=prior_score_maps["contact_energy"],
            secondary_structure_scores=prior_score_maps["secondary_structure"],
            degree_consistency_scores=prior_score_maps["degree_consistency"],
        )
        selected_hash = contact_map_hash(selected_pairs)
        changed = len(selected_pairs ^ current_pairs)
        previous_jaccard = _jaccard(selected_pairs, current_pairs)
        mean_score = _mean([scores[pair] for pair in selected_pairs])
        iterations.append(
            SelfConsistentIteration(
                loop_name=loop_name,
                iteration_index=iteration_index,
                candidate_pair_count=len(candidate_pairs),
                seed_pair_count=len(seed_pairs),
                selected_pair_count=len(selected_pairs),
                previous_selected_pair_count=len(current_pairs),
                changed_pair_count=changed,
                previous_jaccard=previous_jaccard,
                selected_contact_map_hash=selected_hash,
                boundary=boundary,
                mean_selected_score=mean_score,
                mean_selected_contact_energy_score=prior_summary["mean_contact_energy_score"],
                mean_selected_secondary_structure_score=prior_summary["mean_secondary_structure_score"],
                mean_selected_degree_consistency_score=prior_summary["mean_degree_consistency_score"],
                mean_selected_physical_prior_score=prior_summary["mean_physical_prior_score"],
            )
        )
        previous_pairs = current_pairs
        current_pairs = selected_pairs
        last_boundary = boundary
        if selected_hash in seen_hashes or changed == 0:
            break
        seen_hashes.add(selected_hash)

    final_scores, final_prior_score_maps = _score_candidate_pairs(
        candidate_pairs=candidate_pairs,
        current_pairs=current_pairs,
        coupling_scores=coupling_scores,
        event_scores=event_scores,
        event_counts=event_counts,
        static_physical_priors=static_physical_priors,
    )
    _ = final_scores
    final_prior_summary = _summarize_prior_score_maps(
        selected_pairs=current_pairs,
        contact_energy_scores=final_prior_score_maps["contact_energy"],
        secondary_structure_scores=final_prior_score_maps["secondary_structure"],
        degree_consistency_scores=final_prior_score_maps["degree_consistency"],
    )
    components = _loop_strength_components(
        final_pairs=current_pairs,
        previous_pairs=previous_pairs,
        seed_pairs=seed_pairs,
        coupling_scores=coupling_scores,
        boundary=last_boundary,
    )
    components.update(final_prior_summary)
    return current_pairs, tuple(iterations), components


def _candidate_pairs_by_separation(candidate_pairs: set[ContactPair]) -> dict[int, tuple[ContactPair, ...]]:
    buckets: dict[int, list[ContactPair]] = defaultdict(list)
    for pair in sorted(candidate_pairs):
        buckets[pair[1] - pair[0]].append(pair)
    return {separation: tuple(pairs) for separation, pairs in buckets.items()}


def _stable_index(*parts: object, modulus: int) -> int:
    encoded = ":".join(str(part) for part in parts).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest(), 16) % modulus


def _matched_control_seed_and_scores(
    *,
    row: RealCoordinateVisualRow,
    control_index: int,
    candidate_pairs: set[ContactPair],
    real_coupling_scores: Mapping[ContactPair, float],
) -> tuple[set[ContactPair], dict[ContactPair, float]]:
    buckets = _candidate_pairs_by_separation(candidate_pairs)
    selected: set[ContactPair] = set()
    control_scores: dict[ContactPair, float] = {}
    real_items = sorted(real_coupling_scores.items(), key=lambda item: (-item[1], item[0][0], item[0][1]))
    real_pairs = set(real_coupling_scores)
    for pair_index, (real_pair, confidence) in enumerate(real_items):
        separation = real_pair[1] - real_pair[0]
        bucket = buckets.get(separation, ())
        if not bucket:
            continue
        start = _stable_index(
            row.row_id,
            "self_consistent_control",
            control_index,
            pair_index,
            real_pair[0],
            real_pair[1],
            modulus=len(bucket),
        )
        chosen: ContactPair | None = None
        for offset in range(len(bucket)):
            candidate = bucket[(start + offset) % len(bucket)]
            if candidate in real_pairs or candidate in selected:
                continue
            chosen = candidate
            break
        if chosen is None:
            continue
        selected.add(chosen)
        control_scores[chosen] = float(confidence)
    return selected, control_scores


def _negative_control_count(seed_pair_count: int) -> int:
    # Bootstrap breadth is derived from the available seed evidence itself.
    # It is a runtime/audit budget, not an acceptance threshold.
    return max(1, int(sqrt(sqrt(max(1, seed_pair_count)))))


def _decision_from_controls(
    *,
    real_strength: float,
    controls: Sequence[SelfConsistentControlResult],
) -> tuple[str, bool, int, float, float]:
    scored = [("real", real_strength)] + [
        (control.control_id, control.self_consistency_strength)
        for control in controls
    ]
    scored.sort(key=lambda item: (-item[1], item[0]))
    real_rank = next(index + 1 for index, item in enumerate(scored) if item[0] == "real")
    best_control = max((control.self_consistency_strength for control in controls), default=0.0)
    margin = _score(real_strength - best_control)
    if not controls:
        return "rejected_no_negative_controls", False, real_rank, best_control, margin
    if real_rank != 1:
        return "rejected_by_matched_negative_control_rank", False, real_rank, best_control, margin
    values = [score for _label, score in scored]
    if len(values) == 1:
        return "self_consistency_survived_without_control_comparator", True, real_rank, best_control, margin
    gaps = [values[index] - values[index + 1] for index in range(len(values) - 1)]
    largest_gap_index = gaps.index(max(gaps)) if gaps else 0
    if values[0] == values[1]:
        return "rejected_real_control_tie", False, real_rank, best_control, margin
    if largest_gap_index != 0:
        return "rejected_real_not_separated_by_largest_rank_gap", False, real_rank, best_control, margin
    return "self_consistency_survived_matched_negative_controls", True, real_rank, best_control, margin


def self_consistent_evidence_from_pairs(
    *,
    row: RealCoordinateVisualRow,
    pairs: Sequence[ContactPair],
    confidence: float,
    source_id: str = "self_consistent_contact_loop_v0",
) -> tuple[IndependentContactEvidencePair, ...]:
    return tuple(
        IndependentContactEvidencePair(
            row_id=row.row_id,
            source_accession=row.source_accession,
            source_id=source_id,
            source_kind=SELF_CONSISTENT_INTERNAL_CONTACT_SOURCE_KIND,
            source_family="self_consistent_internal",
            i=pair[0],
            j=pair[1],
            confidence=_rounded(confidence),
            sequence_separation=pair[1] - pair[0],
            coordinate_truth_used_before_selection=False,
            native_truth_used_before_selection=False,
            raw_sequence_exposed=False,
        )
        for pair in normalized_contact_pairs(pairs)
    )


def run_self_consistent_contact_loop(
    *,
    row: RealCoordinateVisualRow,
    events: Sequence[NucleusClosureEvent],
    constraints: Sequence[CouplingConstraint],
    event_source: str,
) -> SelfConsistentContactLoopPacket:
    candidate_pairs, event_scores, event_counts = _candidate_pairs_from_events(row, events)
    coupling_scores = _coupling_scores_by_candidate_pair(row, constraints, candidate_pairs)
    seed_pairs = set(coupling_scores)
    static_physical_priors = build_sequence_physical_prior_scores(
        row=row,
        candidate_pairs=candidate_pairs,
        current_pairs=(),
    )

    final_pairs, iterations, components = _run_self_loop(
        loop_name="real",
        row=row,
        candidate_pairs=candidate_pairs,
        seed_pairs=seed_pairs,
        coupling_scores=coupling_scores,
        event_scores=event_scores,
        event_counts=event_counts,
        static_physical_priors=static_physical_priors,
    )

    controls: list[SelfConsistentControlResult] = []
    for control_index in range(1, _negative_control_count(len(seed_pairs)) + 1):
        control_seed, control_scores = _matched_control_seed_and_scores(
            row=row,
            control_index=control_index,
            candidate_pairs=candidate_pairs,
            real_coupling_scores=coupling_scores,
        )
        control_final, control_iterations, control_components = _run_self_loop(
            loop_name=f"control_{control_index:02d}",
            row=row,
            candidate_pairs=candidate_pairs,
            seed_pairs=control_seed,
            coupling_scores=control_scores,
            event_scores=event_scores,
            event_counts=event_counts,
            static_physical_priors=static_physical_priors,
        )
        controls.append(
            SelfConsistentControlResult(
                control_id=f"matched_control_{control_index:02d}",
                seed_pair_count=len(control_seed),
                final_pair_count=len(control_final),
                final_contact_map_hash=contact_map_hash(control_final),
                self_consistency_strength=control_components["self_consistency_strength"],
                final_stability=control_components["final_stability"],
                seed_retention=control_components["seed_retention"],
                coupling_retention=control_components["coupling_retention"],
                graph_closure_coherence=control_components["graph_closure_coherence"],
                boundary_gap=control_components["boundary_gap"],
                mean_final_contact_energy_score=control_components["mean_contact_energy_score"],
                mean_final_secondary_structure_score=control_components["mean_secondary_structure_score"],
                mean_final_degree_consistency_score=control_components["mean_degree_consistency_score"],
                mean_final_physical_prior_score=control_components["mean_physical_prior_score"],
            )
        )
        # The full per-control iteration trace is intentionally not merged into
        # the main report.  The control summary carries the decision evidence;
        # the main trace remains readable and stable.
        _ = control_iterations

    status, allowed, real_rank, best_control, margin = _decision_from_controls(
        real_strength=components["self_consistency_strength"],
        controls=controls,
    )

    selected_pairs = normalized_contact_pairs(final_pairs)
    native_pairs = set(row.native_contact_pairs())
    native_long_pairs = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
    selected_long_pairs = {pair for pair in selected_pairs if pair[1] - pair[0] >= 24}
    long_tp = selected_long_pairs & native_long_pairs
    metric = evaluate_contact_prediction(
        native_pairs=native_pairs,
        predicted_pairs=selected_pairs,
    )
    evidence = self_consistent_evidence_from_pairs(
        row=row,
        pairs=selected_pairs,
        confidence=components["self_consistency_strength"],
    )
    seed_metric = evaluate_contact_prediction(
        native_pairs=native_pairs,
        predicted_pairs=seed_pairs,
    )
    claim_rejection = (
        "global_folding_claim_rejected_self_source_is_not_external_independent"
        if allowed
        else f"global_folding_claim_rejected:{status}"
    )
    report = SelfConsistentContactLoopReport(
        kind=SELF_CONSISTENT_CONTACT_LOOP_KIND,
        row_id=row.row_id,
        source_accession=row.source_accession,
        event_source=event_source,
        candidate_event_count=len([event for event in events if event.row_id == row.row_id]),
        candidate_pair_count=len(candidate_pairs),
        external_coupling_pair_count=len(coupling_scores),
        seed_pair_count=len(seed_pairs),
        final_pair_count=len(selected_pairs),
        final_long_range_pair_count=len(selected_long_pairs),
        final_contact_map_hash=contact_map_hash(selected_pairs),
        iteration_count=len(iterations),
        decision_rule=SELF_CONSISTENT_SCORING_RULE,
        sequence_physical_prior_kind=SEQUENCE_PHYSICAL_PRIOR_KIND,
        self_consistency_status=status,
        self_consistent_internal_claim_allowed=allowed,
        external_independent_claim_allowed=False,
        global_folding_claim_allowed=False,
        folding_problem_solved=False,
        claim_rejection_reason=claim_rejection,
        negative_control_count=len(controls),
        real_rank_among_controls=real_rank,
        real_self_consistency_strength=components["self_consistency_strength"],
        best_control_self_consistency_strength=best_control,
        control_margin=margin,
        final_stability=components["final_stability"],
        seed_retention=components["seed_retention"],
        coupling_retention=components["coupling_retention"],
        graph_closure_coherence=components["graph_closure_coherence"],
        mean_final_contact_energy_score=components["mean_contact_energy_score"],
        mean_final_secondary_structure_score=components["mean_secondary_structure_score"],
        mean_final_degree_consistency_score=components["mean_degree_consistency_score"],
        mean_final_physical_prior_score=components["mean_physical_prior_score"],
        seed_contact_precision_after_native_audit=seed_metric.native_contact_precision,
        seed_contact_recall_after_native_audit=seed_metric.native_contact_recall,
        seed_long_range_recall_after_native_audit=seed_metric.long_range_contact_recall,
        contact_precision_delta_vs_seed_after_native_audit=_score(metric.native_contact_precision - seed_metric.native_contact_precision),
        long_range_recall_delta_vs_seed_after_native_audit=_score(metric.long_range_contact_recall - seed_metric.long_range_contact_recall),
        contact_precision_after_native_audit=metric.native_contact_precision,
        contact_recall_after_native_audit=metric.native_contact_recall,
        long_range_precision_after_native_audit=(
            _rounded(len(long_tp) / len(selected_long_pairs)) if selected_long_pairs else 0.0
        ),
        long_range_recall_after_native_audit=(
            _rounded(len(long_tp) / len(native_long_pairs)) if native_long_pairs else 0.0
        ),
        native_contact_count_after_audit=len(native_pairs),
        native_long_range_contact_count_after_audit=len(native_long_pairs),
        true_positive_contacts_after_audit=metric.true_positive_contacts,
        true_positive_long_range_contacts_after_audit=len(long_tp),
    )
    return SelfConsistentContactLoopPacket(
        report=report,
        iterations=iterations,
        controls=tuple(controls),
        selected_pairs=selected_pairs,
        self_evidence=evidence,
        metric=metric,
    )
