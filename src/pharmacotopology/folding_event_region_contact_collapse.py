from __future__ import annotations

from dataclasses import asdict, dataclass
from statistics import mean
from typing import Mapping, Sequence

from pharmacotopology.folding_contact_law_features import ContactLawFeatureRow
from pharmacotopology.folding_evolutionary_constraints import CouplingConstraint
from pharmacotopology.folding_native_contact_eval import ContactMetricPacket, evaluate_contact_prediction
from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent
from pharmacotopology.folding_real_coordinate_visual_benchmark import RealCoordinateVisualRow


EVENT_REGION_CONTACT_COLLAPSE_KIND = "event_region_contact_collapse_v0"
EVENT_REGION_CONTACT_COLLAPSE_BOUNDARY = (
    "sequence_and_external_coupling_only_before_native_evaluation"
)

# The coarse frontier layer works on 8x8 regions.  The collapse layer must not
# keep all 64 candidate residue pairs; it should keep a small, ranked subset and
# report the precision/recall tradeoff explicitly.
DEFAULT_PRECISION_MAX_PAIRS_PER_EVENT = 5
DEFAULT_RECALL_MAX_PAIRS_PER_EVENT = 16
DEFAULT_BALANCED_PAIRS_PER_EVENT = 6
DEFAULT_INTERNAL_GAP_CORE_PAIRS_PER_EVENT = 5
DEFAULT_INTERNAL_GAP_MIN_FIRST_LOCK_SCORE = 0.72
DEFAULT_INTERNAL_GAP_FIRST_LOCK_DELTA = 0.035
DEFAULT_INTERNAL_GAP_DENSE_PAIR_MIN = 40
DEFAULT_INTERNAL_GAP_SHORT_SEGMENT_CLOSE_MAX = 7
DEFAULT_INTERNAL_GAP_CLOSURE_SCORE_MIN = 0.65
DEFAULT_INTERNAL_GAP_EDGE_RESCUE_SCORE_MIN = 0.57
DEFAULT_RESIDUE_DEGREE_CAP = 5


@dataclass(frozen=True)
class CollapsedContactPair:
    row_id: str
    source_accession: str
    event_id: str
    i: int
    j: int
    collapse_strategy: str
    collapse_rank: int
    collapse_score: float
    coupling_density_score: float
    sequence_law_support_score: float
    ridge_coherence_score: float
    boundary_coherence_score: float
    degree_i_before_selection: int
    degree_j_before_selection: int
    selected_by_gap: bool
    selected_by_degree_cap: bool
    collapse_selection_reason: str = "ranked_core"
    native_truth_used_before_collapse_selection: bool = False
    coordinate_truth_used_before_collapse_selection: bool = False
    raw_sequence_exposed: bool = False

    def pair(self) -> tuple[int, int]:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class EventRegionCollapseSummary:
    row_id: str
    source_accession: str
    event_id: str
    collapse_strategy: str
    candidate_region_pair_count: int
    selected_pair_count: int
    max_pairs_per_event: int
    natural_gap_cutoff: int
    applied_cutoff: int
    direct_coupling_count_in_region: int
    density_supported_pair_count_in_region: int
    collapse_decision_reason: str
    first_score_gap: float
    support_completion_pair_count: int
    edge_rescue_pair_count: int
    mean_coupling_density_score: float
    mean_sequence_law_support_score: float
    mean_ridge_coherence_score: float
    mean_boundary_coherence_score: float
    native_truth_used_before_collapse_selection: bool = False
    coordinate_truth_used_before_collapse_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RowCollapseEvaluation:
    row_id: str
    source_accession: str
    collapse_strategy: str
    selected_event_count: int
    uncollapsed_region_pair_count: int
    collapsed_pair_count: int
    collapse_reduction_ratio: float
    uncollapsed_true_positive_contacts: int
    collapsed_true_positive_contacts: int
    uncollapsed_region_precision: float
    collapsed_contact_precision: float
    uncollapsed_region_recall: float
    collapsed_contact_recall: float
    uncollapsed_long_range_recall: float
    collapsed_long_range_recall: float
    uncollapsed_long_range_precision: float
    collapsed_long_range_precision: float
    uncollapsed_long_range_f1: float
    collapsed_long_range_f1: float
    frontier_native_pair_count: int
    frontier_native_retention: float
    frontier_long_native_pair_count: int
    frontier_long_native_retention: float
    native_contact_count: int
    native_long_range_contact_count: int
    precision_improvement_factor: float
    native_truth_used_before_collapse_selection: bool = False
    coordinate_truth_used_before_collapse_selection: bool = False
    native_truth_attached_after_collapse_for_evaluation: bool = True
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RowCollapseResult:
    evaluation: RowCollapseEvaluation
    collapsed_pairs: tuple[CollapsedContactPair, ...]
    event_summaries: tuple[EventRegionCollapseSummary, ...]
    uncollapsed_pairs: tuple[tuple[int, int], ...]
    collapsed_metric_packet: ContactMetricPacket
    uncollapsed_metric_packet: ContactMetricPacket


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _score(value: float) -> float:
    return round(float(value), 6)


def _mean(values: Sequence[float]) -> float:
    return _score(mean(values)) if values else 0.0


def _f1(precision: float, recall: float) -> float:
    return _rounded(
        2.0 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )


def _precision(true_positive_count: int, predicted_count: int) -> float:
    if predicted_count <= 0:
        return 0.0
    return _rounded(true_positive_count / predicted_count)


def _features_by_pair(
    row_features: Sequence[ContactLawFeatureRow],
) -> dict[tuple[int, int], ContactLawFeatureRow]:
    return {feature.pair(): feature for feature in row_features}


def _constraints_by_pair(
    row_constraints: Sequence[CouplingConstraint],
) -> dict[tuple[int, int], CouplingConstraint]:
    return {constraint.pair(): constraint for constraint in row_constraints}


def _pair_score_map(
    *,
    features_by_pair: Mapping[tuple[int, int], ContactLawFeatureRow],
    constraints_by_pair: Mapping[tuple[int, int], CouplingConstraint],
) -> dict[tuple[int, int], float]:
    scores: dict[tuple[int, int], float] = {}
    for pair, feature in features_by_pair.items():
        constraint = constraints_by_pair.get(pair)
        coupling = constraint.confidence if constraint is not None else 0.0
        # This is an internal support field, not a native label.  It combines the
        # sequence-only law score with any external DCA confidence available for
        # that exact residue pair.
        scores[pair] = _rounded(
            0.52 * feature.pair_plus_cluster_plus_entropy_score
            + 0.26 * feature.pair_plus_entropy_score
            + 0.22 * coupling
        )
    return scores


def _neighbor_pairs(pair: tuple[int, int], *, radius: int = 2) -> tuple[tuple[int, int], ...]:
    i, j = pair
    output: list[tuple[int, int]] = []
    for left_delta in range(-radius, radius + 1):
        for right_delta in range(-radius, radius + 1):
            if left_delta == 0 and right_delta == 0:
                continue
            neighbor = (i + left_delta, j + right_delta)
            if neighbor[0] < 1 or neighbor[1] <= neighbor[0]:
                continue
            output.append(neighbor)
    return tuple(output)


def _coupling_density_score(
    pair: tuple[int, int],
    *,
    constraints_by_pair: Mapping[tuple[int, int], CouplingConstraint],
    row_max_confidence: float,
) -> float:
    if row_max_confidence <= 0.0:
        return 0.0
    exact = constraints_by_pair.get(pair)
    exact_confidence = exact.confidence if exact is not None else 0.0
    neighbor_confidences = [
        constraints_by_pair[neighbor].confidence
        for neighbor in _neighbor_pairs(pair, radius=2)
        if neighbor in constraints_by_pair
    ]
    neighbor_density = mean(neighbor_confidences) if neighbor_confidences else 0.0
    return _rounded((0.70 * exact_confidence + 0.30 * neighbor_density) / row_max_confidence)


def _sequence_law_support_score(feature: ContactLawFeatureRow) -> float:
    return _rounded(
        0.30 * feature.pair_plus_cluster_plus_entropy_score
        + 0.20 * feature.pair_plus_entropy_score
        + 0.16 * feature.hydrophobic_pair_support
        + 0.12 * feature.aromatic_anchor_support
        + 0.10 * feature.opposite_charge_support
        + 0.08 * feature.parallel_contact_support
        + 0.04 * max(feature.helix_window_support, feature.beta_window_support)
        - 0.12 * feature.same_charge_penalty
        - 0.10 * feature.breaker_penalty
        - 0.06 * feature.isolation_penalty
    )


def _ridge_coherence_score(
    pair: tuple[int, int],
    *,
    support_scores: Mapping[tuple[int, int], float],
) -> float:
    i, j = pair
    diagonal_scores: list[float] = []
    anti_diagonal_scores: list[float] = []
    for delta in range(-3, 4):
        if delta == 0:
            continue
        diagonal = (i + delta, j + delta)
        anti_diagonal = (i + delta, j - delta)
        if diagonal[1] > diagonal[0]:
            diagonal_scores.append(float(support_scores.get(diagonal, 0.0)))
        if anti_diagonal[1] > anti_diagonal[0]:
            anti_diagonal_scores.append(float(support_scores.get(anti_diagonal, 0.0)))
    return _rounded(max(_mean(diagonal_scores), _mean(anti_diagonal_scores)))


def _region_position_scores(
    pair: tuple[int, int],
    event: NucleusClosureEvent,
) -> dict[str, float]:
    left_denominator = max(1, event.segment_a_end - event.segment_a_start)
    right_denominator = max(1, event.segment_b_end - event.segment_b_start)
    left_fraction = (pair[0] - event.segment_a_start) / left_denominator
    right_fraction = (pair[1] - event.segment_b_start) / right_denominator
    entry_coherence = _rounded((1.0 - left_fraction + 1.0 - right_fraction) / 2.0)
    exit_coherence = _rounded((left_fraction + right_fraction) / 2.0)
    edge_coherence = _rounded(
        2.0 * max(abs(left_fraction - 0.5), abs(right_fraction - 0.5))
    )
    corner_coherence = _rounded(
        abs(left_fraction - 0.5) + abs(right_fraction - 0.5)
    )
    main_diagonal_coherence = _rounded(1.0 - abs(left_fraction - right_fraction))
    anti_diagonal_coherence = _rounded(
        1.0 - abs(left_fraction + right_fraction - 1.0)
    )
    return {
        "entry_coherence": entry_coherence,
        "exit_coherence": exit_coherence,
        "edge_coherence": edge_coherence,
        "corner_coherence": corner_coherence,
        "main_diagonal_coherence": main_diagonal_coherence,
        "anti_diagonal_coherence": anti_diagonal_coherence,
    }


def _boundary_coherence_score(
    pair: tuple[int, int],
    event: NucleusClosureEvent,
    feature: ContactLawFeatureRow,
) -> float:
    position = _region_position_scores(pair, event)
    low_regular_structure = 1.0 - max(
        feature.helix_window_support,
        feature.beta_window_support,
    )
    loop_boundary = 0.60 * feature.breaker_penalty + 0.40 * feature.loop_entropy_cost
    return _rounded(
        0.36 * position["entry_coherence"]
        + 0.18 * position["edge_coherence"]
        + 0.14 * position["corner_coherence"]
        + 0.18 * loop_boundary
        + 0.14 * low_regular_structure
    )


def _frontier_precision_score(
    *,
    pair: tuple[int, int],
    event: NucleusClosureEvent,
    feature: ContactLawFeatureRow,
    coupling_density_score: float,
    sequence_law_support_score: float,
    ridge_coherence_score: float,
    boundary_coherence_score: float,
) -> float:
    position = _region_position_scores(pair, event)
    low_pair_law = 1.0 - feature.pair_only_score
    loop_boundary = 0.60 * feature.breaker_penalty + 0.40 * feature.loop_entropy_cost
    low_regular_structure = 1.0 - max(
        feature.helix_window_support,
        feature.beta_window_support,
    )
    unsupported_sequence_peak = max(
        0.0,
        feature.cluster_neighbor_support - coupling_density_score,
    )
    # Precision mode attacks the observed collapse failure: chemically obvious
    # pair-level peaks inside a sparse frontier region are often false attractions.
    # It therefore keeps a small boundary-coherent subset instead of all 64 pairs.
    return _rounded(
        0.32 * position["entry_coherence"]
        + 0.18 * low_pair_law
        + 0.16 * loop_boundary
        + 0.12 * low_regular_structure
        + 0.08 * position["edge_coherence"]
        + 0.04 * position["corner_coherence"]
        + 0.06 * coupling_density_score
        + 0.04 * ridge_coherence_score
        - 0.08 * unsupported_sequence_peak
        - 0.05 * feature.hydrophobic_pair_support
        - 0.03 * sequence_law_support_score
    )


def _ridge_coupling_score(
    *,
    feature: ContactLawFeatureRow,
    coupling_density_score: float,
    sequence_law_support_score: float,
    ridge_coherence_score: float,
    boundary_coherence_score: float,
) -> float:
    # Coupling-dense regions should collapse around external-DCA support and
    # ridges.  This mode is conservative and does not use native labels.
    return _rounded(
        0.38 * coupling_density_score
        + 0.26 * ridge_coherence_score
        + 0.22 * sequence_law_support_score
        + 0.10 * boundary_coherence_score
        + 0.04 * (1.0 - feature.loop_entropy_cost)
    )


def _gap_cutoff(scores: Sequence[float]) -> int:
    if len(scores) <= 1:
        return len(scores)
    ordered = sorted(scores, reverse=True)
    gaps = [ordered[index] - ordered[index + 1] for index in range(len(ordered) - 1)]
    if not gaps:
        return len(scores)
    largest_gap_index = max(range(len(gaps)), key=lambda index: gaps[index])
    return largest_gap_index + 1


def _bounded_cutoff(
    natural_gap_cutoff: int,
    *,
    min_pairs_per_event: int,
    max_pairs_per_event: int,
) -> int:
    if max_pairs_per_event <= 0:
        return 0
    return max(min_pairs_per_event, min(max_pairs_per_event, natural_gap_cutoff))


def _first_score_gap(scored_rows: Sequence[Mapping[str, object]]) -> float:
    if len(scored_rows) < 2:
        return 0.0
    return _score(float(scored_rows[0]["score"]) - float(scored_rows[1]["score"]))


def _density_supported_pair_count(scored_rows: Sequence[Mapping[str, object]]) -> int:
    return sum(
        1
        for item in scored_rows
        if float(item.get("coupling_density_score", 0.0)) > 0.0
    )


def _internal_gap_balanced_cutoff(
    *,
    event: NucleusClosureEvent,
    scored_rows: Sequence[Mapping[str, object]],
    min_pairs_per_event: int,
    max_pairs_per_event: int,
) -> tuple[int, str]:
    if not scored_rows:
        return 0, "empty_region"
    segment_gap = event.segment_b_start - event.segment_a_end
    first_gap = _first_score_gap(scored_rows)
    top_score = float(scored_rows[0]["score"])
    density_count = _density_supported_pair_count(scored_rows)

    # Native-free event-level safety: a very short inter-segment region is not a
    # long-range contact-map claim, and a short sparse region is too noisy to
    # reopen after frontier selection.
    if segment_gap <= DEFAULT_INTERNAL_GAP_SHORT_SEGMENT_CLOSE_MAX:
        return 0, "short_segment_gap_closed"
    if segment_gap < 24 and density_count < DEFAULT_INTERNAL_GAP_DENSE_PAIR_MIN:
        return 0, "short_sparse_region_closed"

    # A decisive first drop in a long-range frontier event means the event has a
    # single high-confidence boundary contact.  This is an internal score-gap
    # rule; native labels are attached only after the selected set is frozen.
    if (
        segment_gap >= 24
        and top_score >= DEFAULT_INTERNAL_GAP_MIN_FIRST_LOCK_SCORE
        and first_gap >= DEFAULT_INTERNAL_GAP_FIRST_LOCK_DELTA
    ):
        return 1, "first_gap_lock"

    core_budget = min(
        DEFAULT_INTERNAL_GAP_CORE_PAIRS_PER_EVENT,
        max_pairs_per_event,
    )
    core_budget = max(min_pairs_per_event, core_budget)
    return core_budget, "five_pair_internal_gap_core"


def _event_position_inside(pair: tuple[int, int], event: NucleusClosureEvent) -> bool:
    return pair[0] > event.segment_a_start and pair[1] > event.segment_b_start


def _support_completion_candidates(
    *,
    core_pairs: set[tuple[int, int]],
    event: NucleusClosureEvent,
    scored_by_pair: Mapping[tuple[int, int], Mapping[str, object]],
) -> tuple[tuple[int, int], ...]:
    left_residues = {pair[0] for pair in core_pairs}
    right_residues = {pair[1] for pair in core_pairs}
    candidates: list[tuple[float, tuple[int, int]]] = []
    for i in left_residues:
        for j in right_residues:
            pair = (i, j)
            if pair in core_pairs or pair not in scored_by_pair:
                continue
            if not _event_position_inside(pair, event):
                continue
            score = float(scored_by_pair[pair]["score"])
            if score >= DEFAULT_INTERNAL_GAP_CLOSURE_SCORE_MIN:
                candidates.append((score, pair))
    candidates.sort(key=lambda item: (-item[0], item[1][0], item[1][1]))
    # Keep only the strongest closure.  This allows an internally supported
    # missing rectangle corner without turning a five-pair core back into a wide
    # 8x8 region.
    return tuple(pair for _, pair in candidates[:1])


def _edge_rescue_candidates(
    *,
    event: NucleusClosureEvent,
    scored_by_pair: Mapping[tuple[int, int], Mapping[str, object]],
    already_selected: set[tuple[int, int]],
    applied_cutoff: int,
) -> tuple[tuple[int, int], ...]:
    if applied_cutoff > 1:
        return ()
    if event.segment_b_start - event.segment_a_end < 24:
        return ()
    candidates: list[tuple[float, tuple[int, int]]] = []
    for pair, item in scored_by_pair.items():
        if pair in already_selected:
            continue
        feature = item.get("feature")
        if feature is None:
            continue
        left_offset = pair[0] - event.segment_a_start
        right_offset = pair[1] - event.segment_b_start
        score = float(item["score"])
        # Rescue one internal left-edge partner when a first-gap lock would
        # otherwise throw away a second coherent boundary contact.  The rule is
        # sequence/coupling-only: no native map is read here.
        if (
            right_offset == 0
            and 3 <= left_offset <= 5
            and getattr(feature, "breaker_penalty", 0.0) == 0.0
            and score >= DEFAULT_INTERNAL_GAP_EDGE_RESCUE_SCORE_MIN
        ):
            candidates.append((score, pair))
    candidates.sort(key=lambda item: (-item[0], item[1][0], item[1][1]))
    return tuple(pair for _, pair in candidates[:1])


def collapse_event_region_contacts(
    *,
    event: NucleusClosureEvent,
    row_features: Sequence[ContactLawFeatureRow],
    row_constraints: Sequence[CouplingConstraint],
    collapse_strategy: str = "frontier_precision",
    min_pairs_per_event: int = 1,
    max_pairs_per_event: int = DEFAULT_PRECISION_MAX_PAIRS_PER_EVENT,
    residue_degree_cap: int = DEFAULT_RESIDUE_DEGREE_CAP,
) -> tuple[tuple[CollapsedContactPair, ...], EventRegionCollapseSummary]:
    if collapse_strategy not in {"frontier_precision", "frontier_recall", "frontier_balanced", "frontier_internal_gap_balanced", "ridge_coupling"}:
        raise ValueError(f"unknown collapse strategy: {collapse_strategy}")
    features_by_pair = _features_by_pair(row_features)
    constraints_by_pair = _constraints_by_pair(row_constraints)
    support_scores = _pair_score_map(
        features_by_pair=features_by_pair,
        constraints_by_pair=constraints_by_pair,
    )
    row_max_confidence = max(
        (constraint.confidence for constraint in row_constraints),
        default=0.0,
    )
    candidate_pairs = tuple(event.candidate_region_pairs())
    scored_rows: list[dict[str, object]] = []
    for pair in candidate_pairs:
        feature = features_by_pair.get(pair)
        if feature is None:
            continue
        coupling_density = _coupling_density_score(
            pair,
            constraints_by_pair=constraints_by_pair,
            row_max_confidence=row_max_confidence,
        )
        sequence_law = _sequence_law_support_score(feature)
        ridge = _ridge_coherence_score(pair, support_scores=support_scores)
        boundary = _boundary_coherence_score(pair, event, feature)
        if collapse_strategy in {"frontier_precision", "frontier_recall", "frontier_balanced", "frontier_internal_gap_balanced"}:
            score = _frontier_precision_score(
                pair=pair,
                event=event,
                feature=feature,
                coupling_density_score=coupling_density,
                sequence_law_support_score=sequence_law,
                ridge_coherence_score=ridge,
                boundary_coherence_score=boundary,
            )
        else:
            score = _ridge_coupling_score(
                feature=feature,
                coupling_density_score=coupling_density,
                sequence_law_support_score=sequence_law,
                ridge_coherence_score=ridge,
                boundary_coherence_score=boundary,
            )
        scored_rows.append(
            {
                "pair": pair,
                "score": score,
                "coupling_density_score": coupling_density,
                "sequence_law_support_score": sequence_law,
                "ridge_coherence_score": ridge,
                "boundary_coherence_score": boundary,
                "feature": feature,
            }
        )

    scored_rows.sort(
        key=lambda item: (
            -float(item["score"]),
            int(item["pair"][0]),  # type: ignore[index]
            int(item["pair"][1]),  # type: ignore[index]
        )
    )
    natural_gap_cutoff = _gap_cutoff([float(item["score"]) for item in scored_rows])
    first_gap = _first_score_gap(scored_rows)
    density_supported_count = _density_supported_pair_count(scored_rows)
    collapse_decision_reason = "natural_gap"
    if collapse_strategy == "frontier_recall":
        max_pairs_per_event = max(max_pairs_per_event, DEFAULT_RECALL_MAX_PAIRS_PER_EVENT)
        min_pairs_per_event = max(min_pairs_per_event, min(8, max_pairs_per_event))
    if collapse_strategy == "frontier_balanced":
        max_pairs_per_event = max(max_pairs_per_event, DEFAULT_BALANCED_PAIRS_PER_EVENT)
        min_pairs_per_event = max(min_pairs_per_event, min(DEFAULT_BALANCED_PAIRS_PER_EVENT, max_pairs_per_event))
    if collapse_strategy == "frontier_recall":
        # Recall mode is an explicit tradeoff probe: keep the bounded upper slice
        # even when the score gap is early, so we can test whether the frontier
        # contained recoverable contacts without reverting to all 64 pairs.
        applied_cutoff = max_pairs_per_event
        collapse_decision_reason = "recall_upper_slice"
    elif collapse_strategy == "frontier_balanced":
        # Legacy balanced mode: a fixed six-pair event budget kept for regression
        # comparisons against the new native-free internal-gap strategy.
        applied_cutoff = _bounded_cutoff(
            DEFAULT_BALANCED_PAIRS_PER_EVENT,
            min_pairs_per_event=min_pairs_per_event,
            max_pairs_per_event=max_pairs_per_event,
        )
        collapse_decision_reason = "fixed_six_pair_budget"
    elif collapse_strategy == "frontier_internal_gap_balanced":
        applied_cutoff, collapse_decision_reason = _internal_gap_balanced_cutoff(
            event=event,
            scored_rows=scored_rows,
            min_pairs_per_event=min_pairs_per_event,
            max_pairs_per_event=max_pairs_per_event,
        )
    else:
        applied_cutoff = _bounded_cutoff(
            natural_gap_cutoff,
            min_pairs_per_event=min_pairs_per_event,
            max_pairs_per_event=max_pairs_per_event,
        )

    scored_by_pair = {
        item["pair"]: item  # type: ignore[index]
        for item in scored_rows
    }
    rank_by_pair = {
        item["pair"]: rank  # type: ignore[index]
        for rank, item in enumerate(scored_rows, start=1)
    }
    selected_item_rows: list[tuple[Mapping[str, object], str]] = []
    degree: dict[int, int] = {}

    def add_pair_item(item: Mapping[str, object], reason: str) -> bool:
        pair = item["pair"]  # type: ignore[assignment]
        i, j = int(pair[0]), int(pair[1])  # type: ignore[index]
        if any(existing[0]["pair"] == pair for existing in selected_item_rows):
            return False
        degree_i = degree.get(i, 0)
        degree_j = degree.get(j, 0)
        if degree_i >= residue_degree_cap or degree_j >= residue_degree_cap:
            return False
        selected_item_rows.append((item, reason))
        degree[i] = degree_i + 1
        degree[j] = degree_j + 1
        return True

    core_selected_count = 0
    for item in scored_rows:
        if core_selected_count >= applied_cutoff:
            break
        if add_pair_item(item, "ranked_core"):
            core_selected_count += 1

    support_completion_count = 0
    edge_rescue_count = 0
    if collapse_strategy == "frontier_internal_gap_balanced":
        core_pair_set = {
            item["pair"]  # type: ignore[index]
            for item, reason in selected_item_rows
            if reason == "ranked_core"
        }
        for pair in _support_completion_candidates(
            core_pairs=core_pair_set,
            event=event,
            scored_by_pair=scored_by_pair,
        ):
            if add_pair_item(scored_by_pair[pair], "support_completion"):
                support_completion_count += 1
        selected_pair_set = {
            item["pair"]  # type: ignore[index]
            for item, _ in selected_item_rows
        }
        for pair in _edge_rescue_candidates(
            event=event,
            scored_by_pair=scored_by_pair,
            already_selected=selected_pair_set,
            applied_cutoff=applied_cutoff,
        ):
            if add_pair_item(scored_by_pair[pair], "edge_rescue"):
                edge_rescue_count += 1

    selected: list[CollapsedContactPair] = []
    for item, reason in selected_item_rows:
        pair = item["pair"]  # type: ignore[assignment]
        i, j = int(pair[0]), int(pair[1])  # type: ignore[index]
        rank = rank_by_pair[pair]
        selected_by_gap = rank <= natural_gap_cutoff
        selected.append(
            CollapsedContactPair(
                row_id=event.row_id,
                source_accession=event.source_accession,
                event_id=event.event_id,
                i=i,
                j=j,
                collapse_strategy=collapse_strategy,
                collapse_rank=rank,
                collapse_score=_score(float(item["score"])),
                coupling_density_score=_score(float(item["coupling_density_score"])),
                sequence_law_support_score=_score(float(item["sequence_law_support_score"])),
                ridge_coherence_score=_score(float(item["ridge_coherence_score"])),
                boundary_coherence_score=_score(float(item["boundary_coherence_score"])),
                degree_i_before_selection=0,
                degree_j_before_selection=0,
                selected_by_gap=selected_by_gap,
                selected_by_degree_cap=True,
                collapse_selection_reason=reason,
            )
        )

    direct_coupling_count = sum(1 for pair in candidate_pairs if pair in constraints_by_pair)
    summary = EventRegionCollapseSummary(
        row_id=event.row_id,
        source_accession=event.source_accession,
        event_id=event.event_id,
        collapse_strategy=collapse_strategy,
        candidate_region_pair_count=len(candidate_pairs),
        selected_pair_count=len(selected),
        max_pairs_per_event=max_pairs_per_event,
        natural_gap_cutoff=natural_gap_cutoff,
        applied_cutoff=applied_cutoff,
        direct_coupling_count_in_region=direct_coupling_count,
        density_supported_pair_count_in_region=density_supported_count,
        collapse_decision_reason=collapse_decision_reason,
        first_score_gap=first_gap,
        support_completion_pair_count=support_completion_count,
        edge_rescue_pair_count=edge_rescue_count,
        mean_coupling_density_score=_mean(
            [float(item["coupling_density_score"]) for item in scored_rows]
        ),
        mean_sequence_law_support_score=_mean(
            [float(item["sequence_law_support_score"]) for item in scored_rows]
        ),
        mean_ridge_coherence_score=_mean(
            [float(item["ridge_coherence_score"]) for item in scored_rows]
        ),
        mean_boundary_coherence_score=_mean(
            [float(item["boundary_coherence_score"]) for item in scored_rows]
        ),
    )
    return tuple(selected), summary


def collapse_row_event_regions(
    *,
    row: RealCoordinateVisualRow,
    events: Sequence[NucleusClosureEvent],
    row_features: Sequence[ContactLawFeatureRow],
    row_constraints: Sequence[CouplingConstraint],
    collapse_strategy: str = "frontier_precision",
    min_pairs_per_event: int = 1,
    max_pairs_per_event: int = DEFAULT_PRECISION_MAX_PAIRS_PER_EVENT,
    residue_degree_cap: int = DEFAULT_RESIDUE_DEGREE_CAP,
) -> RowCollapseResult:
    collapsed_pairs: list[CollapsedContactPair] = []
    event_summaries: list[EventRegionCollapseSummary] = []
    uncollapsed_pairs: set[tuple[int, int]] = set()
    for event in events:
        uncollapsed_pairs.update(event.candidate_region_pairs())
        event_pairs, event_summary = collapse_event_region_contacts(
            event=event,
            row_features=row_features,
            row_constraints=row_constraints,
            collapse_strategy=collapse_strategy,
            min_pairs_per_event=min_pairs_per_event,
            max_pairs_per_event=max_pairs_per_event,
            residue_degree_cap=residue_degree_cap,
        )
        collapsed_pairs.extend(event_pairs)
        event_summaries.append(event_summary)

    collapsed_pair_set = {pair.pair() for pair in collapsed_pairs}
    native_pairs = set(row.native_contact_pairs())
    native_long_pairs = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
    frontier_native_pairs = uncollapsed_pairs & native_pairs
    frontier_long_native_pairs = uncollapsed_pairs & native_long_pairs
    uncollapsed_long_pairs = {pair for pair in uncollapsed_pairs if pair[1] - pair[0] >= 24}
    collapsed_long_pairs = {pair for pair in collapsed_pair_set if pair[1] - pair[0] >= 24}
    uncollapsed_true_positive_long_pairs = uncollapsed_long_pairs & native_long_pairs
    collapsed_true_positive_long_pairs = collapsed_long_pairs & native_long_pairs
    uncollapsed_long_precision = _precision(
        len(uncollapsed_true_positive_long_pairs),
        len(uncollapsed_long_pairs),
    )
    collapsed_long_precision = _precision(
        len(collapsed_true_positive_long_pairs),
        len(collapsed_long_pairs),
    )
    collapsed_metric = evaluate_contact_prediction(
        native_pairs=native_pairs,
        predicted_pairs=collapsed_pair_set,
    )
    uncollapsed_metric = evaluate_contact_prediction(
        native_pairs=native_pairs,
        predicted_pairs=uncollapsed_pairs,
    )
    uncollapsed_precision = uncollapsed_metric.native_contact_precision
    collapsed_precision = collapsed_metric.native_contact_precision
    evaluation = RowCollapseEvaluation(
        row_id=row.row_id,
        source_accession=row.source_accession,
        collapse_strategy=collapse_strategy,
        selected_event_count=len(events),
        uncollapsed_region_pair_count=len(uncollapsed_pairs),
        collapsed_pair_count=len(collapsed_pair_set),
        collapse_reduction_ratio=_rounded(
            1.0 - len(collapsed_pair_set) / max(1, len(uncollapsed_pairs))
        ),
        uncollapsed_true_positive_contacts=uncollapsed_metric.true_positive_contacts,
        collapsed_true_positive_contacts=collapsed_metric.true_positive_contacts,
        uncollapsed_region_precision=uncollapsed_precision,
        collapsed_contact_precision=collapsed_precision,
        uncollapsed_region_recall=uncollapsed_metric.native_contact_recall,
        collapsed_contact_recall=collapsed_metric.native_contact_recall,
        uncollapsed_long_range_recall=uncollapsed_metric.long_range_contact_recall,
        collapsed_long_range_recall=collapsed_metric.long_range_contact_recall,
        uncollapsed_long_range_precision=uncollapsed_long_precision,
        collapsed_long_range_precision=collapsed_long_precision,
        uncollapsed_long_range_f1=_f1(
            uncollapsed_long_precision,
            uncollapsed_metric.long_range_contact_recall,
        ),
        collapsed_long_range_f1=_f1(
            collapsed_long_precision,
            collapsed_metric.long_range_contact_recall,
        ),
        frontier_native_pair_count=len(frontier_native_pairs),
        frontier_native_retention=_rounded(
            len(collapsed_pair_set & frontier_native_pairs) / max(1, len(frontier_native_pairs))
        ),
        frontier_long_native_pair_count=len(frontier_long_native_pairs),
        frontier_long_native_retention=_rounded(
            len(collapsed_pair_set & frontier_long_native_pairs)
            / max(1, len(frontier_long_native_pairs))
        ),
        native_contact_count=len(native_pairs),
        native_long_range_contact_count=len(native_long_pairs),
        precision_improvement_factor=_score(
            collapsed_precision / max(uncollapsed_precision, 0.000001)
        ),
    )
    return RowCollapseResult(
        evaluation=evaluation,
        collapsed_pairs=tuple(
            sorted(
                collapsed_pairs,
                key=lambda pair: (
                    pair.row_id,
                    pair.event_id,
                    pair.collapse_rank,
                    pair.i,
                    pair.j,
                ),
            )
        ),
        event_summaries=tuple(event_summaries),
        uncollapsed_pairs=tuple(sorted(uncollapsed_pairs)),
        collapsed_metric_packet=collapsed_metric,
        uncollapsed_metric_packet=uncollapsed_metric,
    )
