from __future__ import annotations

"""Sequence-only five-axis physical challenge layer.

This module deliberately does **not** claim that protein folding is solved.
It adds the five missing axes called out in the current project diagnosis and
keeps the final claim gate locked unless the coordinate audit is strong enough:

* energy/free-energy support from residue chemistry,
* entropy retention/cost for loops and flexible spans,
* cooperative neighbourhood support,
* sequence-only environment/context support,
* a tiny deterministic ensemble proxy for compact/open/hinge states.

All features are computed before native coordinate labels are attached. Native
contacts are used only by the benchmark/audit functions after selection.
"""

from collections import defaultdict
from dataclasses import asdict, dataclass
from math import sqrt
from statistics import mean
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_contact_topology import BREAKERS
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
    POLAR_UNCHARGED,
    POSITIVE_CHARGED,
    contact_energy_kcal,
    contact_energy_score,
    predict_lightweight_secondary_structure,
    secondary_structure_pair_score,
)
from pharmacotopology.folding_topology import (
    AROMATIC_AMINO_ACIDS,
    HYDROPHOBIC_AMINO_ACIDS,
)


FIVE_AXIS_PHYSICS_KIND = "sequence_only_five_axis_physics_challenge_v0"
FIVE_AXIS_DECISION_KIND = "five_axis_free_energy_cooperative_dynamic_context_contact_v0"
FIVE_AXIS_CONTROL_KIND = "five_axis_separation_matched_negative_control_v0"
FIVE_AXIS_SCORING_RULE = (
    "sequence_only_energy_entropy_cooperativity_context_dynamic_ensemble;"
    "row_local_largest_gap_after_sqrt_length_floor;native_audit_after_selection_only"
)
FIVE_AXIS_CLAIM_RULE = (
    "universal_physical_law_claim_requires_all_rows_precision_ge_0_60_and_"
    "all_rows_long_range_recall_ge_0_60_and_mean_f1_ge_0_60_and_matched_control_margin_ge_0_15"
)


@dataclass(frozen=True)
class FiveAxisContactDecision:
    kind: str
    row_id: str
    source_accession: str
    sequence_hash: str
    sequence_length: int
    i: int
    j: int
    sequence_separation: int
    contact_energy_kcal: float
    energy_support_score: float
    entropy_retention_score: float
    cooperative_neighbour_score: float
    environmental_context_score: float
    compact_state_score: float
    open_state_score: float
    hinge_state_score: float
    dynamic_ensemble_score: float
    dynamic_state_ambiguity: float
    free_energy_support_score: float
    final_score: float
    selected: bool
    selection_reason: str
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    learned_prior_used_before_selection: bool = False
    msa_used_before_selection: bool = False
    template_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def pair(self) -> ContactPair:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class FiveAxisBoundary:
    selected_count: int
    boundary_score: float
    largest_gap: float
    positive_score_count: int
    selection_reason: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class FiveAxisControlResult:
    control_id: str
    selected_contact_count: int
    selected_contact_map_hash: str
    native_contact_precision_after_audit: float
    native_contact_recall_after_audit: float
    long_range_contact_recall_after_audit: float
    contact_map_f1_after_audit: float
    control_kind: str = FIVE_AXIS_CONTROL_KIND
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class FiveAxisRowReport:
    row_id: str
    source_accession: str
    sequence_hash: str
    sequence_length: int
    candidate_pair_count: int
    selected_contact_count: int
    selected_long_range_contact_count: int
    selected_contact_map_hash: str
    boundary: FiveAxisBoundary
    metric_after_native_audit: ContactMetricPacket
    mean_selected_energy_support_score: float
    mean_selected_entropy_retention_score: float
    mean_selected_cooperative_neighbour_score: float
    mean_selected_environmental_context_score: float
    mean_selected_dynamic_ensemble_score: float
    mean_selected_free_energy_support_score: float
    matched_control_count: int
    best_control_f1_after_audit: float
    best_control_long_range_recall_after_audit: float
    f1_margin_vs_best_control: float
    long_range_recall_margin_vs_best_control: float
    row_physical_law_claim_allowed: bool
    row_claim_rejection_reason: str
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    native_truth_attached_after_selection_for_evaluation: bool = True
    learned_prior_used_before_selection: bool = False
    msa_used_before_selection: bool = False
    template_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "row_id": self.row_id,
            "source_accession": self.source_accession,
            "sequence_hash": self.sequence_hash,
            "sequence_length": self.sequence_length,
            "candidate_pair_count": self.candidate_pair_count,
            "selected_contact_count": self.selected_contact_count,
            "selected_long_range_contact_count": self.selected_long_range_contact_count,
            "selected_contact_map_hash": self.selected_contact_map_hash,
            "boundary": self.boundary.to_dict(),
            "metric_after_native_audit": self.metric_after_native_audit.to_dict(),
            "mean_selected_energy_support_score": self.mean_selected_energy_support_score,
            "mean_selected_entropy_retention_score": self.mean_selected_entropy_retention_score,
            "mean_selected_cooperative_neighbour_score": self.mean_selected_cooperative_neighbour_score,
            "mean_selected_environmental_context_score": self.mean_selected_environmental_context_score,
            "mean_selected_dynamic_ensemble_score": self.mean_selected_dynamic_ensemble_score,
            "mean_selected_free_energy_support_score": self.mean_selected_free_energy_support_score,
            "matched_control_count": self.matched_control_count,
            "best_control_f1_after_audit": self.best_control_f1_after_audit,
            "best_control_long_range_recall_after_audit": self.best_control_long_range_recall_after_audit,
            "f1_margin_vs_best_control": self.f1_margin_vs_best_control,
            "long_range_recall_margin_vs_best_control": self.long_range_recall_margin_vs_best_control,
            "row_physical_law_claim_allowed": self.row_physical_law_claim_allowed,
            "row_claim_rejection_reason": self.row_claim_rejection_reason,
            "coordinate_truth_used_before_selection": self.coordinate_truth_used_before_selection,
            "native_truth_used_before_selection": self.native_truth_used_before_selection,
            "native_truth_attached_after_selection_for_evaluation": self.native_truth_attached_after_selection_for_evaluation,
            "learned_prior_used_before_selection": self.learned_prior_used_before_selection,
            "msa_used_before_selection": self.msa_used_before_selection,
            "template_used_before_selection": self.template_used_before_selection,
            "raw_sequence_exposed": self.raw_sequence_exposed,
        }


@dataclass(frozen=True)
class FiveAxisChallengePacket:
    kind: str
    row_count: int
    decision_rule: str
    claim_rule: str
    entropy_axis_included: bool
    cooperativity_axis_included: bool
    dynamics_axis_included: bool
    context_axis_included: bool
    independent_physics_axis_included: bool
    mean_native_contact_precision_after_audit: float
    mean_native_contact_recall_after_audit: float
    mean_long_range_contact_recall_after_audit: float
    mean_contact_map_f1_after_audit: float
    min_row_precision_after_audit: float
    min_row_long_range_recall_after_audit: float
    mean_f1_margin_vs_best_control: float
    mean_long_range_recall_margin_vs_best_control: float
    universal_physical_law_claim_allowed: bool
    folding_problem_solved: bool
    claim_rejection_reason: str
    rows: tuple[FiveAxisRowReport, ...]
    controls: tuple[FiveAxisControlResult, ...]
    decisions: tuple[FiveAxisContactDecision, ...]
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    learned_prior_used_before_selection: bool = False
    msa_used_before_selection: bool = False
    template_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "row_count": self.row_count,
            "decision_rule": self.decision_rule,
            "claim_rule": self.claim_rule,
            "entropy_axis_included": self.entropy_axis_included,
            "cooperativity_axis_included": self.cooperativity_axis_included,
            "dynamics_axis_included": self.dynamics_axis_included,
            "context_axis_included": self.context_axis_included,
            "independent_physics_axis_included": self.independent_physics_axis_included,
            "mean_native_contact_precision_after_audit": self.mean_native_contact_precision_after_audit,
            "mean_native_contact_recall_after_audit": self.mean_native_contact_recall_after_audit,
            "mean_long_range_contact_recall_after_audit": self.mean_long_range_contact_recall_after_audit,
            "mean_contact_map_f1_after_audit": self.mean_contact_map_f1_after_audit,
            "min_row_precision_after_audit": self.min_row_precision_after_audit,
            "min_row_long_range_recall_after_audit": self.min_row_long_range_recall_after_audit,
            "mean_f1_margin_vs_best_control": self.mean_f1_margin_vs_best_control,
            "mean_long_range_recall_margin_vs_best_control": self.mean_long_range_recall_margin_vs_best_control,
            "universal_physical_law_claim_allowed": self.universal_physical_law_claim_allowed,
            "folding_problem_solved": self.folding_problem_solved,
            "claim_rejection_reason": self.claim_rejection_reason,
            "rows": [row.to_dict() for row in self.rows],
            "controls": [control.to_dict() for control in self.controls],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "coordinate_truth_used_before_selection": self.coordinate_truth_used_before_selection,
            "native_truth_used_before_selection": self.native_truth_used_before_selection,
            "learned_prior_used_before_selection": self.learned_prior_used_before_selection,
            "msa_used_before_selection": self.msa_used_before_selection,
            "template_used_before_selection": self.template_used_before_selection,
            "raw_sequence_exposed": self.raw_sequence_exposed,
        }


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, float(value)))


def _rounded(value: float) -> float:
    return round(_clamp(value), 6)


def _score(value: float) -> float:
    return round(float(value), 6)


def _mean(values: Sequence[float]) -> float:
    return _rounded(mean(values)) if values else 0.0


def _window(sequence: str, index: int, radius: int = 4) -> str:
    left = max(0, index - radius - 1)
    right = min(len(sequence), index + radius)
    return sequence[left:right]


def _fraction(values: Iterable[bool]) -> float:
    collected = tuple(values)
    if not collected:
        return 0.0
    return sum(1 for value in collected if value) / len(collected)


def _hydrophobic_fraction(sequence: str, index: int) -> float:
    window = _window(sequence, index)
    return _fraction(aa in HYDROPHOBIC_AMINO_ACIDS for aa in window)


def _breaker_fraction(sequence: str, index: int) -> float:
    window = _window(sequence, index)
    return _fraction(aa in BREAKERS for aa in window)


def _charged(aa: str) -> bool:
    return aa in POSITIVE_CHARGED or aa in NEGATIVE_CHARGED


def _opposite_charge(left: str, right: str) -> bool:
    return (left in POSITIVE_CHARGED and right in NEGATIVE_CHARGED) or (
        right in POSITIVE_CHARGED and left in NEGATIVE_CHARGED
    )


def _same_charge(left: str, right: str) -> bool:
    return (left in POSITIVE_CHARGED and right in POSITIVE_CHARGED) or (
        left in NEGATIVE_CHARGED and right in NEGATIVE_CHARGED
    )


def _pair_energy_support(sequence: str, pair: ContactPair) -> tuple[float, float]:
    energy = contact_energy_kcal(sequence, pair)
    return energy, contact_energy_score(energy)


def _secondary_support(sequence: str, pair: ContactPair, ss: Sequence[str]) -> float:
    return float(secondary_structure_pair_score(pair, ss))


def _entropy_retention_score(
    sequence: str,
    pair: ContactPair,
    *,
    secondary_support: float,
    energy_support: float,
) -> float:
    left, right = pair
    length = max(1, len(sequence))
    separation = right - left
    loop_cost = _clamp((separation / length) ** 0.65)
    local_breaker_pressure = (_breaker_fraction(sequence, left) + _breaker_fraction(sequence, right)) / 2
    order_support = 0.60 * secondary_support + 0.40 * energy_support
    return _rounded(
        1.0
        - 0.72 * loop_cost
        - 0.18 * local_breaker_pressure
        + 0.34 * order_support
    )


def _environmental_context_score(sequence: str, pair: ContactPair) -> float:
    left_index, right_index = pair
    left = sequence[left_index - 1]
    right = sequence[right_index - 1]
    left_h = _hydrophobic_fraction(sequence, left_index)
    right_h = _hydrophobic_fraction(sequence, right_index)
    hydrophobic_context = (left_h + right_h) / 2
    hydrophobic_pair = left in HYDROPHOBIC_AMINO_ACIDS and right in HYDROPHOBIC_AMINO_ACIDS
    aromatic_anchor = (
        left in AROMATIC_AMINO_ACIDS
        and right in HYDROPHOBIC_AMINO_ACIDS
        or right in AROMATIC_AMINO_ACIDS
        and left in HYDROPHOBIC_AMINO_ACIDS
    )
    polar_pair = left in POLAR_UNCHARGED and right in POLAR_UNCHARGED

    score = 0.30
    if hydrophobic_pair:
        score += 0.48 * hydrophobic_context + 0.18
    if aromatic_anchor:
        score += 0.16 * max(hydrophobic_context, 0.45)
    if _opposite_charge(left, right):
        # Salt bridges are more plausible when not buried deeply in the coarse
        # hydrophobic-context proxy.
        score += 0.34 * (1.0 - hydrophobic_context) + 0.12
    if _same_charge(left, right):
        score -= 0.34 + 0.18 * hydrophobic_context
    if _charged(left) != _charged(right) and hydrophobic_context >= 0.58:
        score -= 0.16
    if polar_pair:
        score += 0.10 * (1.0 - hydrophobic_context)
    if left in BREAKERS or right in BREAKERS:
        score -= 0.10
    return _rounded(score)


def _neighbour_pairs(pair: ContactPair, length: int) -> tuple[ContactPair, ...]:
    left, right = pair
    neighbours: list[ContactPair] = []
    for d_left, d_right in ((-2, -2), (-1, -1), (1, 1), (2, 2), (-1, 1), (1, -1), (0, -1), (0, 1)):
        candidate = (left + d_left, right + d_right)
        if candidate[0] < 1 or candidate[1] > length:
            continue
        if candidate[1] - candidate[0] < MIN_SEQUENCE_SEPARATION:
            continue
        neighbours.append(candidate)
    return tuple(dict.fromkeys(neighbours))


def _cooperative_neighbour_score(
    sequence: str,
    pair: ContactPair,
    *,
    energy_cache: Mapping[ContactPair, float],
    ss: Sequence[str],
) -> float:
    neighbours = _neighbour_pairs(pair, len(sequence))
    if not neighbours:
        return 0.0
    local_scores: list[float] = []
    for neighbour in neighbours:
        energy_support = energy_cache.get(neighbour)
        if energy_support is None:
            _energy, energy_support = _pair_energy_support(sequence, neighbour)
        ss_support = _secondary_support(sequence, neighbour, ss)
        context_support = _environmental_context_score(sequence, neighbour)
        local_scores.append(0.48 * energy_support + 0.30 * ss_support + 0.22 * context_support)
    density = sum(1 for score in local_scores if score >= 0.55) / len(local_scores)
    return _rounded(0.72 * mean(local_scores) + 0.28 * density)


def _dynamic_state_scores(
    *,
    energy_support: float,
    entropy_retention: float,
    cooperative_support: float,
    context_support: float,
    secondary_support: float,
    sequence: str,
    pair: ContactPair,
) -> tuple[float, float, float, float, float]:
    left, right = pair
    length = max(1, len(sequence))
    separation = right - left
    normalized_separation = separation / length
    breaker_pressure = (_breaker_fraction(sequence, left) + _breaker_fraction(sequence, right)) / 2
    hydrophobic_context = (_hydrophobic_fraction(sequence, left) + _hydrophobic_fraction(sequence, right)) / 2

    compact_state = _rounded(
        0.38 * energy_support
        + 0.28 * context_support
        + 0.24 * cooperative_support
        + 0.10 * secondary_support
    )
    open_state = _rounded(
        0.34 * entropy_retention
        + 0.26 * secondary_support
        + 0.22 * energy_support
        + 0.18 * (1.0 - hydrophobic_context)
    )
    hinge_window = _clamp(1.0 - abs(normalized_separation - 0.34) / 0.34)
    hinge_state = _rounded(
        0.30 * entropy_retention
        + 0.25 * cooperative_support
        + 0.20 * hinge_window
        + 0.15 * breaker_pressure
        + 0.10 * context_support
    )
    values = (compact_state, open_state, hinge_state)
    dynamic_score = _rounded(0.62 * max(values) + 0.38 * mean(values))
    ambiguity = _rounded(max(values) - min(values))
    return compact_state, open_state, hinge_state, dynamic_score, ambiguity


def _candidate_pairs(row: RealCoordinateVisualRow, *, max_sequence_separation: int | None) -> tuple[ContactPair, ...]:
    length = row.sequence_length
    max_sep = max_sequence_separation or length
    max_sep = max(MIN_SEQUENCE_SEPARATION, min(max_sep, length - 1))
    stride = 1 if length <= 240 else 2
    pairs: list[ContactPair] = []
    for left in range(1, length + 1, stride):
        right_stop = min(length, left + max_sep)
        for right in range(left + MIN_SEQUENCE_SEPARATION, right_stop + 1):
            pairs.append((left, right))
    return tuple(pairs)


def _boundary_from_scores(scores: Mapping[ContactPair, float], *, sequence_length: int) -> FiveAxisBoundary:
    ordered = sorted(scores.items(), key=lambda item: (-item[1], item[0][0], item[0][1]))
    positive = tuple(item for item in ordered if item[1] > 0.0)
    if not positive:
        return FiveAxisBoundary(
            selected_count=0,
            boundary_score=0.0,
            largest_gap=0.0,
            positive_score_count=0,
            selection_reason="empty_positive_five_axis_distribution",
        )
    if len(positive) == 1:
        return FiveAxisBoundary(
            selected_count=1,
            boundary_score=_score(positive[0][1]),
            largest_gap=0.0,
            positive_score_count=1,
            selection_reason="single_positive_five_axis_pair",
        )

    values = [score for _pair, score in positive]
    gaps = [values[index] - values[index + 1] for index in range(len(values) - 1)]
    floor = max(8, int(sqrt(max(1, sequence_length))) + 2)
    cap = max(18, min(120, sequence_length // 2))
    allowed_indices = [index for index in range(len(gaps)) if floor <= index + 1 <= cap]
    if not allowed_indices:
        selected_count = min(cap, len(positive))
        largest_gap = 0.0
        reason = "bounded_top_score_cap_no_gap_window"
    else:
        largest_gap_index = max(allowed_indices, key=lambda index: gaps[index])
        selected_count = largest_gap_index + 1
        largest_gap = gaps[largest_gap_index]
        if largest_gap <= 0.004:
            selected_count = min(cap, len(positive))
            reason = "weak_gap_bounded_top_score_cap"
        else:
            reason = "row_local_largest_gap_after_sqrt_length_floor"
    selected_count = max(0, min(selected_count, len(positive)))
    return FiveAxisBoundary(
        selected_count=selected_count,
        boundary_score=_score(positive[selected_count - 1][1]) if selected_count else 0.0,
        largest_gap=_score(largest_gap),
        positive_score_count=len(positive),
        selection_reason=reason,
    )


def build_five_axis_contact_decisions(
    row: RealCoordinateVisualRow,
    *,
    max_sequence_separation: int | None = 160,
) -> tuple[FiveAxisContactDecision, ...]:
    """Build sequence-only five-axis contact decisions for one row.

    Native labels and coordinates are not read by this function. The returned
    decisions only contain sequence-derived scores and a row-local selected flag.
    """

    sequence = row.sequence
    pairs = _candidate_pairs(row, max_sequence_separation=max_sequence_separation)
    ss = predict_lightweight_secondary_structure(sequence)
    energy_cache: dict[ContactPair, float] = {}
    raw_decisions: list[dict[str, object]] = []
    scores: dict[ContactPair, float] = {}

    for pair in pairs:
        energy, energy_support = _pair_energy_support(sequence, pair)
        energy_cache[pair] = energy_support
        secondary_support = _secondary_support(sequence, pair, ss)
        entropy_support = _entropy_retention_score(
            sequence,
            pair,
            secondary_support=secondary_support,
            energy_support=energy_support,
        )
        context_support = _environmental_context_score(sequence, pair)
        cooperative_support = _cooperative_neighbour_score(
            sequence,
            pair,
            energy_cache=energy_cache,
            ss=ss,
        )
        compact, open_state, hinge, dynamic_score, ambiguity = _dynamic_state_scores(
            energy_support=energy_support,
            entropy_retention=entropy_support,
            cooperative_support=cooperative_support,
            context_support=context_support,
            secondary_support=secondary_support,
            sequence=sequence,
            pair=pair,
        )
        free_energy = _rounded(0.62 * energy_support + 0.38 * entropy_support)
        final_score = _rounded(
            0.25 * free_energy
            + 0.21 * cooperative_support
            + 0.20 * context_support
            + 0.18 * dynamic_score
            + 0.10 * secondary_support
            + 0.06 * (1.0 - ambiguity)
        )
        scores[pair] = final_score
        raw_decisions.append(
            {
                "pair": pair,
                "contact_energy_kcal": energy,
                "energy_support_score": energy_support,
                "entropy_retention_score": entropy_support,
                "cooperative_neighbour_score": cooperative_support,
                "environmental_context_score": context_support,
                "compact_state_score": compact,
                "open_state_score": open_state,
                "hinge_state_score": hinge,
                "dynamic_ensemble_score": dynamic_score,
                "dynamic_state_ambiguity": ambiguity,
                "free_energy_support_score": free_energy,
                "final_score": final_score,
            }
        )

    boundary = _boundary_from_scores(scores, sequence_length=row.sequence_length)
    selected_pairs = {
        pair
        for pair, _value in sorted(scores.items(), key=lambda item: (-item[1], item[0][0], item[0][1]))[
            : boundary.selected_count
        ]
    }
    output: list[FiveAxisContactDecision] = []
    for raw in raw_decisions:
        pair = raw["pair"]
        assert isinstance(pair, tuple)
        selected = pair in selected_pairs
        output.append(
            FiveAxisContactDecision(
                kind=FIVE_AXIS_DECISION_KIND,
                row_id=row.row_id,
                source_accession=row.source_accession,
                sequence_hash=row.sequence_sha256,
                sequence_length=row.sequence_length,
                i=pair[0],
                j=pair[1],
                sequence_separation=pair[1] - pair[0],
                contact_energy_kcal=float(raw["contact_energy_kcal"]),
                energy_support_score=float(raw["energy_support_score"]),
                entropy_retention_score=float(raw["entropy_retention_score"]),
                cooperative_neighbour_score=float(raw["cooperative_neighbour_score"]),
                environmental_context_score=float(raw["environmental_context_score"]),
                compact_state_score=float(raw["compact_state_score"]),
                open_state_score=float(raw["open_state_score"]),
                hinge_state_score=float(raw["hinge_state_score"]),
                dynamic_ensemble_score=float(raw["dynamic_ensemble_score"]),
                dynamic_state_ambiguity=float(raw["dynamic_state_ambiguity"]),
                free_energy_support_score=float(raw["free_energy_support_score"]),
                final_score=float(raw["final_score"]),
                selected=selected,
                selection_reason=boundary.selection_reason if selected else "below_five_axis_boundary",
            )
        )
    output.sort(key=lambda item: (-item.final_score, item.i, item.j))
    return tuple(output)


def selected_pairs_from_decisions(decisions: Sequence[FiveAxisContactDecision]) -> tuple[ContactPair, ...]:
    return normalized_contact_pairs(decision.pair() for decision in decisions if decision.selected)


def _stable_index(*parts: object, modulus: int) -> int:
    import hashlib

    encoded = ":".join(str(part) for part in parts).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest(), 16) % max(1, modulus)


def _pairs_by_separation(pairs: Iterable[ContactPair]) -> dict[int, tuple[ContactPair, ...]]:
    buckets: dict[int, list[ContactPair]] = defaultdict(list)
    for pair in sorted(set(pairs)):
        buckets[pair[1] - pair[0]].append(pair)
    return {sep: tuple(values) for sep, values in buckets.items()}


def matched_control_pairs(
    *,
    row: RealCoordinateVisualRow,
    selected_pairs: Sequence[ContactPair],
    candidate_pairs: Sequence[ContactPair],
    control_index: int,
) -> tuple[ContactPair, ...]:
    buckets = _pairs_by_separation(candidate_pairs)
    selected_set = set(selected_pairs)
    control: list[ContactPair] = []
    used: set[ContactPair] = set()
    for pair_index, selected_pair in enumerate(sorted(selected_set)):
        separation = selected_pair[1] - selected_pair[0]
        bucket = buckets.get(separation, ())
        if not bucket:
            continue
        start = _stable_index(
            row.row_id,
            "five_axis_control",
            control_index,
            pair_index,
            selected_pair[0],
            selected_pair[1],
            modulus=len(bucket),
        )
        chosen: ContactPair | None = None
        for offset in range(len(bucket)):
            candidate = bucket[(start + offset) % len(bucket)]
            if candidate in selected_set or candidate in used:
                continue
            chosen = candidate
            break
        if chosen is not None:
            used.add(chosen)
            control.append(chosen)
    return normalized_contact_pairs(control)


def _control_count(selected_pair_count: int) -> int:
    return max(2, min(6, int(sqrt(max(1, selected_pair_count)))))


def _selected_mean(decisions: Sequence[FiveAxisContactDecision], attribute: str) -> float:
    values = [float(getattr(decision, attribute)) for decision in decisions if decision.selected]
    return _mean(values)


def run_five_axis_row_challenge(
    row: RealCoordinateVisualRow,
    *,
    max_sequence_separation: int | None = 160,
) -> tuple[FiveAxisRowReport, tuple[FiveAxisControlResult, ...], tuple[FiveAxisContactDecision, ...]]:
    decisions = build_five_axis_contact_decisions(
        row,
        max_sequence_separation=max_sequence_separation,
    )
    selected_pairs = selected_pairs_from_decisions(decisions)
    candidate_pairs = tuple(decision.pair() for decision in decisions)
    boundary = _boundary_from_scores(
        {decision.pair(): decision.final_score for decision in decisions},
        sequence_length=row.sequence_length,
    )
    native_pairs = row.native_contact_pairs()
    metric = evaluate_contact_prediction(native_pairs=native_pairs, predicted_pairs=selected_pairs)
    selected_long = tuple(pair for pair in selected_pairs if pair[1] - pair[0] >= 24)

    controls: list[FiveAxisControlResult] = []
    for control_index in range(1, _control_count(len(selected_pairs)) + 1):
        control_pairs = matched_control_pairs(
            row=row,
            selected_pairs=selected_pairs,
            candidate_pairs=candidate_pairs,
            control_index=control_index,
        )
        control_metric = evaluate_contact_prediction(
            native_pairs=native_pairs,
            predicted_pairs=control_pairs,
        )
        controls.append(
            FiveAxisControlResult(
                control_id=f"{row.row_id}:five_axis_matched_control_{control_index:02d}",
                selected_contact_count=len(control_pairs),
                selected_contact_map_hash=contact_map_hash(control_pairs),
                native_contact_precision_after_audit=control_metric.native_contact_precision,
                native_contact_recall_after_audit=control_metric.native_contact_recall,
                long_range_contact_recall_after_audit=control_metric.long_range_contact_recall,
                contact_map_f1_after_audit=control_metric.contact_map_f1,
            )
        )
    best_control_f1 = max((control.contact_map_f1_after_audit for control in controls), default=0.0)
    best_control_lr = max((control.long_range_contact_recall_after_audit for control in controls), default=0.0)
    f1_margin = _score(metric.contact_map_f1 - best_control_f1)
    lr_margin = _score(metric.long_range_contact_recall - best_control_lr)
    row_allowed = (
        metric.native_contact_precision >= 0.60
        and metric.long_range_contact_recall >= 0.60
        and metric.contact_map_f1 >= 0.60
        and f1_margin >= 0.15
        and lr_margin >= 0.15
    )
    if row_allowed:
        rejection = "row_claim_survived_strict_five_axis_gate"
    elif metric.native_contact_precision < 0.60:
        rejection = "row_claim_rejected_precision_below_0_60"
    elif metric.long_range_contact_recall < 0.60:
        rejection = "row_claim_rejected_long_range_recall_below_0_60"
    elif metric.contact_map_f1 < 0.60:
        rejection = "row_claim_rejected_f1_below_0_60"
    elif f1_margin < 0.15 or lr_margin < 0.15:
        rejection = "row_claim_rejected_matched_control_margin_below_0_15"
    else:
        rejection = "row_claim_rejected_unknown_gate_failure"

    report = FiveAxisRowReport(
        row_id=row.row_id,
        source_accession=row.source_accession,
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        candidate_pair_count=len(candidate_pairs),
        selected_contact_count=len(selected_pairs),
        selected_long_range_contact_count=len(selected_long),
        selected_contact_map_hash=contact_map_hash(selected_pairs),
        boundary=boundary,
        metric_after_native_audit=metric,
        mean_selected_energy_support_score=_selected_mean(decisions, "energy_support_score"),
        mean_selected_entropy_retention_score=_selected_mean(decisions, "entropy_retention_score"),
        mean_selected_cooperative_neighbour_score=_selected_mean(decisions, "cooperative_neighbour_score"),
        mean_selected_environmental_context_score=_selected_mean(decisions, "environmental_context_score"),
        mean_selected_dynamic_ensemble_score=_selected_mean(decisions, "dynamic_ensemble_score"),
        mean_selected_free_energy_support_score=_selected_mean(decisions, "free_energy_support_score"),
        matched_control_count=len(controls),
        best_control_f1_after_audit=best_control_f1,
        best_control_long_range_recall_after_audit=best_control_lr,
        f1_margin_vs_best_control=f1_margin,
        long_range_recall_margin_vs_best_control=lr_margin,
        row_physical_law_claim_allowed=row_allowed,
        row_claim_rejection_reason=rejection,
    )
    return report, tuple(controls), decisions


def run_five_axis_challenge(
    rows: Sequence[RealCoordinateVisualRow],
    *,
    max_sequence_separation: int | None = 160,
) -> FiveAxisChallengePacket:
    row_reports: list[FiveAxisRowReport] = []
    controls: list[FiveAxisControlResult] = []
    decisions: list[FiveAxisContactDecision] = []
    for row in rows:
        row_report, row_controls, row_decisions = run_five_axis_row_challenge(
            row,
            max_sequence_separation=max_sequence_separation,
        )
        row_reports.append(row_report)
        controls.extend(row_controls)
        # Keep all selected decisions and the top rejected tail for audit without
        # writing every possible pair into the default report artifact.
        selected = [decision for decision in row_decisions if decision.selected]
        rejected = [decision for decision in row_decisions if not decision.selected][:20]
        decisions.extend(selected + rejected)

    precision_values = [row.metric_after_native_audit.native_contact_precision for row in row_reports]
    recall_values = [row.metric_after_native_audit.native_contact_recall for row in row_reports]
    long_recall_values = [row.metric_after_native_audit.long_range_contact_recall for row in row_reports]
    f1_values = [row.metric_after_native_audit.contact_map_f1 for row in row_reports]
    f1_margins = [row.f1_margin_vs_best_control for row in row_reports]
    lr_margins = [row.long_range_recall_margin_vs_best_control for row in row_reports]
    universal_allowed = (
        bool(row_reports)
        and all(row.row_physical_law_claim_allowed for row in row_reports)
        and min(precision_values) >= 0.60
        and min(long_recall_values) >= 0.60
        and _mean(f1_values) >= 0.60
        and _mean(f1_margins) >= 0.15
        and _mean(lr_margins) >= 0.15
    )
    if universal_allowed:
        claim_rejection = "five_axis_physical_law_claim_survived_strict_gate"
    else:
        failed_rows = [row.source_accession for row in row_reports if not row.row_physical_law_claim_allowed]
        claim_rejection = (
            "global_folding_claim_rejected_five_axis_gate_failed_for_rows:"
            + ",".join(failed_rows[:12])
        )
    return FiveAxisChallengePacket(
        kind=FIVE_AXIS_PHYSICS_KIND,
        row_count=len(row_reports),
        decision_rule=FIVE_AXIS_SCORING_RULE,
        claim_rule=FIVE_AXIS_CLAIM_RULE,
        entropy_axis_included=True,
        cooperativity_axis_included=True,
        dynamics_axis_included=True,
        context_axis_included=True,
        independent_physics_axis_included=True,
        mean_native_contact_precision_after_audit=_mean(precision_values),
        mean_native_contact_recall_after_audit=_mean(recall_values),
        mean_long_range_contact_recall_after_audit=_mean(long_recall_values),
        mean_contact_map_f1_after_audit=_mean(f1_values),
        min_row_precision_after_audit=_rounded(min(precision_values)) if precision_values else 0.0,
        min_row_long_range_recall_after_audit=_rounded(min(long_recall_values)) if long_recall_values else 0.0,
        mean_f1_margin_vs_best_control=_mean(f1_margins),
        mean_long_range_recall_margin_vs_best_control=_mean(lr_margins),
        universal_physical_law_claim_allowed=universal_allowed,
        folding_problem_solved=False if not universal_allowed else True,
        claim_rejection_reason=claim_rejection,
        rows=tuple(row_reports),
        controls=tuple(controls),
        decisions=tuple(decisions),
    )
