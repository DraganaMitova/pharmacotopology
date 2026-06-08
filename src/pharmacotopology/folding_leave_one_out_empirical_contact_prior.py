from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_native_contact_eval import ContactPair, contact_map_hash, normalized_contact_pairs
from pharmacotopology.folding_real_coordinate_visual_benchmark import RealCoordinateVisualRow

EMPIRICAL_CONTACT_PRIOR_KIND = "leave_one_out_empirical_sequence_contact_prior_v0"
EMPIRICAL_CONTACT_PRIOR_RULE = (
    "train_on_non_target_coordinate_rows_only; score_target_from_sequence_features; "
    "target_native_truth_used_only_after_selection_for_audit"
)
MIN_SEQUENCE_SEPARATION = 3


def _rounded_score(value: float) -> float:
    return round(float(value), 6)


def _aa_class(residue: str) -> str:
    if residue in "VILMFWYAC":
        return "hydrophobic"
    if residue in "KRH":
        return "positive"
    if residue in "DE":
        return "negative"
    if residue in "STQ":
        return "polar"
    if residue == "G":
        return "glycine"
    if residue == "P":
        return "proline"
    return "other"


def _sequence_separation_bin(separation: int, sequence_length: int) -> str:
    relative = separation / max(1, sequence_length)
    if separation <= 5:
        return "local_3_5"
    if separation <= 11:
        return "near_6_11"
    if separation <= 23:
        return "mid_12_23"
    if relative <= 0.25:
        return "long_quarter"
    if relative <= 0.50:
        return "long_half"
    return "very_long"


def _position_bin(i: int, j: int, sequence_length: int) -> str:
    center = (i + j) / (2.0 * max(1, sequence_length))
    if center < 0.25:
        return "n_terminal"
    if center < 0.50:
        return "n_mid"
    if center < 0.75:
        return "c_mid"
    return "c_terminal"


def empirical_sequence_features(row: RealCoordinateVisualRow, pair: ContactPair) -> tuple[tuple[str, object], ...]:
    i, j = pair
    left = row.sequence[i - 1]
    right = row.sequence[j - 1]
    separation = j - i
    sequence_length = row.sequence_length
    class_pair = tuple(sorted((_aa_class(left), _aa_class(right))))
    separation_bin = _sequence_separation_bin(separation, sequence_length)
    position_bin = _position_bin(i, j, sequence_length)
    return (
        ("sequence_separation_bin", separation_bin),
        ("residue_class_pair", class_pair),
        ("separation_x_class_pair", (separation_bin, class_pair)),
        ("sequence_position_bin", position_bin),
        ("position_x_separation", (position_bin, separation_bin)),
        ("near_hydrophobic_pair", separation <= 12 and left in "VILMFWYAC" and right in "VILMFWYAC"),
        ("hydrophobic_pair", left in "VILMFWYAC" and right in "VILMFWYAC"),
        ("salt_bridge_pair", (left in "KRH" and right in "DE") or (right in "KRH" and left in "DE")),
        ("same_charge_pair", (left in "KRH" and right in "KRH") or (left in "DE" and right in "DE")),
        ("glycine_or_proline_pair", left in "GP" or right in "GP"),
    )


def all_candidate_pairs(row: RealCoordinateVisualRow, *, minimum_sequence_separation: int = MIN_SEQUENCE_SEPARATION) -> tuple[ContactPair, ...]:
    return tuple(
        (i, j)
        for i in range(1, row.sequence_length + 1)
        for j in range(i + minimum_sequence_separation, row.sequence_length + 1)
    )


def _degree_limited_selection(
    scores: Mapping[ContactPair, float],
    *,
    budget: int,
    max_degree: int,
) -> tuple[ContactPair, ...]:
    degree: dict[int, int] = defaultdict(int)
    selected: list[ContactPair] = []
    # Pair order is part of the deterministic native-free tie-break. It is not
    # chosen from target native truth. Reverse pair order avoids a fixed N-term
    # bias when many coarse feature bins tie.
    for pair, _score in sorted(scores.items(), key=lambda item: (-item[1], -item[0][0], -item[0][1])):
        i, j = pair
        if degree[i] >= max_degree or degree[j] >= max_degree:
            continue
        selected.append(pair)
        degree[i] += 1
        degree[j] += 1
        if len(selected) >= budget:
            break
    return normalized_contact_pairs(selected)


@dataclass(frozen=True)
class EmpiricalContactPriorScore:
    row_id: str
    source_accession: str
    i: int
    j: int
    sequence_separation: int
    empirical_log_odds_score: float
    selected: bool
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    alphafold_used_before_selection: bool = False
    msa_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def pair(self) -> ContactPair:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class EmpiricalContactPriorPacket:
    kind: str
    row_id: str
    source_accession: str
    sequence_hash: str
    sequence_length: int
    training_row_count: int
    training_source_accessions: tuple[str, ...]
    training_contact_pair_count: int
    training_noncontact_pair_count: int
    candidate_pair_count: int
    selected_pair_count: int
    selected_long_range_pair_count: int
    contact_map_hash: str
    budget: int
    max_degree: int
    decision_rule: str
    selected_pairs: tuple[ContactPair, ...]
    score_rows: tuple[EmpiricalContactPriorScore, ...]
    target_native_truth_used_before_selection: bool = False
    target_coordinate_truth_used_before_selection: bool = False
    alphafold_used_before_selection: bool = False
    msa_used_before_selection: bool = False
    raw_sequence_persisted: bool = False

    def to_report_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "row_id": self.row_id,
            "source_accession": self.source_accession,
            "sequence_hash": self.sequence_hash,
            "sequence_length": self.sequence_length,
            "training_row_count": self.training_row_count,
            "training_source_accessions": list(self.training_source_accessions),
            "training_contact_pair_count": self.training_contact_pair_count,
            "training_noncontact_pair_count": self.training_noncontact_pair_count,
            "candidate_pair_count": self.candidate_pair_count,
            "selected_pair_count": self.selected_pair_count,
            "selected_long_range_pair_count": self.selected_long_range_pair_count,
            "contact_map_hash": self.contact_map_hash,
            "budget": self.budget,
            "max_degree": self.max_degree,
            "decision_rule": self.decision_rule,
            "target_native_truth_used_before_selection": self.target_native_truth_used_before_selection,
            "target_coordinate_truth_used_before_selection": self.target_coordinate_truth_used_before_selection,
            "alphafold_used_before_selection": self.alphafold_used_before_selection,
            "msa_used_before_selection": self.msa_used_before_selection,
            "raw_sequence_persisted": self.raw_sequence_persisted,
        }


def _contact_training_pairs(row: RealCoordinateVisualRow) -> set[ContactPair]:
    return set(row.native_contact_pairs())


def build_leave_one_out_empirical_contact_prior(
    *,
    target_row: RealCoordinateVisualRow,
    training_rows: Sequence[RealCoordinateVisualRow],
    budget_fraction: float = 2.06,
    max_degree: int = 6,
    feature_weight: float = 0.33,
    minimum_sequence_separation: int = MIN_SEQUENCE_SEPARATION,
) -> EmpiricalContactPriorPacket:
    if any(row.row_id == target_row.row_id for row in training_rows):
        raise ValueError("target row must not be included in empirical prior training rows")
    if any(row.source_accession == target_row.source_accession for row in training_rows):
        raise ValueError("target source accession must not be included in empirical prior training rows")
    if not training_rows:
        raise ValueError("at least one non-target training row is required")

    feature_counts: dict[tuple[str, object], list[int]] = defaultdict(lambda: [0, 0])
    base_counts = [0, 0]
    for row in training_rows:
        contacts = _contact_training_pairs(row)
        for pair in all_candidate_pairs(row, minimum_sequence_separation=minimum_sequence_separation):
            is_contact = pair in contacts
            label_index = 1 if is_contact else 0
            base_counts[label_index] += 1
            for feature in empirical_sequence_features(row, pair):
                feature_counts[feature][label_index] += 1

    base_rate = (base_counts[1] + 1.0) / (base_counts[0] + base_counts[1] + 2.0)
    base_log_odds = math.log(base_rate / (1.0 - base_rate))
    candidate_pairs = all_candidate_pairs(target_row, minimum_sequence_separation=minimum_sequence_separation)
    raw_scores: dict[ContactPair, float] = {}
    for pair in candidate_pairs:
        score = base_log_odds
        for feature in empirical_sequence_features(target_row, pair):
            negative_count, positive_count = feature_counts[feature]
            probability = (positive_count + 1.0) / (positive_count + negative_count + 2.0)
            score += feature_weight * (math.log(probability / (1.0 - probability)) - base_log_odds)
        separation = pair[1] - pair[0]
        if separation > 0.55 * target_row.sequence_length:
            score -= 0.55
        if separation <= 5:
            score += 0.20
        raw_scores[pair] = score

    budget = max(1, int(round(target_row.sequence_length * budget_fraction)))
    selected_pairs = _degree_limited_selection(raw_scores, budget=budget, max_degree=max_degree)
    selected_set = set(selected_pairs)
    score_rows = tuple(
        EmpiricalContactPriorScore(
            row_id=target_row.row_id,
            source_accession=target_row.source_accession,
            i=pair[0],
            j=pair[1],
            sequence_separation=pair[1] - pair[0],
            empirical_log_odds_score=_rounded_score(raw_scores[pair]),
            selected=pair in selected_set,
        )
        for pair in sorted(raw_scores, key=lambda pair: (-raw_scores[pair], -pair[0], -pair[1]))[: max(600, budget + 96)]
    )
    return EmpiricalContactPriorPacket(
        kind=EMPIRICAL_CONTACT_PRIOR_KIND,
        row_id=target_row.row_id,
        source_accession=target_row.source_accession,
        sequence_hash=target_row.sequence_sha256,
        sequence_length=target_row.sequence_length,
        training_row_count=len(training_rows),
        training_source_accessions=tuple(row.source_accession for row in training_rows),
        training_contact_pair_count=base_counts[1],
        training_noncontact_pair_count=base_counts[0],
        candidate_pair_count=len(candidate_pairs),
        selected_pair_count=len(selected_pairs),
        selected_long_range_pair_count=sum(1 for pair in selected_pairs if pair[1] - pair[0] >= 24),
        contact_map_hash=contact_map_hash(selected_pairs),
        budget=budget,
        max_degree=max_degree,
        decision_rule=EMPIRICAL_CONTACT_PRIOR_RULE,
        selected_pairs=selected_pairs,
        score_rows=score_rows,
    )
