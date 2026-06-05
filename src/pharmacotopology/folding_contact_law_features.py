from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Mapping, Sequence

from pharmacotopology.folding_contact_topology import (
    BETA_RESIDUES,
    BREAKERS,
    HELIX_RESIDUES,
    _candidate_score,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
)
from pharmacotopology.folding_topology import (
    AROMATIC_AMINO_ACIDS,
    CHARGED_AMINO_ACIDS,
    HYDROPHOBIC_AMINO_ACIDS,
)


CONTACT_LAW_FEATURE_KIND = "sequence_only_contact_law_pair_features_v1"
CONTACT_LAW_FEATURE_BOUNDARY = "sequence_only_features_before_native_label_attach"


@dataclass(frozen=True)
class ContactLawFeatureRow:
    row_id: str
    source_accession: str
    sequence_hash: str
    sequence_length: int
    i: int
    j: int
    sequence_separation: int
    normalized_separation: float
    local_i_to_i4_support: float
    helix_window_support: float
    beta_window_support: float
    hydrophobic_pair_support: float
    aromatic_anchor_support: float
    opposite_charge_support: float
    same_charge_penalty: float
    breaker_penalty: float
    loop_entropy_cost: float
    cluster_neighbor_support: float
    parallel_contact_support: float
    isolation_penalty: float
    native_contact: bool
    current_scalar_score: float
    pair_only_score: float
    pair_plus_cluster_score: float
    pair_plus_entropy_score: float
    pair_plus_cluster_plus_entropy_score: float

    def pair(self) -> tuple[int, int]:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def rounded_score(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _window_fraction(sequence: str, index: int, residues: set[str] | frozenset[str]) -> float:
    left = max(0, index - 3)
    right = min(len(sequence), index + 2)
    window = sequence[left:right]
    if not window:
        return 0.0
    return sum(1 for aa in window if aa in residues) / len(window)


def _opposite_charge(left: str, right: str) -> bool:
    return (left in "KRH" and right in "DE") or (right in "KRH" and left in "DE")


def _same_charge(left: str, right: str) -> bool:
    return (left in "KRH" and right in "KRH") or (left in "DE" and right in "DE")


def _pair_chemistry_support(sequence: str, i: int, j: int) -> float:
    left = sequence[i - 1]
    right = sequence[j - 1]
    score = 0.0
    if left in HYDROPHOBIC_AMINO_ACIDS and right in HYDROPHOBIC_AMINO_ACIDS:
        score += 0.34
    if (
        left in AROMATIC_AMINO_ACIDS
        and right in HYDROPHOBIC_AMINO_ACIDS
        or right in AROMATIC_AMINO_ACIDS
        and left in HYDROPHOBIC_AMINO_ACIDS
    ):
        score += 0.14
    if _opposite_charge(left, right):
        score += 0.18
    if _same_charge(left, right):
        score -= 0.14
    if left in BREAKERS or right in BREAKERS:
        score -= 0.18
    return rounded_score(score)


def _cluster_neighbor_support(sequence: str, i: int, j: int) -> float:
    supports: list[float] = []
    for left_delta, right_delta in ((-1, -1), (-1, 1), (1, -1), (1, 1), (0, -1), (0, 1)):
        left = i + left_delta
        right = j + right_delta
        if left < 1 or right > len(sequence) or right - left < 3:
            continue
        supports.append(_pair_chemistry_support(sequence, left, right))
    if not supports:
        return 0.0
    return rounded_score(sum(supports) / len(supports))


def _feature_scores(
    *,
    current_scalar_score: float,
    local_i_to_i4_support: float,
    helix_window_support: float,
    beta_window_support: float,
    hydrophobic_pair_support: float,
    aromatic_anchor_support: float,
    opposite_charge_support: float,
    same_charge_penalty: float,
    breaker_penalty: float,
    loop_entropy_cost: float,
    cluster_neighbor_support: float,
    parallel_contact_support: float,
    isolation_penalty: float,
) -> dict[str, float]:
    pair_only = rounded_score(
        0.06
        + 0.28 * local_i_to_i4_support
        + 0.20 * helix_window_support
        + 0.18 * beta_window_support
        + 0.16 * hydrophobic_pair_support
        + 0.08 * aromatic_anchor_support
        + 0.08 * opposite_charge_support
        - 0.12 * same_charge_penalty
        - 0.14 * breaker_penalty
    )
    pair_plus_cluster = rounded_score(
        pair_only
        + 0.22 * cluster_neighbor_support
        + 0.10 * parallel_contact_support
        - 0.08 * isolation_penalty
    )
    pair_plus_entropy = rounded_score(
        pair_only
        - 0.20 * loop_entropy_cost
        + 0.10 * local_i_to_i4_support
    )
    cluster_entropy = rounded_score(
        pair_only
        + 0.26 * cluster_neighbor_support
        + 0.12 * parallel_contact_support
        - 0.18 * loop_entropy_cost
        - 0.08 * isolation_penalty
    )
    return {
        "current_scalar_score": round(current_scalar_score, 6),
        "pair_only_score": pair_only,
        "pair_plus_cluster_score": pair_plus_cluster,
        "pair_plus_entropy_score": pair_plus_entropy,
        "pair_plus_cluster_plus_entropy_score": cluster_entropy,
    }


def contact_law_feature_rows_for_row(
    row: RealCoordinateVisualRow,
) -> tuple[ContactLawFeatureRow, ...]:
    sequence = row.sequence
    native_pairs = set(row.native_contact_pairs())
    output: list[ContactLawFeatureRow] = []
    for i in range(1, len(sequence) + 1):
        for j in range(i + 3, len(sequence) + 1):
            left = sequence[i - 1]
            right = sequence[j - 1]
            separation = j - i
            local_support = 1.0 if separation == 4 and left not in BREAKERS and right not in BREAKERS else 0.0
            helix_support = round(
                (
                    _window_fraction(sequence, i, HELIX_RESIDUES)
                    + _window_fraction(sequence, j, HELIX_RESIDUES)
                )
                / 2,
                6,
            )
            beta_support = round(
                (
                    _window_fraction(sequence, i, BETA_RESIDUES)
                    + _window_fraction(sequence, j, BETA_RESIDUES)
                )
                / 2,
                6,
            )
            hydrophobic_support = (
                1.0
                if left in HYDROPHOBIC_AMINO_ACIDS
                and right in HYDROPHOBIC_AMINO_ACIDS
                else 0.0
            )
            aromatic_support = (
                1.0
                if (
                    left in AROMATIC_AMINO_ACIDS
                    and right in HYDROPHOBIC_AMINO_ACIDS
                )
                or (
                    right in AROMATIC_AMINO_ACIDS
                    and left in HYDROPHOBIC_AMINO_ACIDS
                )
                else 0.0
            )
            opposite_support = 1.0 if _opposite_charge(left, right) else 0.0
            same_penalty = 1.0 if _same_charge(left, right) else 0.0
            breaker_penalty = 1.0 if left in BREAKERS or right in BREAKERS else 0.0
            normalized_separation = round(separation / max(len(sequence), 1), 6)
            loop_entropy = rounded_score(normalized_separation ** 0.65)
            cluster_support = _cluster_neighbor_support(sequence, i, j)
            parallel_support = (
                1.0
                if separation >= 8
                and (i + j) % 2 == 0
                and beta_support >= 0.45
                else 0.0
            )
            isolation_penalty = rounded_score(1.0 - cluster_support)
            current_scalar, _, _ = _candidate_score(sequence, i, j)
            scores = _feature_scores(
                current_scalar_score=current_scalar,
                local_i_to_i4_support=local_support,
                helix_window_support=helix_support,
                beta_window_support=beta_support,
                hydrophobic_pair_support=hydrophobic_support,
                aromatic_anchor_support=aromatic_support,
                opposite_charge_support=opposite_support,
                same_charge_penalty=same_penalty,
                breaker_penalty=breaker_penalty,
                loop_entropy_cost=loop_entropy,
                cluster_neighbor_support=cluster_support,
                parallel_contact_support=parallel_support,
                isolation_penalty=isolation_penalty,
            )
            output.append(
                ContactLawFeatureRow(
                    row_id=row.row_id,
                    source_accession=row.source_accession,
                    sequence_hash=row.sequence_sha256,
                    sequence_length=row.sequence_length,
                    i=i,
                    j=j,
                    sequence_separation=separation,
                    normalized_separation=normalized_separation,
                    local_i_to_i4_support=local_support,
                    helix_window_support=helix_support,
                    beta_window_support=beta_support,
                    hydrophobic_pair_support=hydrophobic_support,
                    aromatic_anchor_support=aromatic_support,
                    opposite_charge_support=opposite_support,
                    same_charge_penalty=same_penalty,
                    breaker_penalty=breaker_penalty,
                    loop_entropy_cost=loop_entropy,
                    cluster_neighbor_support=cluster_support,
                    parallel_contact_support=parallel_support,
                    isolation_penalty=isolation_penalty,
                    native_contact=(i, j) in native_pairs,
                    **scores,
                )
            )
    return tuple(output)


def contact_law_feature_rows(
    rows: Sequence[RealCoordinateVisualRow],
) -> tuple[ContactLawFeatureRow, ...]:
    output: list[ContactLawFeatureRow] = []
    for row in rows:
        output.extend(contact_law_feature_rows_for_row(row))
    return tuple(output)


def feature_rows_by_row_id(
    rows: Sequence[ContactLawFeatureRow],
) -> dict[str, tuple[ContactLawFeatureRow, ...]]:
    grouped: dict[str, list[ContactLawFeatureRow]] = {}
    for row in rows:
        grouped.setdefault(row.row_id, []).append(row)
    return {
        row_id: tuple(row_values)
        for row_id, row_values in sorted(grouped.items())
    }


def feature_rows_to_dicts(
    rows: Sequence[ContactLawFeatureRow],
) -> list[dict[str, object]]:
    return [row.to_dict() for row in rows]


def score_for_model(row: ContactLawFeatureRow, model_id: str) -> float:
    value = getattr(row, model_id)
    if not isinstance(value, (int, float)):
        raise ValueError(f"unknown contact law score model: {model_id}")
    return float(value)


def native_pairs_from_feature_rows(
    rows: Sequence[ContactLawFeatureRow],
) -> tuple[tuple[int, int], ...]:
    return tuple((row.i, row.j) for row in rows if row.native_contact)


def predicted_pairs_from_threshold(
    rows: Sequence[ContactLawFeatureRow],
    *,
    model_id: str,
    threshold: float,
) -> tuple[tuple[int, int], ...]:
    return tuple(
        (row.i, row.j)
        for row in rows
        if score_for_model(row, model_id) >= threshold
    )


def feature_generation_certificate(
    rows: Sequence[ContactLawFeatureRow],
) -> Mapping[str, object]:
    return {
        "feature_kind": CONTACT_LAW_FEATURE_KIND,
        "feature_boundary": CONTACT_LAW_FEATURE_BOUNDARY,
        "pair_feature_row_count": len(rows),
        "native_truth_used_before_feature_generation": False,
        "native_label_attached_after_feature_generation": True,
        "raw_sequence_exposed": False,
    }
