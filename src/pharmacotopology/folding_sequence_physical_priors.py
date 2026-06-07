from __future__ import annotations

"""MSA-independent sequence/physics priors for contact scoring.

The functions in this module deliberately avoid native coordinates, predicted
structures, AlphaFold outputs, and MSA-derived coupling evidence.  They use only
coarse residue chemistry, a lightweight sequence-window secondary-structure
prior, and the current contact graph degree budget.
"""

from collections import defaultdict
from dataclasses import asdict, dataclass
from math import exp
from statistics import mean
from typing import Mapping, Sequence

from pharmacotopology.folding_contact_topology import (
    BETA_RESIDUES,
    BREAKERS,
    HELIX_RESIDUES,
)
from pharmacotopology.folding_native_contact_eval import ContactPair
from pharmacotopology.folding_real_coordinate_visual_benchmark import RealCoordinateVisualRow
from pharmacotopology.folding_topology import (
    AROMATIC_AMINO_ACIDS,
    HYDROPHOBIC_AMINO_ACIDS,
)


SEQUENCE_PHYSICAL_PRIOR_KIND = "sequence_only_energy_secondary_structure_degree_prior_v0"
SECONDARY_STRUCTURE_PRIOR_KIND = "lightweight_sequence_window_secondary_structure_prior_v0"
CONTACT_ENERGY_PRIOR_KIND = "sequence_pair_contact_energy_prior_v0"
CONTACT_DEGREE_PRIOR_KIND = "current_graph_contact_degree_consistency_prior_v0"

POSITIVE_CHARGED = frozenset("KRH")
NEGATIVE_CHARGED = frozenset("DE")
POLAR_UNCHARGED = frozenset("STNQ")


@dataclass(frozen=True)
class SequenceContactPhysicalPrior:
    row_id: str
    source_accession: str
    sequence_hash: str
    sequence_length: int
    i: int
    j: int
    sequence_separation: int
    contact_energy_kcal: float
    contact_energy_score: float
    secondary_structure_i: str
    secondary_structure_j: str
    secondary_structure_score: float
    degree_i_after_contact: int
    degree_j_after_contact: int
    degree_consistency_score: float
    physical_prior_score: float
    prior_kind: str = SEQUENCE_PHYSICAL_PRIOR_KIND
    contact_energy_prior_kind: str = CONTACT_ENERGY_PRIOR_KIND
    secondary_structure_prior_kind: str = SECONDARY_STRUCTURE_PRIOR_KIND
    contact_degree_prior_kind: str = CONTACT_DEGREE_PRIOR_KIND
    coordinate_truth_used_before_selection: bool = False
    native_truth_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def pair(self) -> ContactPair:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, float(value))), 6)


def _score(value: float) -> float:
    return round(float(value), 6)


def predict_lightweight_secondary_structure(sequence: str) -> tuple[str, ...]:
    """Return H/E/C calls from sequence windows only.

    This is intentionally lightweight and deterministic.  It is not a PSIPRED
    replacement; it is a native-free prior that separates obvious helix-like,
    strand-like, and breaker/coil-like local chemistry.
    """

    normalized = sequence.upper()
    calls: list[str] = []
    for index in range(1, len(normalized) + 1):
        window = normalized[max(0, index - 4) : min(len(normalized), index + 3)]
        if not window:
            calls.append("C")
            continue
        helix_pressure = mean(
            1.00 if aa in HELIX_RESIDUES else 0.45 if aa not in BREAKERS else 0.05
            for aa in window
        )
        beta_pressure = mean(
            1.00 if aa in BETA_RESIDUES else 0.40 if aa not in BREAKERS else 0.10
            for aa in window
        )
        breaker_fraction = sum(1 for aa in window if aa in BREAKERS) / len(window)
        if breaker_fraction >= 0.22:
            calls.append("C")
        elif helix_pressure >= beta_pressure + 0.12:
            calls.append("H")
        elif beta_pressure >= helix_pressure + 0.05:
            calls.append("E")
        else:
            calls.append("C")
    return tuple(calls)


def contact_energy_kcal(sequence: str, pair: ContactPair) -> float:
    """Coarse pair-contact energy from sequence chemistry only.

    Negative values are favorable.  The constants are deliberately coarse; they
    encode directionality for hydrophobic burial, salt bridges, charge clashes,
    polar compatibility, breaker penalties, and long-loop entropy.
    """

    left_index, right_index = pair
    left = sequence[left_index - 1].upper()
    right = sequence[right_index - 1].upper()
    separation = right_index - left_index
    energy = 0.0

    if left in HYDROPHOBIC_AMINO_ACIDS and right in HYDROPHOBIC_AMINO_ACIDS:
        energy -= 2.5
    if (
        left in AROMATIC_AMINO_ACIDS
        and right in HYDROPHOBIC_AMINO_ACIDS
        or right in AROMATIC_AMINO_ACIDS
        and left in HYDROPHOBIC_AMINO_ACIDS
    ):
        energy -= 0.8
    if (left in POSITIVE_CHARGED and right in NEGATIVE_CHARGED) or (
        right in POSITIVE_CHARGED and left in NEGATIVE_CHARGED
    ):
        energy -= 4.0
    if (left in POSITIVE_CHARGED and right in POSITIVE_CHARGED) or (
        left in NEGATIVE_CHARGED and right in NEGATIVE_CHARGED
    ):
        energy += 4.0
    if left in POLAR_UNCHARGED and right in POLAR_UNCHARGED:
        energy -= 0.3
    if left in BREAKERS or right in BREAKERS:
        energy += 0.9

    length = max(1, len(sequence))
    if separation > max(40, length // 3):
        energy += 0.35
    if separation > max(80, length // 2):
        energy += 0.35
    return _score(energy)


def contact_energy_score(energy_kcal: float) -> float:
    """Map lower/favorable contact energy onto a bounded support score."""

    return _rounded(1.0 / (1.0 + exp((float(energy_kcal) + 0.25) / 1.6)))


def secondary_structure_pair_score(
    pair: ContactPair,
    secondary_structure: Sequence[str],
) -> float:
    left_index, right_index = pair
    left = secondary_structure[left_index - 1]
    right = secondary_structure[right_index - 1]
    separation = right_index - left_index

    if separation in (3, 4, 5):
        return 0.92 if left == "H" and right == "H" else 0.48
    if left == "E" and right == "E":
        return 0.90 if separation >= 8 else 0.55
    if left == "H" and right == "H":
        # Long H/H contacts can be helix-packing contacts, but they should not
        # be over-promoted by a sequence-only prior.
        return 0.42 if separation >= 8 else 0.55
    if left == "C" and right == "C":
        return 0.38
    if "C" in (left, right):
        return 0.48
    return 0.58


def _adjacency(pairs: Sequence[ContactPair] | set[ContactPair]) -> dict[int, set[int]]:
    adjacency: dict[int, set[int]] = defaultdict(set)
    for left, right in pairs:
        adjacency[left].add(right)
        adjacency[right].add(left)
    return adjacency


def _single_degree_consistency(degree: int) -> float:
    if 2 <= degree <= 5:
        return 1.00
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


def degree_consistency_for_pair(
    pair: ContactPair,
    current_pairs: Sequence[ContactPair] | set[ContactPair],
) -> tuple[int, int, float]:
    adjacency = _adjacency(current_pairs)
    left, right = pair
    left_degree = len(adjacency[left]) + (0 if right in adjacency[left] else 1)
    right_degree = len(adjacency[right]) + (0 if left in adjacency[right] else 1)
    score = mean(
        (
            _single_degree_consistency(left_degree),
            _single_degree_consistency(right_degree),
        )
    )
    return left_degree, right_degree, _rounded(score)


def build_sequence_physical_prior_scores(
    *,
    row: RealCoordinateVisualRow,
    candidate_pairs: Sequence[ContactPair] | set[ContactPair],
    current_pairs: Sequence[ContactPair] | set[ContactPair] = (),
) -> dict[ContactPair, SequenceContactPhysicalPrior]:
    secondary_structure = predict_lightweight_secondary_structure(row.sequence)
    output: dict[ContactPair, SequenceContactPhysicalPrior] = {}
    for pair in candidate_pairs:
        if pair[0] < 1 or pair[1] > row.sequence_length or pair[1] <= pair[0]:
            continue
        energy = contact_energy_kcal(row.sequence, pair)
        energy_support = contact_energy_score(energy)
        ss_support = secondary_structure_pair_score(pair, secondary_structure)
        degree_i, degree_j, degree_support = degree_consistency_for_pair(pair, current_pairs)
        physical_score = _rounded(
            0.45 * energy_support
            + 0.30 * ss_support
            + 0.25 * degree_support
        )
        output[pair] = SequenceContactPhysicalPrior(
            row_id=row.row_id,
            source_accession=row.source_accession,
            sequence_hash=row.sequence_sha256,
            sequence_length=row.sequence_length,
            i=pair[0],
            j=pair[1],
            sequence_separation=pair[1] - pair[0],
            contact_energy_kcal=energy,
            contact_energy_score=energy_support,
            secondary_structure_i=secondary_structure[pair[0] - 1],
            secondary_structure_j=secondary_structure[pair[1] - 1],
            secondary_structure_score=ss_support,
            degree_i_after_contact=degree_i,
            degree_j_after_contact=degree_j,
            degree_consistency_score=degree_support,
            physical_prior_score=physical_score,
        )
    return output


def summarize_physical_priors(
    priors: Mapping[ContactPair, SequenceContactPhysicalPrior],
    pairs: Sequence[ContactPair] | set[ContactPair],
) -> dict[str, float]:
    selected = [priors[pair] for pair in pairs if pair in priors]
    if not selected:
        return {
            "mean_contact_energy_score": 0.0,
            "mean_secondary_structure_score": 0.0,
            "mean_degree_consistency_score": 0.0,
            "mean_physical_prior_score": 0.0,
        }
    return {
        "mean_contact_energy_score": _rounded(mean(item.contact_energy_score for item in selected)),
        "mean_secondary_structure_score": _rounded(mean(item.secondary_structure_score for item in selected)),
        "mean_degree_consistency_score": _rounded(mean(item.degree_consistency_score for item in selected)),
        "mean_physical_prior_score": _rounded(mean(item.physical_prior_score for item in selected)),
    }
