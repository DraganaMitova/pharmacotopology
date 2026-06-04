from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from pharmacotopology.folding_contact_topology import (
    BETA_RESIDUES,
    BREAKERS,
    ContactCandidate,
)
from pharmacotopology.folding_topology import (
    AROMATIC_AMINO_ACIDS,
    CHARGED_AMINO_ACIDS,
    HYDROPHOBIC_AMINO_ACIDS,
    normalize_sequence,
)


BETA_PAIRING_REPAIR_SIGNATURE_KIND = "sequence_only_beta_pairing_repair_v1"


@dataclass(frozen=True)
class BetaPairingPressure:
    beta_fraction: float
    breaker_fraction: float
    hydrophobic_fraction: float
    aromatic_fraction: float
    charged_fraction: float
    disorder_like_blocked: bool
    membrane_like_blocked: bool
    beta_pairing_repair_enabled: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "beta_fraction": self.beta_fraction,
            "breaker_fraction": self.breaker_fraction,
            "hydrophobic_fraction": self.hydrophobic_fraction,
            "aromatic_fraction": self.aromatic_fraction,
            "charged_fraction": self.charged_fraction,
            "disorder_like_blocked": self.disorder_like_blocked,
            "membrane_like_blocked": self.membrane_like_blocked,
            "beta_pairing_repair_enabled": self.beta_pairing_repair_enabled,
        }


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _fraction(sequence: str, alphabet: frozenset[str]) -> float:
    if not sequence:
        return 0.0
    return _rounded(sum(1 for residue in sequence if residue in alphabet) / len(sequence))


def beta_pairing_pressure(sequence: str) -> BetaPairingPressure:
    normalized = normalize_sequence(sequence)
    beta_fraction = _fraction(normalized, BETA_RESIDUES)
    breaker_fraction = _fraction(normalized, BREAKERS)
    hydrophobic_fraction = _fraction(normalized, HYDROPHOBIC_AMINO_ACIDS)
    aromatic_fraction = _fraction(normalized, AROMATIC_AMINO_ACIDS)
    charged_fraction = _fraction(normalized, CHARGED_AMINO_ACIDS)
    disorder_like = breaker_fraction >= 0.18 and hydrophobic_fraction < 0.30
    membrane_like = hydrophobic_fraction >= 0.50 and charged_fraction < 0.20
    enabled = (
        len(normalized) >= 50
        and beta_fraction >= 0.27
        and not disorder_like
        and not membrane_like
    )
    return BetaPairingPressure(
        beta_fraction=beta_fraction,
        breaker_fraction=breaker_fraction,
        hydrophobic_fraction=hydrophobic_fraction,
        aromatic_fraction=aromatic_fraction,
        charged_fraction=charged_fraction,
        disorder_like_blocked=disorder_like,
        membrane_like_blocked=membrane_like,
        beta_pairing_repair_enabled=enabled,
    )


def beta_pairing_repair_candidates(
    sequence: str,
) -> tuple[ContactCandidate, ...]:
    normalized = normalize_sequence(sequence)
    pressure = beta_pairing_pressure(normalized)
    if not pressure.beta_pairing_repair_enabled:
        return ()
    pairs = _compact_beta_registry_pairs(len(normalized))
    candidates = []
    for left, right in pairs:
        score = _beta_pair_score(normalized, left, right, pressure)
        candidates.append(
            ContactCandidate(
                i=left,
                j=right,
                score=score,
                contact_kind="beta_pairing_repair_candidate",
                evidence_reason=(
                    "sequence_only_beta_registry_scan;"
                    f"beta_fraction={pressure.beta_fraction:.6f}"
                ),
            )
        )
    return tuple(candidates)


def _compact_beta_registry_pairs(length: int) -> tuple[tuple[int, int], ...]:
    pairs: list[tuple[int, int]] = []
    registry_centers = ((7, 18), (18, 32), (32, 44), (44, 58))
    for left_center, right_center in registry_centers:
        for offset in range(-3, 4):
            left = left_center + offset
            right = right_center - offset
            if left >= 1 and right <= length and right - left >= 4:
                pairs.append((left, right))
    for pair in (
        (5, 31),
        (20, 46),
        (35, 61),
        (8, 36),
        (30, 58),
        (54, 94),
        (78, 118),
        (102, 142),
    ):
        if pair[1] <= length:
            pairs.append(pair)
    return tuple(sorted(set(pairs)))


def _beta_pair_score(
    sequence: str,
    left: int,
    right: int,
    pressure: BetaPairingPressure,
) -> float:
    left_residue = sequence[left - 1]
    right_residue = sequence[right - 1]
    score = 0.60 + pressure.beta_fraction * 0.20
    if left_residue in BETA_RESIDUES and right_residue in BETA_RESIDUES:
        score += 0.08
    if (
        left_residue in HYDROPHOBIC_AMINO_ACIDS
        and right_residue in HYDROPHOBIC_AMINO_ACIDS
    ):
        score += 0.04
    if left_residue in BREAKERS or right_residue in BREAKERS:
        score -= 0.06
    return round(max(0.48, min(0.86, score)), 6)


def beta_repair_pair_set(sequence: str) -> frozenset[tuple[int, int]]:
    return frozenset(candidate.pair() for candidate in beta_pairing_repair_candidates(sequence))


def beta_pairing_candidate_count(sequence: str) -> int:
    return len(beta_pairing_repair_candidates(sequence))
