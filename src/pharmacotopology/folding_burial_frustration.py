from __future__ import annotations

from dataclasses import asdict, dataclass

from pharmacotopology.folding_nucleus_closure_search import NucleusClosureEvent
from pharmacotopology.folding_topology import (
    AROMATIC_AMINO_ACIDS,
    CHARGED_AMINO_ACIDS,
    DISORDER_PROMOTING_AMINO_ACIDS,
    HYDROPHOBIC_AMINO_ACIDS,
    POLAR_AMINO_ACIDS,
)


BURIAL_FRUSTRATION_KIND = "sequence_only_burial_frustration_proxy_v1"


@dataclass(frozen=True)
class BurialFrustrationPacket:
    hydrophobic_closure_fraction: float
    aromatic_anchor_fraction: float
    polar_or_charged_closure_fraction: float
    loop_disorder_fraction: float
    loop_breaker_fraction: float
    charge_pair_satisfaction: float
    loop_strain: float
    steric_clash_score: float
    burial_gain: float
    unsatisfied_polar_penalty: float
    future_frustration_score: float
    native_truth_used_before_physical_scoring: bool = False
    raw_sequence_exposed: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _fraction(sequence: str, alphabet: frozenset[str] | set[str]) -> float:
    if not sequence:
        return 0.0
    return round(sum(1 for residue in sequence if residue in alphabet) / len(sequence), 6)


def _charge_pair_satisfaction(left: str, right: str) -> float:
    residues = left + right
    negative = sum(1 for residue in residues if residue in "DE")
    positive = sum(1 for residue in residues if residue in "KRH")
    return _rounded(min(negative, positive) / max(1, negative + positive))


def burial_frustration_for_event(
    *,
    event: NucleusClosureEvent,
    sequence: str,
) -> BurialFrustrationPacket:
    left_segment = sequence[event.segment_a_start - 1 : event.segment_a_end]
    right_segment = sequence[event.segment_b_start - 1 : event.segment_b_end]
    loop = sequence[event.segment_a_end : event.segment_b_start - 1]
    closure = left_segment + right_segment

    hydrophobic_fraction = _fraction(closure, HYDROPHOBIC_AMINO_ACIDS)
    aromatic_fraction = _fraction(closure, AROMATIC_AMINO_ACIDS)
    polar_fraction = _fraction(closure, POLAR_AMINO_ACIDS | CHARGED_AMINO_ACIDS)
    loop_disorder_fraction = _fraction(loop, DISORDER_PROMOTING_AMINO_ACIDS)
    loop_breaker_fraction = _fraction(loop, set("PG"))
    charge_satisfaction = _charge_pair_satisfaction(left_segment, right_segment)
    loop_length = max(0, event.segment_b_start - event.segment_a_end - 1)
    normalized_loop = loop_length / max(1, event.sequence_length)

    loop_strain = _rounded(
        max(0.0, event.normalized_span - 0.35) * 0.55
        + max(0.0, 0.10 - normalized_loop) * 0.30
        + loop_breaker_fraction * 0.15
        - loop_disorder_fraction * 0.08
    )
    steric_clash_score = _rounded(
        max(0.0, hydrophobic_fraction - 0.55) * 0.25
        + max(0.0, event.normalized_span - 0.65) * 0.35
        + aromatic_fraction * 0.12
        + event.geometry_violation_cost * 0.28
    )
    burial_gain = _rounded(
        hydrophobic_fraction * 0.55
        + aromatic_fraction * 0.25
        + event.hydrophobic_burial_gain * 0.35
        + charge_satisfaction * 0.10
        - polar_fraction * 0.12
    )
    unsatisfied_polar_penalty = _rounded(
        max(0.0, polar_fraction - charge_satisfaction * 0.70) * 0.55
        + max(0.0, hydrophobic_fraction - 0.62) * polar_fraction * 0.25
    )
    future_frustration_score = _rounded(
        event.loop_entropy_cost * 0.30
        + event.isolation_penalty * 0.18
        + event.frustration_cost * 0.22
        + max(0.0, event.normalized_span - 0.55) * 0.18
        + max(0.0, unsatisfied_polar_penalty - 0.25) * 0.12
    )
    return BurialFrustrationPacket(
        hydrophobic_closure_fraction=hydrophobic_fraction,
        aromatic_anchor_fraction=aromatic_fraction,
        polar_or_charged_closure_fraction=polar_fraction,
        loop_disorder_fraction=loop_disorder_fraction,
        loop_breaker_fraction=loop_breaker_fraction,
        charge_pair_satisfaction=charge_satisfaction,
        loop_strain=loop_strain,
        steric_clash_score=steric_clash_score,
        burial_gain=burial_gain,
        unsatisfied_polar_penalty=unsatisfied_polar_penalty,
        future_frustration_score=future_frustration_score,
    )
