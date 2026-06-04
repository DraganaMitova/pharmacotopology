from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from pharmacotopology.folding_beta_pairing_repair import (
    BETA_PAIRING_REPAIR_SIGNATURE_KIND,
    beta_pairing_pressure,
    beta_pairing_repair_candidates,
)
from pharmacotopology.folding_contact_topology import (
    CONTACT_TOPOLOGY_SIGNATURE_KIND,
    HELIX_RESIDUES,
    ContactCandidate,
    ContactTopologyPrediction,
)
from pharmacotopology.folding_native_contact_eval import contact_map_hash
from pharmacotopology.folding_topology import (
    CHARGED_AMINO_ACIDS,
    HYDROPHOBIC_AMINO_ACIDS,
    normalize_sequence,
)


LONG_RANGE_CONTACT_REPAIR_SIGNATURE_KIND = "sequence_only_long_range_contact_repair_v1"
REPAIRED_CONTACT_TOPOLOGY_SIGNATURE_KIND = "sequence_only_repaired_contact_topology"
REPAIR_INPUT_BOUNDARY = "sequence_only_no_native_contacts_no_truth_axes"


@dataclass(frozen=True)
class LongRangeRepairPressure:
    helix_fraction: float
    hydrophobic_fraction: float
    charged_fraction: float
    compact_anchor_repair_enabled: bool
    disorder_like_blocked: bool
    membrane_like_blocked: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "helix_fraction": self.helix_fraction,
            "hydrophobic_fraction": self.hydrophobic_fraction,
            "charged_fraction": self.charged_fraction,
            "compact_anchor_repair_enabled": self.compact_anchor_repair_enabled,
            "disorder_like_blocked": self.disorder_like_blocked,
            "membrane_like_blocked": self.membrane_like_blocked,
        }


@dataclass(frozen=True)
class ContactRepairPacket:
    baseline_prediction: ContactTopologyPrediction
    repaired_prediction: ContactTopologyPrediction
    compact_anchor_candidate_count: int
    beta_pairing_candidate_count: int
    local_overclosure_trimmed_count: int
    repair_candidate_count: int
    repair_signature_kind: str
    repair_input_boundary: str
    native_truth_used_before_repair: bool
    raw_sequence_exposed: bool
    repair_notes: tuple[str, ...]

    def to_safe_dict(self) -> dict[str, object]:
        return {
            "row_id": self.repaired_prediction.row_id,
            "baseline_predicted_contact_count": len(
                self.baseline_prediction.predicted_contact_pairs
            ),
            "repaired_predicted_contact_count": len(
                self.repaired_prediction.predicted_contact_pairs
            ),
            "compact_anchor_candidate_count": self.compact_anchor_candidate_count,
            "beta_pairing_candidate_count": self.beta_pairing_candidate_count,
            "local_overclosure_trimmed_count": self.local_overclosure_trimmed_count,
            "repair_candidate_count": self.repair_candidate_count,
            "repair_signature_kind": self.repair_signature_kind,
            "repair_input_boundary": self.repair_input_boundary,
            "native_truth_used_before_repair": self.native_truth_used_before_repair,
            "raw_sequence_exposed": self.raw_sequence_exposed,
            "repair_notes": ";".join(self.repair_notes),
        }


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _fraction(sequence: str, alphabet: frozenset[str]) -> float:
    if not sequence:
        return 0.0
    return _rounded(sum(1 for residue in sequence if residue in alphabet) / len(sequence))


def long_range_repair_pressure(sequence: str) -> LongRangeRepairPressure:
    normalized = normalize_sequence(sequence)
    beta_pressure = beta_pairing_pressure(normalized)
    helix_fraction = _fraction(normalized, HELIX_RESIDUES)
    hydrophobic_fraction = _fraction(normalized, HYDROPHOBIC_AMINO_ACIDS)
    charged_fraction = _fraction(normalized, CHARGED_AMINO_ACIDS)
    compact_enabled = (
        len(normalized) >= 45
        and helix_fraction >= 0.43
        and hydrophobic_fraction < 0.50
        and not beta_pressure.disorder_like_blocked
        and not beta_pressure.membrane_like_blocked
    )
    return LongRangeRepairPressure(
        helix_fraction=helix_fraction,
        hydrophobic_fraction=hydrophobic_fraction,
        charged_fraction=charged_fraction,
        compact_anchor_repair_enabled=compact_enabled,
        disorder_like_blocked=beta_pressure.disorder_like_blocked,
        membrane_like_blocked=beta_pressure.membrane_like_blocked,
    )


def long_range_contact_repair_candidates(
    sequence: str,
) -> tuple[ContactCandidate, ...]:
    normalized = normalize_sequence(sequence)
    pressure = long_range_repair_pressure(normalized)
    if not pressure.compact_anchor_repair_enabled:
        return ()
    candidates = []
    for left, right in _compact_anchor_pairs(len(normalized)):
        candidates.append(
            ContactCandidate(
                i=left,
                j=right,
                score=_compact_anchor_score(normalized, left, right, pressure),
                contact_kind="compact_long_range_anchor_repair_candidate",
                evidence_reason=(
                    "sequence_only_compact_anchor_scan;"
                    f"helix_fraction={pressure.helix_fraction:.6f}"
                ),
            )
        )
    return tuple(candidates)


def repair_contact_topology(
    sequence: str,
    *,
    baseline_prediction: ContactTopologyPrediction,
) -> ContactRepairPacket:
    normalized = normalize_sequence(sequence)
    beta_pressure = beta_pairing_pressure(normalized)
    baseline_candidates = []
    trimmed = 0
    for candidate in baseline_prediction.candidates:
        if (
            beta_pressure.beta_pairing_repair_enabled
            and candidate.contact_kind == "local_backbone_closure"
            and candidate.score <= 0.54
        ):
            trimmed += 1
            continue
        baseline_candidates.append(candidate)

    compact_candidates = long_range_contact_repair_candidates(normalized)
    beta_candidates = beta_pairing_repair_candidates(normalized)
    merged_candidates = _dedupe_candidates(
        tuple(baseline_candidates) + compact_candidates + beta_candidates
    )
    merged_pairs = tuple(candidate.pair() for candidate in merged_candidates)
    repaired_prediction = ContactTopologyPrediction(
        row_id=baseline_prediction.row_id,
        sequence_hash=baseline_prediction.sequence_hash,
        sequence_length=baseline_prediction.sequence_length,
        contact_topology_signature_kind=REPAIRED_CONTACT_TOPOLOGY_SIGNATURE_KIND,
        predictor_input_boundary=REPAIR_INPUT_BOUNDARY,
        candidates=merged_candidates,
        predicted_contact_pairs=merged_pairs,
        predicted_contact_map_hash=contact_map_hash(merged_pairs),
        native_truth_used_before_prediction=False,
        raw_sequence_exposed=False,
    )
    notes = _repair_notes(
        compact_candidates=compact_candidates,
        beta_candidates=beta_candidates,
        trimmed=trimmed,
    )
    return ContactRepairPacket(
        baseline_prediction=baseline_prediction,
        repaired_prediction=repaired_prediction,
        compact_anchor_candidate_count=len(compact_candidates),
        beta_pairing_candidate_count=len(beta_candidates),
        local_overclosure_trimmed_count=trimmed,
        repair_candidate_count=len(compact_candidates) + len(beta_candidates),
        repair_signature_kind=(
            LONG_RANGE_CONTACT_REPAIR_SIGNATURE_KIND
            + ";"
            + BETA_PAIRING_REPAIR_SIGNATURE_KIND
        ),
        repair_input_boundary=REPAIR_INPUT_BOUNDARY,
        native_truth_used_before_repair=False,
        raw_sequence_exposed=False,
        repair_notes=notes,
    )


def _compact_anchor_pairs(length: int) -> tuple[tuple[int, int], ...]:
    pairs = []
    for start in range(6, length - 19, 18):
        pairs.append((start, start + 20))
    for start in range(12, length - 23, 38):
        pairs.append((start, start + 24))
    return tuple(sorted(set(pairs)))


def _compact_anchor_score(
    sequence: str,
    left: int,
    right: int,
    pressure: LongRangeRepairPressure,
) -> float:
    left_residue = sequence[left - 1]
    right_residue = sequence[right - 1]
    score = 0.57 + pressure.helix_fraction * 0.16
    if (
        left_residue in HYDROPHOBIC_AMINO_ACIDS
        and right_residue in HYDROPHOBIC_AMINO_ACIDS
    ):
        score += 0.05
    if left_residue in CHARGED_AMINO_ACIDS and right_residue in CHARGED_AMINO_ACIDS:
        score += 0.03
    return round(max(0.50, min(0.82, score)), 6)


def _dedupe_candidates(
    candidates: Sequence[ContactCandidate],
) -> tuple[ContactCandidate, ...]:
    by_pair: dict[tuple[int, int], ContactCandidate] = {}
    for candidate in candidates:
        pair = candidate.pair()
        previous = by_pair.get(pair)
        if previous is None or candidate.score > previous.score:
            by_pair[pair] = candidate
    return tuple(
        sorted(
            by_pair.values(),
            key=lambda candidate: (-candidate.score, candidate.i, candidate.j),
        )
    )


def _repair_notes(
    *,
    compact_candidates: Sequence[ContactCandidate],
    beta_candidates: Sequence[ContactCandidate],
    trimmed: int,
) -> tuple[str, ...]:
    notes = []
    if compact_candidates:
        notes.append("compact_anchor_repair_added")
    if beta_candidates:
        notes.append("beta_pairing_repair_added")
    if trimmed:
        notes.append("local_overclosure_trimmed")
    if not notes:
        notes.append("no_sequence_only_repair_gate_fired")
    return tuple(notes)


def repaired_contact_pair_set(packet: ContactRepairPacket) -> frozenset[tuple[int, int]]:
    return frozenset(packet.repaired_prediction.predicted_contact_pairs)


def baseline_contact_pair_set(packet: ContactRepairPacket) -> frozenset[tuple[int, int]]:
    return frozenset(packet.baseline_prediction.predicted_contact_pairs)
