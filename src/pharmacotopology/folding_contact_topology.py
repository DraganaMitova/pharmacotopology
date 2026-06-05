from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from pharmacotopology.folding_native_contact_eval import (
    ContactPair,
    contact_map_hash,
    normalized_contact_pairs,
)
from pharmacotopology.folding_topology import (
    AROMATIC_AMINO_ACIDS,
    CHARGED_AMINO_ACIDS,
    HYDROPHOBIC_AMINO_ACIDS,
    normalize_sequence,
)


VISUAL_MECHANISM_BENCHMARK_KIND = "visual_folding_mechanism_benchmark_v1"
VISUAL_MECHANISM_SPLIT = "visual_mechanism_12"
CONTACT_TOPOLOGY_SIGNATURE_KIND = "sequence_only_contact_topology_candidates"
PREDICTOR_INPUT_BOUNDARY = "sequence_only_no_native_contacts_no_truth_axes"

BETA_RESIDUES = frozenset("VIFYWTC")
HELIX_RESIDUES = frozenset("ALEMQKRH")
BREAKERS = frozenset("PG")


@dataclass(frozen=True)
class VisualMechanismRow:
    row_id: str
    source_id: str
    source_kind: str
    sequence: str
    sequence_sha256: str
    length: int
    truth_axes: dict[str, str]
    native_contact_pairs: tuple[ContactPair, ...]
    native_contact_map_hash: str
    native_scope: str
    mechanism_expected_difficulty: str

    def to_safe_dict(self) -> dict[str, object]:
        data = asdict(self)
        data.pop("sequence", None)
        return data


@dataclass(frozen=True)
class ContactCandidate:
    i: int
    j: int
    score: float
    contact_kind: str
    evidence_reason: str

    def pair(self) -> ContactPair:
        return (self.i, self.j)

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ContactTopologyPrediction:
    row_id: str
    sequence_hash: str
    sequence_length: int
    contact_topology_signature_kind: str
    predictor_input_boundary: str
    candidates: tuple[ContactCandidate, ...]
    predicted_contact_pairs: tuple[ContactPair, ...]
    predicted_contact_map_hash: str
    native_truth_used_before_prediction: bool
    raw_sequence_exposed: bool

    def to_safe_dict(self) -> dict[str, object]:
        return {
            "row_id": self.row_id,
            "sequence_hash": self.sequence_hash,
            "sequence_length": self.sequence_length,
            "contact_topology_signature_kind": self.contact_topology_signature_kind,
            "predictor_input_boundary": self.predictor_input_boundary,
            "predicted_contact_count": len(self.predicted_contact_pairs),
            "predicted_contact_map_hash": self.predicted_contact_map_hash,
            "native_truth_used_before_prediction": (
                self.native_truth_used_before_prediction
            ),
            "raw_sequence_exposed": self.raw_sequence_exposed,
        }


def sha256_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def short_sequence_hash(sequence: str) -> str:
    return sha256_sequence(sequence)[:16]


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _load_json(path: Path) -> Mapping[str, Any]:
    parsed = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(parsed, Mapping):
        raise ValueError(f"{path} must contain a JSON object")
    return parsed


def load_visual_mechanism_rows(path: Path) -> tuple[VisualMechanismRow, ...]:
    data = _load_json(path)
    rows = data.get("references")
    if not isinstance(rows, list):
        raise ValueError("visual mechanism benchmark must include references")
    output: list[VisualMechanismRow] = []
    for index, raw in enumerate(rows, start=1):
        if not isinstance(raw, Mapping):
            raise ValueError(f"row[{index}] must be an object")
        sequence = normalize_sequence(str(raw.get("sequence", "")))
        expected_sha = str(raw.get("sequence_sha256", ""))
        actual_sha = sha256_sequence(sequence)
        if expected_sha != actual_sha:
            raise ValueError(f"row[{index}].sequence_sha256_mismatch")
        length = int(raw.get("length", 0))
        if length != len(sequence):
            raise ValueError(f"row[{index}].length_mismatch")
        native_pairs = normalized_contact_pairs(raw.get("native_contact_pairs", ()))
        expected_hash = str(raw.get("native_contact_map_hash", ""))
        actual_hash = contact_map_hash(native_pairs)
        if expected_hash != actual_hash:
            raise ValueError(f"row[{index}].native_contact_map_hash_mismatch")
        truth_axes = raw.get("truth_axes", {})
        if not isinstance(truth_axes, Mapping):
            raise ValueError(f"row[{index}].truth_axes_must_be_object")
        output.append(
            VisualMechanismRow(
                row_id=str(raw["row_id"]),
                source_id=str(raw["source_id"]),
                source_kind=str(raw.get("source_kind", "locked_visual_mechanism_row")),
                sequence=sequence,
                sequence_sha256=actual_sha,
                length=length,
                truth_axes={str(key): str(value) for key, value in truth_axes.items()},
                native_contact_pairs=native_pairs,
                native_contact_map_hash=actual_hash,
                native_scope=str(raw["native_scope"]),
                mechanism_expected_difficulty=str(
                    raw["mechanism_expected_difficulty"]
                ),
            )
        )
    return tuple(output)


def validate_visual_mechanism_lock(rows: Sequence[VisualMechanismRow]) -> dict[str, object]:
    row_ids = [row.row_id for row in rows]
    sequence_hashes = [row.sequence_sha256 for row in rows]
    distribution: dict[str, int] = {}
    for row in rows:
        key = row.truth_axes.get("secondary_structure_axis", "unknown")
        if row.truth_axes.get("environment_axis") == "membrane_like":
            key = "membrane_like_control"
        elif row.truth_axes.get("order_axis") == "disordered_flexible":
            key = "disordered_flexible_control"
        elif row.truth_axes.get("architecture_axis") == "multidomain_or_segmented":
            key = "architecture_stress_control"
        distribution[key] = distribution.get(key, 0) + 1
    violations = []
    if len(rows) != 12:
        violations.append("visual_mechanism_row_count_not_12")
    if len(set(row_ids)) != len(row_ids):
        violations.append("duplicate_row_id")
    if len(set(sequence_hashes)) != len(sequence_hashes):
        violations.append("duplicate_sequence_sha256")
    return {
        "visual_mechanism_row_count": len(rows),
        "visual_mechanism_unique_sequence_count": len(set(sequence_hashes)),
        "visual_mechanism_distribution": distribution,
        "visual_mechanism_lock_valid": not violations,
        "visual_mechanism_lock_violations": tuple(violations),
    }


def _candidate_score(sequence: str, i: int, j: int) -> tuple[float, str, str]:
    left = sequence[i - 1]
    right = sequence[j - 1]
    separation = j - i
    if separation < 4:
        return 0.0, "ignored_neighbor", "near-neighbor contacts are ignored"

    score = 0.0
    reasons: list[str] = []
    kind = "weak_sequence_contact"
    if separation == 4 and left not in BREAKERS and right not in BREAKERS:
        score += 0.52
        kind = "local_backbone_closure"
        reasons.append("i_to_i_plus_4_nonbreaker_closure")
    if separation == 4 and left in HELIX_RESIDUES and right in HELIX_RESIDUES:
        score += 0.78
        kind = "local_helix_closure"
        reasons.append("i_to_i_plus_4_helix_pressure")
    if separation in {3, 5} and left in HELIX_RESIDUES and right in HELIX_RESIDUES:
        score += 0.45
        kind = "local_secondary_closure"
        reasons.append("near_helical_local_closure")
    if separation >= 14 and left in BETA_RESIDUES and right in BETA_RESIDUES:
        score += 0.38
        kind = "beta_pairing_candidate"
        reasons.append("long_range_beta_residue_pair")
        if (i + j) % 2 == 0:
            score += 0.18
            reasons.append("alternating_beta_register")
    if separation >= 16 and left in HYDROPHOBIC_AMINO_ACIDS and right in HYDROPHOBIC_AMINO_ACIDS:
        score += 0.30
        kind = "hydrophobic_core_candidate"
        reasons.append("long_range_hydrophobic_closure")
    if separation >= 10 and left in AROMATIC_AMINO_ACIDS and right in HYDROPHOBIC_AMINO_ACIDS:
        score += 0.16
        reasons.append("aromatic_hydrophobic_anchor")
    if separation >= 10 and right in AROMATIC_AMINO_ACIDS and left in HYDROPHOBIC_AMINO_ACIDS:
        score += 0.16
        reasons.append("aromatic_hydrophobic_anchor")
    if separation >= 8 and left in CHARGED_AMINO_ACIDS and right in CHARGED_AMINO_ACIDS:
        opposite = (left in "KRH" and right in "DE") or (right in "KRH" and left in "DE")
        if opposite:
            score += 0.20
            reasons.append("opposite_charge_bridge")
        else:
            score -= 0.12
            reasons.append("same_charge_penalty")
    if left in BREAKERS or right in BREAKERS:
        score -= 0.20
        reasons.append("breaker_penalty")
    if separation > 80:
        score -= 0.05
    return _rounded(score), kind, ";".join(reasons) or "weak sequence contact"


def predict_contact_topology(
    sequence: str,
    *,
    row_id: str = "sequence",
    max_contacts: int | None = None,
) -> ContactTopologyPrediction:
    normalized = normalize_sequence(sequence)
    length = len(normalized)
    candidates: list[ContactCandidate] = []
    stride = 1 if length <= 120 else 2
    for i in range(1, length + 1, stride):
        upper = min(length, i + 92)
        for j in range(i + 4, upper + 1):
            score, kind, reason = _candidate_score(normalized, i, j)
            if score >= 0.48:
                candidates.append(
                    ContactCandidate(
                        i=i,
                        j=j,
                        score=score,
                        contact_kind=kind,
                        evidence_reason=reason,
                    )
                )
    candidates.sort(key=lambda item: (-item.score, item.i, item.j))
    contact_limit = max_contacts or max(18, min(90, length // 2))
    selected = tuple(candidates[:contact_limit])
    predicted_pairs = tuple(candidate.pair() for candidate in selected)
    return ContactTopologyPrediction(
        row_id=row_id,
        sequence_hash=short_sequence_hash(normalized),
        sequence_length=length,
        contact_topology_signature_kind=CONTACT_TOPOLOGY_SIGNATURE_KIND,
        predictor_input_boundary=PREDICTOR_INPUT_BOUNDARY,
        candidates=selected,
        predicted_contact_pairs=predicted_pairs,
        predicted_contact_map_hash=contact_map_hash(predicted_pairs),
        native_truth_used_before_prediction=False,
        raw_sequence_exposed=False,
    )
