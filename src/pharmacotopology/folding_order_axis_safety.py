from __future__ import annotations

import hashlib
from dataclasses import asdict, dataclass
from typing import Mapping

from pharmacotopology.folding_axis_adjudication import UNKNOWN_BY_AXIS
from pharmacotopology.folding_hierarchical_gates import predict_hierarchical_gate
from pharmacotopology.folding_order_aware import extract_order_aware_features
from pharmacotopology.folding_regime_analysis import detect_protein_regime
from pharmacotopology.folding_topology import normalize_sequence


ORDER_AXIS_SAFETY_SIGNATURE_KIND = "external_safe_order_axis_quarantine_packet"
ORDER_AXIS_UNKNOWN = UNKNOWN_BY_AXIS["order_axis"]


@dataclass(frozen=True)
class OrderAxisSafetyPacket:
    row_id: str
    sequence_hash: str
    sequence_length: int
    protein_regime: str
    source_predicted_fold_class: str
    source_confidence: float
    disorder_pressure: float
    disorder_run_evidence: float
    breaker_density: float
    local_disorder_pressure: float
    charge_frustration_pressure: float
    folded_beta_mimic_pressure: float
    folded_mixed_mimic_pressure: float
    compact_closure_pressure: float
    beta_pairing_support: float
    long_range_closure_evidence: float
    contact_prior_density: float
    order_axis_prediction: str
    order_axis_claim_allowed: bool
    order_axis_confidence: float
    order_axis_abstention_reason: str
    order_axis_decision_reason: str

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _rounded(value: float) -> float:
    return round(max(0.0, min(1.0, value)), 6)


def _sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:16]


def order_axis_safety_packet_from_source(
    sequence: str,
    source_row: Mapping[str, object],
    *,
    row_id: str = "sequence",
) -> OrderAxisSafetyPacket:
    normalized = normalize_sequence(sequence)
    regime = detect_protein_regime(normalized, protein_id=row_id)
    hierarchy = predict_hierarchical_gate(normalized, protein_id=row_id)
    order_features, _edges = extract_order_aware_features(
        normalized,
        protein_id=row_id,
    )
    evidence = hierarchy.motif_evidence
    scores = hierarchy.gate_scores
    source_class = str(source_row["predicted_fold_class"])
    source_confidence = float(source_row["confidence"])
    protein_regime = str(source_row["protein_regime"])
    forced = bool(source_row["forced_prediction"])

    compact_closure_pressure = _rounded(
        scores.compactness_gate_score * 0.56
        + evidence.long_range_closure_evidence * 0.28
        + float(order_features["contact_prior_mean_weight"]) * 0.16
    )
    folded_beta_mimic_pressure = _rounded(
        evidence.local_beta_pressure_evidence * 0.28
        + scores.beta_pairing_support_score * 0.30
        + evidence.long_range_closure_evidence * 0.24
        + float(order_features["contact_prior_density"]) * 0.18
    )
    folded_mixed_mimic_pressure = _rounded(
        min(scores.alpha_structure_score, scores.beta_structure_score) * 0.40
        + scores.mixed_structure_score * 0.30
        + evidence.long_range_closure_evidence * 0.18
        + float(order_features["contact_prior_density"]) * 0.12
    )

    prediction = ORDER_AXIS_UNKNOWN
    claim_allowed = False
    confidence = 0.0
    abstention_reason = "order axis evidence is insufficient or mixed"
    decision_reason = abstention_reason

    if forced and source_class in {
        "alpha_rich",
        "beta_rich",
        "alpha_beta_mixed",
        "multidomain_boundary",
    }:
        prediction = "ordered"
        claim_allowed = True
        confidence = _rounded(source_confidence)
        abstention_reason = ""
        decision_reason = "forced folded source class is allowed as ordered axis"
    elif forced and source_class == "disordered_flexible":
        strong_disorder = (
            regime.intrinsically_disordered_pressure >= 0.80
            and evidence.disorder_run_evidence >= 0.10
            and float(order_features["breaker_density"]) >= 0.24
            and evidence.local_disorder_pressure_evidence >= 0.70
            and compact_closure_pressure <= 0.30
            and scores.beta_pairing_support_score <= 0.22
            and folded_beta_mimic_pressure <= 0.22
            and folded_mixed_mimic_pressure <= 0.22
        )
        if strong_disorder:
            prediction = "disordered_flexible"
            claim_allowed = True
            confidence = _rounded(
                regime.intrinsically_disordered_pressure * 0.44
                + evidence.disorder_run_evidence * 0.22
                + float(order_features["breaker_density"]) * 0.20
                + (1.0 - compact_closure_pressure) * 0.14
            )
            abstention_reason = ""
            decision_reason = (
                "disorder order-axis claim has strong run/breaker evidence and "
                "weak folded-domain mimic evidence"
            )
        else:
            abstention_reason = "external_order_axis_folded_mimic_quarantine"
            decision_reason = (
                "disordered_flexible was not projected to order_axis because "
                "folded-domain mimic or weak disorder-run evidence remains"
            )
    elif not forced:
        abstention_reason = "source did not force a fold class for order-axis projection"
        decision_reason = abstention_reason

    return OrderAxisSafetyPacket(
        row_id=row_id,
        sequence_hash=_sequence_hash(normalized),
        sequence_length=len(normalized),
        protein_regime=protein_regime,
        source_predicted_fold_class=source_class,
        source_confidence=source_confidence,
        disorder_pressure=regime.intrinsically_disordered_pressure,
        disorder_run_evidence=evidence.disorder_run_evidence,
        breaker_density=_rounded(float(order_features["breaker_density"])),
        local_disorder_pressure=evidence.local_disorder_pressure_evidence,
        charge_frustration_pressure=evidence.charge_frustration_evidence,
        folded_beta_mimic_pressure=folded_beta_mimic_pressure,
        folded_mixed_mimic_pressure=folded_mixed_mimic_pressure,
        compact_closure_pressure=compact_closure_pressure,
        beta_pairing_support=scores.beta_pairing_support_score,
        long_range_closure_evidence=evidence.long_range_closure_evidence,
        contact_prior_density=_rounded(float(order_features["contact_prior_density"])),
        order_axis_prediction=prediction,
        order_axis_claim_allowed=claim_allowed,
        order_axis_confidence=confidence,
        order_axis_abstention_reason=abstention_reason,
        order_axis_decision_reason=decision_reason,
    )
