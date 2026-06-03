from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_motif_alignment import (
    ABSTAINED_CLASS,
    MOTIF_ALIGNMENT_BENCHMARK_KIND,
    MOTIF_SIGNATURE_KIND,
    MotifEvidenceVector,
    load_motif_alignment_inputs,
    predict_motif_alignment,
)
from pharmacotopology.folding_structure_benchmark import (
    LABEL_BENCHMARK_KIND_V0,
    StructureEvidenceRow,
)
from pharmacotopology.folding_topology import (
    CONTACT_PROXY_DIMENSIONS,
    FoldingReferenceExample,
    FoldingTopologySignature,
    contact_map_proxy_similarity,
    normalize_fold_class,
    signature_to_dict,
)


HIERARCHICAL_GATE_BENCHMARK_KIND = "hierarchical_folding_decision_gate_benchmark"
HIERARCHICAL_GATE_SIGNATURE_KIND = "motif_evidence_hierarchical_gate_interpretation"


@dataclass(frozen=True)
class HierarchicalGateScores:
    disorder_gate_score: float
    compactness_gate_score: float
    segmentation_gate_score: float
    secondary_structure_gate_score: float
    alpha_structure_score: float
    beta_structure_score: float
    mixed_structure_score: float
    beta_pairing_support_score: float

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class HierarchicalGatePrediction:
    protein_id: str
    sequence_length: int
    motif_evidence: MotifEvidenceVector
    topology_signature: FoldingTopologySignature
    raw_motif_predicted_fold_class: str
    motif_predicted_fold_class: str
    predicted_fold_class: str
    gate_scores: HierarchicalGateScores
    gate_path: tuple[str, ...]
    gate_decision_reason: str
    flexible_segmentation_warning: bool
    beta_evidence_requires_pairing: bool
    alpha_evidence_requires_periodicity: bool
    hierarchy_prediction_used: bool
    hierarchy_changed_raw_prediction: bool
    hierarchy_prevented_false_multidomain: bool
    hierarchy_prevented_false_beta: bool
    hierarchy_prevented_false_mixed: bool
    evidence_conflict_score: float
    confidence: float
    uncertainty_radius: float
    claim_strength: str
    forced_prediction: bool
    abstained: bool
    motif_conflict_score: float


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


def _rounded(value: float) -> float:
    return round(_clamp(value), 6)


def _mean(values: Iterable[float]) -> float:
    values = tuple(values)
    if not values:
        return 0.0
    return sum(values) / len(values)


def _gate_scores(evidence: MotifEvidenceVector) -> HierarchicalGateScores:
    compactness = _rounded(
        evidence.compact_core_evidence * 0.56
        + evidence.long_range_closure_evidence * 0.34
        + (1.0 - evidence.local_disorder_pressure_evidence) * 0.10
        - evidence.breaker_turn_evidence * 0.08
    )
    disorder = _rounded(
        evidence.local_disorder_pressure_evidence * 0.52
        + evidence.breaker_turn_evidence * 0.22
        + evidence.disorder_run_evidence * 0.16
        + evidence.charge_frustration_evidence * 0.05
        + evidence.domain_boundary_evidence * 0.05
        - evidence.compact_core_evidence * 0.19
        - evidence.long_range_closure_evidence * 0.10
    )
    segmentation = _rounded(
        evidence.domain_boundary_evidence * 0.46
        + compactness * 0.24
        + evidence.long_range_closure_evidence * 0.20
        - evidence.local_disorder_pressure_evidence * 0.18
        + evidence.charge_frustration_evidence * 0.04
    )
    beta_pairing = _rounded(
        evidence.long_range_closure_evidence * 0.38
        + evidence.compact_core_evidence * 0.26
        + evidence.local_beta_pressure_evidence * 0.22
        + evidence.beta_alternation_evidence * 0.08
        - evidence.local_disorder_pressure_evidence * 0.14
    )
    alpha = _rounded(
        evidence.local_alpha_pressure_evidence * 0.48
        + evidence.alpha_periodicity_evidence * 0.30
        + evidence.compact_core_evidence * 0.14
        + evidence.long_range_closure_evidence * 0.04
        - evidence.local_beta_pressure_evidence * 0.07
        - evidence.local_disorder_pressure_evidence * 0.10
    )
    beta = _rounded(
        evidence.local_beta_pressure_evidence * 0.34
        + evidence.beta_alternation_evidence * 0.21
        + beta_pairing * 0.30
        + evidence.long_range_closure_evidence * 0.08
        - evidence.local_alpha_pressure_evidence * 0.08
        - evidence.local_disorder_pressure_evidence * 0.12
    )
    mixed = _rounded(
        min(alpha, beta) * 0.55
        + evidence.mixed_motif_evidence * 0.24
        + compactness * 0.12
        - abs(alpha - beta) * 0.12
    )
    return HierarchicalGateScores(
        disorder_gate_score=disorder,
        compactness_gate_score=compactness,
        segmentation_gate_score=segmentation,
        secondary_structure_gate_score=max(alpha, beta, mixed),
        alpha_structure_score=alpha,
        beta_structure_score=beta,
        mixed_structure_score=mixed,
        beta_pairing_support_score=beta_pairing,
    )


def _claim_strength(confidence: float, abstained: bool) -> str:
    if abstained:
        return "abstained"
    if confidence >= 0.72:
        return "strong"
    if confidence >= 0.58:
        return "medium"
    return "weak"


def _hierarchy_conflict(confidence: float, abstained: bool) -> float:
    if abstained:
        return 0.86
    return _rounded(0.66 - confidence * 0.60)


def _append_confidence_gate(
    gate_path: list[str],
    *,
    candidate: str,
    confidence: float,
    reason: str,
) -> tuple[str, bool, str]:
    if candidate == ABSTAINED_CLASS:
        gate_path.append("confidence_gate:abstained")
        return candidate, True, reason
    if confidence < 0.34:
        gate_path.append("confidence_gate:abstained_low_decision_margin")
        return (
            ABSTAINED_CLASS,
            True,
            f"{reason}; confidence gate abstained on low decision margin",
        )
    gate_path.append("confidence_gate:forced_low_claim_strength")
    return candidate, False, reason


def predict_hierarchical_gate(
    sequence: str,
    *,
    protein_id: str = "sequence",
) -> HierarchicalGatePrediction:
    motif_prediction = predict_motif_alignment(sequence, protein_id=protein_id)
    evidence = motif_prediction.motif_evidence
    scores = _gate_scores(evidence)
    gate_path: list[str] = []
    candidate = ABSTAINED_CLASS
    confidence = 0.0
    reason = "no hierarchical gate reached a class interpretation"
    flexible_segmentation_warning = (
        scores.disorder_gate_score >= 0.30
        and evidence.domain_boundary_evidence >= 0.36
        and scores.compactness_gate_score < 0.45
    )

    if (
        scores.disorder_gate_score >= 0.30
        and scores.compactness_gate_score < 0.45
    ):
        candidate = "disordered_flexible"
        confidence = _rounded(
            0.38
            + (scores.disorder_gate_score - 0.30) * 1.10
            + (0.45 - scores.compactness_gate_score) * 0.80
            + evidence.local_disorder_pressure_evidence * 0.16
        )
        gate_path.extend(
            (
                "disorder_gate:disordered_or_flexible",
                "compactness_gate:not_compact_enough_for_folded_domain",
                "segmentation_gate:flexible_segmentation_not_multidomain",
                "secondary_structure_gate:skipped_by_disorder_gate",
            )
        )
        reason = (
            "disorder pressure and breaker/turn evidence are prioritized over "
            "secondary-structure interpretation because compact closure is weak"
        )
    elif scores.compactness_gate_score < 0.42:
        gate_path.extend(
            (
                "disorder_gate:foldable_not_ruled_out",
                "compactness_gate:weak_closure_abstain",
                "segmentation_gate:not_evaluated",
                "secondary_structure_gate:not_evaluated",
            )
        )
        confidence = 0.0
        reason = "compactness and closure are too weak to force a fold class"
    elif (
        evidence.domain_boundary_evidence >= 0.43
        and scores.compactness_gate_score >= 0.49
        and evidence.long_range_closure_evidence >= 0.50
        and scores.disorder_gate_score < 0.24
    ):
        candidate = "multidomain_boundary"
        confidence = _rounded(
            0.38
            + (evidence.domain_boundary_evidence - 0.43) * 1.40
            + (scores.compactness_gate_score - 0.49) * 0.70
            + (evidence.long_range_closure_evidence - 0.50) * 0.55
        )
        gate_path.extend(
            (
                "disorder_gate:foldable_candidate",
                "compactness_gate:compact_closure_supported",
                "segmentation_gate:folded_multidomain_supported",
                "secondary_structure_gate:skipped_by_segmentation_gate",
            )
        )
        reason = (
            "domain-boundary evidence is treated as folded multidomain only "
            "because compact core and long-range closure are also present"
        )
    else:
        gate_path.extend(
            (
                "disorder_gate:foldable_candidate",
                "compactness_gate:compact_or_borderline_supported",
                "segmentation_gate:not_folded_multidomain",
            )
        )
        alpha_periodic = (
            evidence.local_alpha_pressure_evidence >= 0.54
            and evidence.alpha_periodicity_evidence >= 0.42
            and scores.compactness_gate_score >= 0.50
            and evidence.local_disorder_pressure_evidence < 0.34
        )
        beta_pairable = (
            scores.beta_structure_score >= 0.37
            and scores.beta_pairing_support_score >= 0.39
            and evidence.local_beta_pressure_evidence >= 0.47
            and scores.alpha_structure_score < 0.34
            and evidence.local_disorder_pressure_evidence < 0.50
        )
        mixed_supported = (
            scores.alpha_structure_score >= 0.32
            and scores.beta_structure_score >= 0.34
            and scores.compactness_gate_score >= 0.47
            and evidence.local_disorder_pressure_evidence < 0.42
        )
        if alpha_periodic:
            candidate = "alpha_rich"
            confidence = _rounded(
                0.42
                + max(
                    evidence.local_alpha_pressure_evidence
                    - evidence.local_beta_pressure_evidence,
                    0.0,
                )
                * 0.75
                + evidence.alpha_periodicity_evidence * 0.10
                + max(scores.compactness_gate_score - 0.50, 0.0) * 0.20
            )
            gate_path.append("secondary_structure_gate:alpha_periodic_compact")
            reason = (
                "alpha interpretation requires local alpha pressure, periodic "
                "support, compact core, and limited disorder pressure"
            )
        elif beta_pairable:
            candidate = "beta_rich"
            confidence = _rounded(
                0.36
                + max(
                    scores.beta_structure_score - scores.alpha_structure_score,
                    0.0,
                )
                * 0.60
                + max(scores.beta_pairing_support_score - 0.39, 0.0) * 0.40
            )
            gate_path.append("secondary_structure_gate:beta_pairing_supported")
            reason = (
                "beta interpretation requires alternation plus local beta "
                "pressure and sheet/contact-pairing support"
            )
        elif mixed_supported:
            candidate = "alpha_beta_mixed"
            confidence = _rounded(
                min(scores.alpha_structure_score, scores.beta_structure_score) * 0.50
                + scores.compactness_gate_score * 0.20
                - abs(scores.alpha_structure_score - scores.beta_structure_score)
                * 0.50
            )
            gate_path.append("secondary_structure_gate:mixed_supported_but_weak")
            reason = (
                "mixed alpha/beta evidence exists, but the confidence gate must "
                "decide whether the local signals are stable enough"
            )
        else:
            top_secondary = max(
                (
                    ("alpha_rich", scores.alpha_structure_score),
                    ("beta_rich", scores.beta_structure_score),
                    ("alpha_beta_mixed", scores.mixed_structure_score),
                ),
                key=lambda item: item[1],
            )
            candidate = top_secondary[0]
            confidence = _rounded(top_secondary[1] * 0.82)
            gate_path.append("secondary_structure_gate:weak_secondary_fallback")
            reason = (
                "secondary structure evidence was too weakly separated for a "
                "strong gate; fallback remains subject to abstention"
            )

    predicted, abstained, reason = _append_confidence_gate(
        gate_path,
        candidate=candidate,
        confidence=confidence,
        reason=reason,
    )
    forced = not abstained
    hierarchy_changed_raw = predicted != motif_prediction.raw_predicted_fold_class
    return HierarchicalGatePrediction(
        protein_id=protein_id,
        sequence_length=motif_prediction.sequence_length,
        motif_evidence=evidence,
        topology_signature=motif_prediction.topology_signature,
        raw_motif_predicted_fold_class=motif_prediction.raw_predicted_fold_class,
        motif_predicted_fold_class=motif_prediction.predicted_fold_class,
        predicted_fold_class=predicted,
        gate_scores=scores,
        gate_path=tuple(gate_path),
        gate_decision_reason=reason,
        flexible_segmentation_warning=flexible_segmentation_warning,
        beta_evidence_requires_pairing=True,
        alpha_evidence_requires_periodicity=True,
        hierarchy_prediction_used=True,
        hierarchy_changed_raw_prediction=hierarchy_changed_raw,
        hierarchy_prevented_false_multidomain=False,
        hierarchy_prevented_false_beta=False,
        hierarchy_prevented_false_mixed=False,
        evidence_conflict_score=_hierarchy_conflict(confidence, abstained),
        confidence=confidence,
        uncertainty_radius=_rounded(
            0.18
            + _hierarchy_conflict(confidence, abstained) * 0.52
            + (0.0 if forced else 0.16)
        ),
        claim_strength=_claim_strength(confidence, abstained),
        forced_prediction=forced,
        abstained=abstained,
        motif_conflict_score=motif_prediction.evidence_conflict_score,
    )


def _with_failure_mode_flags(
    prediction: HierarchicalGatePrediction,
    *,
    structure_class: str,
) -> HierarchicalGatePrediction:
    prevented_multidomain = (
        prediction.raw_motif_predicted_fold_class == "multidomain_boundary"
        and structure_class != "multidomain_boundary"
        and prediction.predicted_fold_class != "multidomain_boundary"
    )
    prevented_beta = (
        prediction.raw_motif_predicted_fold_class == "beta_rich"
        and structure_class == "disordered_flexible"
        and prediction.predicted_fold_class != "beta_rich"
    )
    prevented_mixed = (
        prediction.raw_motif_predicted_fold_class == "alpha_beta_mixed"
        and structure_class == "alpha_rich"
        and prediction.predicted_fold_class != "alpha_beta_mixed"
    )
    return HierarchicalGatePrediction(
        **{
            **asdict(prediction),
            "motif_evidence": prediction.motif_evidence,
            "topology_signature": prediction.topology_signature,
            "gate_scores": prediction.gate_scores,
            "gate_path": prediction.gate_path,
            "hierarchy_prevented_false_multidomain": prevented_multidomain,
            "hierarchy_prevented_false_beta": prevented_beta,
            "hierarchy_prevented_false_mixed": prevented_mixed,
        }
    )


def _failure_reason(row: Mapping[str, object]) -> str:
    if row["prediction_structure_class_match"]:
        return ""
    if row["abstained"]:
        return str(row["gate_decision_reason"])
    if (
        row["predicted_fold_class"] == "multidomain_boundary"
        and row["structure_fold_class"] == "disordered_flexible"
    ):
        return "flexible_segmentation_still_misread_as_multidomain"
    if (
        row["predicted_fold_class"] == "beta_rich"
        and row["structure_fold_class"] == "disordered_flexible"
    ):
        return "disorder_still_misread_as_beta"
    if (
        row["predicted_fold_class"] == "alpha_beta_mixed"
        and row["structure_fold_class"] == "alpha_rich"
    ):
        return "alpha_rich_still_misread_as_mixed"
    return "hierarchical_gate_prediction_structure_mismatch"


def hierarchical_gate_rows(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
) -> list[dict[str, object]]:
    evidence_by_id = {row.protein_id: row for row in evidence_rows}
    rows: list[dict[str, object]] = []
    for reference in references:
        structure = evidence_by_id[reference.protein_id]
        label_class = normalize_fold_class(reference.reference_fold_class)
        structure_class = structure.structure_fold_class
        prediction = _with_failure_mode_flags(
            predict_hierarchical_gate(
                reference.sequence,
                protein_id=reference.protein_id,
            ),
            structure_class=structure_class,
        )
        class_match = prediction.predicted_fold_class == structure_class
        label_match = prediction.predicted_fold_class == label_class
        similarity = contact_map_proxy_similarity(
            prediction.topology_signature,
            structure.structure_topology_signature,
            CONTACT_PROXY_DIMENSIONS,
        )
        label_similarity = contact_map_proxy_similarity(
            prediction.topology_signature,
            reference.reference_topology_signature,
            CONTACT_PROXY_DIMENSIONS,
        )
        high_confidence_wrong = (
            prediction.forced_prediction
            and not class_match
            and prediction.confidence >= 0.58
        )
        evidence = prediction.motif_evidence.to_dict()
        scores = prediction.gate_scores.to_dict()
        row: dict[str, object] = {
            "protein_id": reference.protein_id,
            "sequence_length": prediction.sequence_length,
            "topology_evidence_vector_kind": MOTIF_SIGNATURE_KIND,
            "hierarchical_gate_signature_kind": HIERARCHICAL_GATE_SIGNATURE_KIND,
            "alpha_periodicity_evidence": evidence["alpha_periodicity_evidence"],
            "beta_alternation_evidence": evidence["beta_alternation_evidence"],
            "compact_core_evidence": evidence["compact_core_evidence"],
            "disorder_run_evidence": evidence["disorder_run_evidence"],
            "domain_boundary_evidence": evidence["domain_boundary_evidence"],
            "long_range_closure_evidence": evidence["long_range_closure_evidence"],
            "breaker_turn_evidence": evidence["breaker_turn_evidence"],
            "charge_frustration_evidence": evidence[
                "charge_frustration_evidence"
            ],
            "local_alpha_pressure_evidence": evidence[
                "local_alpha_pressure_evidence"
            ],
            "local_beta_pressure_evidence": evidence[
                "local_beta_pressure_evidence"
            ],
            "local_disorder_pressure_evidence": evidence[
                "local_disorder_pressure_evidence"
            ],
            "mixed_motif_evidence": evidence["mixed_motif_evidence"],
            "disorder_gate_score": scores["disorder_gate_score"],
            "compactness_gate_score": scores["compactness_gate_score"],
            "segmentation_gate_score": scores["segmentation_gate_score"],
            "secondary_structure_gate_score": scores[
                "secondary_structure_gate_score"
            ],
            "alpha_structure_score": scores["alpha_structure_score"],
            "beta_structure_score": scores["beta_structure_score"],
            "mixed_structure_score": scores["mixed_structure_score"],
            "beta_pairing_support_score": scores["beta_pairing_support_score"],
            "gate_path": " | ".join(prediction.gate_path),
            "gate_decision_reason": prediction.gate_decision_reason,
            "flexible_segmentation_warning": (
                prediction.flexible_segmentation_warning
            ),
            "beta_evidence_requires_pairing": (
                prediction.beta_evidence_requires_pairing
            ),
            "alpha_evidence_requires_periodicity": (
                prediction.alpha_evidence_requires_periodicity
            ),
            "hierarchy_prediction_used": prediction.hierarchy_prediction_used,
            "hierarchy_changed_raw_prediction": (
                prediction.hierarchy_changed_raw_prediction
            ),
            "hierarchy_prevented_false_multidomain": (
                prediction.hierarchy_prevented_false_multidomain
            ),
            "hierarchy_prevented_false_beta": (
                prediction.hierarchy_prevented_false_beta
            ),
            "hierarchy_prevented_false_mixed": (
                prediction.hierarchy_prevented_false_mixed
            ),
            "motif_predicted_fold_class": prediction.motif_predicted_fold_class,
            "raw_motif_predicted_fold_class": (
                prediction.raw_motif_predicted_fold_class
            ),
            "predicted_fold_class": prediction.predicted_fold_class,
            "structure_fold_class": structure_class,
            "label_fold_class": label_class,
            "prediction_vs_structure_score": similarity,
            "prediction_vs_label_score": label_similarity,
            "prediction_structure_class_match": class_match,
            "prediction_label_class_match": label_match,
            "evidence_conflict_score": prediction.evidence_conflict_score,
            "motif_conflict_score": prediction.motif_conflict_score,
            "confidence": prediction.confidence,
            "uncertainty_radius": prediction.uncertainty_radius,
            "claim_strength": prediction.claim_strength,
            "forced_prediction": prediction.forced_prediction,
            "abstained": prediction.abstained,
            "high_confidence_wrong": high_confidence_wrong,
            "folding_problem_solved": False,
            "folding_solution_claim_created": False,
        }
        row["likely_failure_reason"] = _failure_reason(row)
        rows.append(row)
    return rows


def gate_path_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    return [
        {
            "protein_id": row["protein_id"],
            "raw_motif_predicted_fold_class": row[
                "raw_motif_predicted_fold_class"
            ],
            "predicted_fold_class": row["predicted_fold_class"],
            "structure_fold_class": row["structure_fold_class"],
            "disorder_gate_score": row["disorder_gate_score"],
            "compactness_gate_score": row["compactness_gate_score"],
            "segmentation_gate_score": row["segmentation_gate_score"],
            "secondary_structure_gate_score": row[
                "secondary_structure_gate_score"
            ],
            "gate_path": row["gate_path"],
            "gate_decision_reason": row["gate_decision_reason"],
            "forced_prediction": row["forced_prediction"],
            "abstained": row["abstained"],
        }
        for row in rows
    ]


def gate_failure_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    return [
        {
            "protein_id": row["protein_id"],
            "predicted_fold_class": row["predicted_fold_class"],
            "raw_motif_predicted_fold_class": row[
                "raw_motif_predicted_fold_class"
            ],
            "structure_fold_class": row["structure_fold_class"],
            "label_fold_class": row["label_fold_class"],
            "gate_path": row["gate_path"],
            "gate_decision_reason": row["gate_decision_reason"],
            "likely_failure_reason": row["likely_failure_reason"],
            "confidence": row["confidence"],
            "uncertainty_radius": row["uncertainty_radius"],
            "claim_strength": row["claim_strength"],
            "high_confidence_wrong": row["high_confidence_wrong"],
            "abstained": row["abstained"],
        }
        for row in rows
        if not bool(row["prediction_structure_class_match"])
    ]


def _mean_from_rows(rows: Sequence[Mapping[str, object]], key: str) -> float:
    return round(_mean(float(row.get(key, 0.0)) for row in rows), 6)


def _bool_mean(values: Iterable[bool]) -> float:
    values = tuple(values)
    if not values:
        return 0.0
    return round(sum(1 for value in values if value) / len(values), 6)


def build_hierarchical_gate_report(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
    *,
    source_benchmark_file: Path,
    structure_evidence_file: Path,
) -> dict[str, object]:
    rows = hierarchical_gate_rows(references, evidence_rows)
    forced_rows = [row for row in rows if bool(row["forced_prediction"])]
    structure_matches = sum(
        1 for row in rows if bool(row["prediction_structure_class_match"])
    )
    label_matches = sum(1 for row in rows if bool(row["prediction_label_class_match"]))
    disorder_gate_accuracy = _bool_mean(
        (
            str(row["structure_fold_class"]) == "disordered_flexible"
        )
        == str(row["gate_path"]).startswith("disorder_gate:disordered_or_flexible")
        for row in rows
    )
    compactness_gate_accuracy = _bool_mean(
        (
            str(row["structure_fold_class"]) != "disordered_flexible"
        )
        == (
            float(row["compactness_gate_score"]) >= 0.42
            and not str(row["gate_path"]).startswith(
                "disorder_gate:disordered_or_flexible"
            )
        )
        for row in rows
    )
    segmentation_gate_accuracy = _bool_mean(
        (
            str(row["structure_fold_class"]) == "multidomain_boundary"
        )
        == (str(row["predicted_fold_class"]) == "multidomain_boundary")
        for row in rows
    )
    secondary_rows = [
        row
        for row in rows
        if row["structure_fold_class"]
        in {"alpha_rich", "beta_rich", "alpha_beta_mixed"}
    ]
    secondary_structure_gate_accuracy = _bool_mean(
        row["prediction_structure_class_match"] for row in secondary_rows
    )
    high_confidence_wrong = sum(
        1 for row in rows if bool(row["high_confidence_wrong"])
    )
    return {
        "benchmark_kind": HIERARCHICAL_GATE_BENCHMARK_KIND,
        "source_motif_alignment_benchmark_kind": MOTIF_ALIGNMENT_BENCHMARK_KIND,
        "source_label_benchmark_kind": LABEL_BENCHMARK_KIND_V0,
        "topology_evidence_vector_kind": MOTIF_SIGNATURE_KIND,
        "hierarchical_gate_signature_kind": HIERARCHICAL_GATE_SIGNATURE_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "structure_evidence_file": str(structure_evidence_file),
        "predictor_input_boundary": "sequence_only_no_labels_no_structure_answers",
        "truth_channels_used_only_after_prediction": True,
        "benchmark_size": len(rows),
        "prediction_vs_structure_accuracy": round(
            structure_matches / max(len(rows), 1),
            6,
        ),
        "prediction_vs_label_accuracy": round(label_matches / max(len(rows), 1), 6),
        "forced_prediction_accuracy": round(
            sum(1 for row in forced_rows if bool(row["prediction_structure_class_match"]))
            / max(len(forced_rows), 1),
            6,
        ),
        "forced_prediction_count": len(forced_rows),
        "abstained_prediction_count": sum(1 for row in rows if bool(row["abstained"])),
        "high_confidence_wrong_count": high_confidence_wrong,
        "evidence_conflict_mean": _mean_from_rows(rows, "evidence_conflict_score"),
        "motif_conflict_mean": _mean_from_rows(rows, "motif_conflict_score"),
        "disorder_gate_accuracy": disorder_gate_accuracy,
        "compactness_gate_accuracy": compactness_gate_accuracy,
        "segmentation_gate_accuracy": segmentation_gate_accuracy,
        "secondary_structure_gate_accuracy": secondary_structure_gate_accuracy,
        "flexible_segmentation_false_multidomain_count": sum(
            1
            for row in rows
            if row["structure_fold_class"] == "disordered_flexible"
            and row["predicted_fold_class"] == "multidomain_boundary"
        ),
        "false_beta_from_disorder_count": sum(
            1
            for row in rows
            if row["structure_fold_class"] == "disordered_flexible"
            and row["predicted_fold_class"] == "beta_rich"
        ),
        "false_mixed_from_alpha_count": sum(
            1
            for row in rows
            if row["structure_fold_class"] == "alpha_rich"
            and row["predicted_fold_class"] == "alpha_beta_mixed"
        ),
        "hierarchy_changed_raw_prediction_count": sum(
            1 for row in rows if bool(row["hierarchy_changed_raw_prediction"])
        ),
        "hierarchy_prevented_false_multidomain_count": sum(
            1 for row in rows if bool(row["hierarchy_prevented_false_multidomain"])
        ),
        "hierarchy_prevented_false_beta_count": sum(
            1 for row in rows if bool(row["hierarchy_prevented_false_beta"])
        ),
        "hierarchy_prevented_false_mixed_count": sum(
            1 for row in rows if bool(row["hierarchy_prevented_false_mixed"])
        ),
        "revision_required": True,
        "claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "boundary_statement": (
            "This layer interprets the sequence-only motif evidence through "
            "staged disorder, compactness, segmentation, secondary-structure, "
            "and confidence gates. Reference labels and structure classes are "
            "used only after prediction for scoring and failure diagnosis."
        ),
        "rows": rows,
    }


def write_hierarchical_gate_outputs(
    *,
    report: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    path_rows: Sequence[Mapping[str, object]],
    failure_rows: Sequence[Mapping[str, object]],
    report_path: Path,
    rows_path: Path,
    gate_paths_path: Path,
    gate_failures_path: Path,
    dashboard_path: Path,
) -> tuple[Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(
        rows,
        rows_path,
        nested_json_fields=set(),
    )
    _write_csv_rows(path_rows, gate_paths_path, nested_json_fields=set())
    _write_csv_rows(failure_rows, gate_failures_path, nested_json_fields=set())
    dashboard_path.write_text(
        render_hierarchical_gate_dashboard(report),
        encoding="utf-8",
    )
    return report_path, rows_path, gate_paths_path, gate_failures_path, dashboard_path


def _write_csv_rows(
    rows: Sequence[Mapping[str, object]],
    path: Path,
    *,
    nested_json_fields: set[str],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as file:
        if not rows:
            file.write("")
            return path
        fieldnames = list(rows[0])
        writer = csv.DictWriter(file, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            output = {}
            for key, value in row.items():
                if key in nested_json_fields:
                    output[key] = json.dumps(value, sort_keys=True, separators=(",", ":"))
                else:
                    output[key] = value
            writer.writerow(output)
    return path


def _escape(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _metric_cards(report: Mapping[str, object]) -> str:
    labels = (
        "prediction_vs_structure_accuracy",
        "forced_prediction_count",
        "abstained_prediction_count",
        "high_confidence_wrong_count",
        "accuracy_delta_from_10",
        "stability_status",
        "false_beta_from_disorder_count",
        "false_mixed_from_alpha_count",
        "flexible_segmentation_false_multidomain_count",
    )
    return "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(report.get(label))}</strong>"
        "</div>"
        for label in labels
    )


def _mapping_table(title: str, mapping: object) -> str:
    if not isinstance(mapping, Mapping) or not mapping:
        return ""
    rows = []
    for key, value in mapping.items():
        if isinstance(value, Mapping):
            nested = ", ".join(
                f"{_escape(nested_key)}: {_escape(nested_value)}"
                for nested_key, nested_value in value.items()
            )
            rendered_value = nested
        else:
            rendered_value = _escape(value)
        rows.append(
            "<tr>"
            f"<td>{_escape(key)}</td>"
            f"<td>{rendered_value}</td>"
            "</tr>"
        )
    return (
        f"<section><h2>{_escape(title)}</h2>"
        "<table><thead><tr><th>key</th><th>value</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></section>"
    )


def _confusion_matrix_table(report: Mapping[str, object]) -> str:
    matrix = report.get("confusion_matrix")
    if not isinstance(matrix, Mapping) or not matrix:
        return ""
    predicted = sorted(
        {
            str(predicted_class)
            for row in matrix.values()
            if isinstance(row, Mapping)
            for predicted_class in row
        }
    )
    header = "".join(f"<th>{_escape(label)}</th>" for label in predicted)
    body = []
    for actual, row in matrix.items():
        if not isinstance(row, Mapping):
            continue
        cells = "".join(
            f"<td>{_escape(row.get(predicted_class, 0))}</td>"
            for predicted_class in predicted
        )
        body.append(f"<tr><td>{_escape(actual)}</td>{cells}</tr>")
    return (
        "<section><h2>Confusion Matrix</h2>"
        "<table><thead><tr><th>structure</th>"
        + header
        + "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></section>"
    )


def _high_confidence_wrong_table(report: Mapping[str, object]) -> str:
    cases = report.get("high_confidence_wrong_cases")
    if not isinstance(cases, Sequence) or not cases:
        return (
            "<section><h2>High-Confidence Wrong Cases</h2>"
            "<p>No high-confidence wrong cases recorded.</p></section>"
        )
    body = []
    for case in cases:
        if not isinstance(case, Mapping):
            continue
        body.append(
            "<tr>"
            f"<td>{_escape(case.get('protein_id', ''))}</td>"
            f"<td>{_escape(case.get('predicted_fold_class', ''))}</td>"
            f"<td>{_escape(case.get('structure_fold_class', ''))}</td>"
            f"<td>{_escape(case.get('confidence', ''))}</td>"
            f"<td>{_escape(case.get('gate_path', ''))}</td>"
            "</tr>"
        )
    return (
        "<section><h2>High-Confidence Wrong Cases</h2>"
        "<table><thead><tr><th>protein</th><th>predicted</th>"
        "<th>structure</th><th>confidence</th><th>gate path</th>"
        "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></section>"
    )


def _gate_table(rows: Sequence[Mapping[str, object]]) -> str:
    body = []
    for row in rows:
        body.append(
            "<tr>"
            f"<td>{_escape(row['protein_id'])}</td>"
            f"<td>{_escape(row['raw_motif_predicted_fold_class'])}</td>"
            f"<td>{_escape(row['predicted_fold_class'])}</td>"
            f"<td>{_escape(row['structure_fold_class'])}</td>"
            f"<td>{_escape(row['gate_path'])}</td>"
            f"<td>{_escape(row['confidence'])}</td>"
            f"<td>{_escape(row['likely_failure_reason'])}</td>"
            "</tr>"
        )
    return (
        "<table><thead><tr><th>protein</th><th>raw motif</th>"
        "<th>hierarchy</th><th>structure</th><th>gate path</th>"
        "<th>confidence</th><th>failure reason</th></tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table>"
    )


def _score_bars(rows: Sequence[Mapping[str, object]]) -> str:
    sections = []
    keys = (
        "disorder_gate_score",
        "compactness_gate_score",
        "segmentation_gate_score",
        "secondary_structure_gate_score",
    )
    for row in rows:
        bars = []
        for key in keys:
            width = round(float(row[key]) * 100, 3)
            bars.append(
                "<div class=\"bar-row\">"
                f"<span>{_escape(key)}</span>"
                "<div class=\"bar-track\">"
                f"<div class=\"bar\" style=\"width:{width}%\"></div>"
                "</div>"
                f"<strong>{_escape(row[key])}</strong>"
                "</div>"
            )
        sections.append(
            "<section class=\"protein\">"
            f"<h3>{_escape(row['protein_id'])}</h3>"
            + "".join(bars)
            + "</section>"
        )
    return "".join(sections)


def render_hierarchical_gate_dashboard(report: Mapping[str, object]) -> str:
    rows = [row for row in report.get("rows", []) if isinstance(row, Mapping)]
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Hierarchical Folding Decision Gates</title>
<style>
body {{ margin: 0; font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #17211c; background: #fff; }}
header {{ padding: 28px clamp(18px, 4vw, 54px); background: #f1f5ee; border-bottom: 1px solid #cad6d0; }}
main {{ padding: 24px clamp(18px, 4vw, 54px) 52px; }}
h1 {{ margin: 0 0 10px; font-size: clamp(28px, 4vw, 44px); letter-spacing: 0; }}
h2 {{ font-size: 22px; letter-spacing: 0; }}
h3 {{ margin: 18px 0 8px; font-size: 16px; }}
p {{ max-width: 980px; line-height: 1.55; }}
.warning {{ display: inline-block; padding: 8px 10px; border: 1px solid #9b3d2e; color: #9b3d2e; background: #fff8f6; font-weight: 700; }}
.metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 10px; margin-top: 18px; }}
.metric {{ border: 1px solid #cad6d0; background: #f8faf7; padding: 12px; min-height: 82px; }}
.metric span {{ display: block; color: #5d6762; font-size: 13px; }}
.metric strong {{ display: block; margin-top: 8px; font-size: 22px; }}
.protein {{ border-top: 1px solid #cad6d0; padding: 12px 0; }}
.bar-row {{ display: grid; grid-template-columns: minmax(220px, 300px) 1fr 70px; gap: 10px; align-items: center; margin: 8px 0; }}
.bar-track {{ height: 16px; border: 1px solid #cad6d0; background: #fff; }}
.bar {{ height: 100%; background: #3c7a42; }}
table {{ width: 100%; border-collapse: collapse; font-size: 14px; }}
th, td {{ border: 1px solid #cad6d0; padding: 8px; text-align: left; vertical-align: top; }}
th {{ background: #f1f5ee; }}
@media (max-width: 760px) {{ .bar-row {{ grid-template-columns: 1fr; }} }}
</style>
</head>
<body>
<header>
  <p class="warning">NOT A FOLDING SOLUTION / HIERARCHICAL GATE REVIEW</p>
  <h1>Hierarchical Folding Decision Gates</h1>
  <p>This page routes motif evidence through disorder, compactness,
  segmentation, secondary-structure, and confidence gates before making any
  broad fold interpretation.</p>
  <div class="metrics">{_metric_cards(report)}</div>
</header>
<main>
  {_mapping_table("10 vs 50 Stability", report.get("result_compared_to_10_row_benchmark"))}
  {_mapping_table("Accuracy By Class", report.get("accuracy_by_class"))}
  {_mapping_table("Abstention Rate By Class", report.get("abstention_rate_by_class"))}
  {_mapping_table("Gate Path Distribution By True Class", report.get("gate_path_distribution_by_true_class"))}
  {_mapping_table("Per-Class Stability Status", report.get("per_class_stability_status"))}
  {_confusion_matrix_table(report)}
  {_high_confidence_wrong_table(report)}
  <section>
    <h2>Gate Paths</h2>
    {_gate_table(rows)}
  </section>
  <section>
    <h2>Gate Scores</h2>
    {_score_bars(rows)}
  </section>
</main>
</body>
</html>
"""


def load_hierarchical_gate_inputs(
    benchmark_file: Path,
    structure_evidence_file: Path,
) -> tuple[tuple[FoldingReferenceExample, ...], tuple[StructureEvidenceRow, ...]]:
    return load_motif_alignment_inputs(benchmark_file, structure_evidence_file)
