from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_order_aware import (
    ORDER_AWARE_BENCHMARK_KIND,
    extract_order_aware_features,
    load_order_aware_inputs,
    order_aware_control_separation_rows,
)
from pharmacotopology.folding_structure_benchmark import (
    LABEL_BENCHMARK_KIND_V0,
    StructureEvidenceRow,
)
from pharmacotopology.folding_topology import (
    CONTACT_PROXY_DIMENSIONS,
    HYDROPHOBIC_AMINO_ACIDS,
    FoldingReferenceExample,
    FoldingTopologySignature,
    contact_map_proxy_similarity,
    normalize_fold_class,
    normalize_sequence,
    sequence_features,
    signature_to_dict,
)


MOTIF_ALIGNMENT_BENCHMARK_KIND = "motif_to_structure_alignment_benchmark"
MOTIF_SIGNATURE_KIND = "order_aware_motif_evidence_vector"
ABSTAINED_CLASS = "insufficient_topology_evidence"


@dataclass(frozen=True)
class MotifEvidenceVector:
    alpha_periodicity_evidence: float
    beta_alternation_evidence: float
    compact_core_evidence: float
    disorder_run_evidence: float
    domain_boundary_evidence: float
    long_range_closure_evidence: float
    breaker_turn_evidence: float
    charge_frustration_evidence: float
    local_alpha_pressure_evidence: float
    local_beta_pressure_evidence: float
    local_disorder_pressure_evidence: float
    mixed_motif_evidence: float

    def to_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class MotifPrediction:
    protein_id: str
    sequence_length: int
    motif_evidence: MotifEvidenceVector
    topology_signature: FoldingTopologySignature
    predicted_fold_class: str
    raw_predicted_fold_class: str
    dominant_features: tuple[str, ...]
    conflicting_features: tuple[str, ...]
    evidence_conflict_score: float
    confidence: float
    uncertainty_radius: float
    claim_strength: str
    forced_prediction: bool
    abstained: bool
    motif_signal_seen: bool


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


def _rounded(value: float) -> float:
    return round(_clamp(value), 6)


def _mean(values: Iterable[float]) -> float:
    values = tuple(values)
    if not values:
        return 0.0
    return sum(values) / len(values)


def _max_run(sequence: str, alphabet: frozenset[str]) -> int:
    current = 0
    best = 0
    for residue in sequence:
        if residue in alphabet:
            current += 1
            best = max(best, current)
        else:
            current = 0
    return best


def _alpha_periodicity(sequence: str) -> float:
    hydrophobic = {index for index, residue in enumerate(sequence) if residue in HYDROPHOBIC_AMINO_ACIDS}
    if len(hydrophobic) < 2:
        return 0.0
    periodic_hits = 0
    reviewed = 0
    for index in range(len(sequence)):
        if index not in hydrophobic:
            continue
        for offset in (3, 4):
            reviewed += 1
            if index + offset in hydrophobic:
                periodic_hits += 1
    return round(periodic_hits / max(reviewed, 1), 6)


def _beta_alternation(sequence: str, window_size: int = 7) -> float:
    if len(sequence) < window_size:
        return 0.0
    scores: list[float] = []
    for start in range(0, len(sequence) - window_size + 1):
        window = sequence[start : start + window_size]
        states = [residue in HYDROPHOBIC_AMINO_ACIDS for residue in window]
        alternations = sum(1 for left, right in zip(states, states[1:]) if left != right)
        scores.append(alternations / (window_size - 1))
    return round(_mean(scores), 6)


def _domain_split_signal(features: Mapping[str, float], sequence_length: int) -> float:
    length_pressure = _clamp((sequence_length - 120) / 160)
    return _rounded(
        float(features.get("segment_boundary_contrast", 0.0)) * 1.55
        + float(features.get("contact_prior_crossing_fraction", 0.0)) * 0.52
        + length_pressure * 0.22
        + float(features.get("contact_prior_density", 0.0)) * 0.12
        - float(features.get("hydrophobic_interruption_rate", 0.0)) * 0.16
    )


def motif_evidence_from_sequence(sequence: str) -> tuple[MotifEvidenceVector, dict[str, float]]:
    normalized = normalize_sequence(sequence)
    order_features, _ = extract_order_aware_features(normalized)
    alpha_periodicity = _alpha_periodicity(normalized)
    beta_alternation = _beta_alternation(normalized)
    disorder_run = _max_run(normalized, frozenset("DEKRPQSNG")) / max(len(normalized), 1)
    compact_core = _rounded(
        float(order_features["local_collapse_pressure_mean"]) * 0.62
        + float(order_features["contact_prior_mean_weight"]) * 1.18
        + float(order_features["contact_prior_max_weight"]) * 0.34
        + float(order_features["hydrophobic_cluster_periodicity"]) * 0.11
        - float(order_features["hydrophobic_interruption_rate"]) * 0.16
    )
    local_alpha = _rounded(
        float(order_features["local_helix_pressure_mean"]) * 0.70
        + alpha_periodicity * 0.42
        + compact_core * 0.22
        - float(order_features["breaker_density"]) * 0.46
        - max(
            float(order_features["local_beta_pressure_mean"])
            - float(order_features["local_helix_pressure_mean"]),
            0.0,
        )
        * 0.28
    )
    local_beta = _rounded(
        float(order_features["local_beta_pressure_mean"]) * 0.76
        + beta_alternation * 0.34
        + float(order_features["long_range_closure_potential"]) * 0.34
        + float(order_features["contact_prior_average_order"]) * 1.10
        - float(order_features["breaker_density"]) * 0.30
    )
    local_disorder = _rounded(
        float(order_features["local_disorder_pressure_mean"]) * 0.62
        + float(order_features["breaker_density"]) * 1.10
        + float(order_features["hydrophobic_interruption_rate"]) * 0.28
        + disorder_run * 1.20
        - compact_core * 0.26
        - float(order_features["contact_prior_mean_weight"]) * 0.26
    )
    domain = _domain_split_signal(order_features, len(normalized))
    long_range = _rounded(
        float(order_features["long_range_closure_potential"]) * 1.55
        + float(order_features["contact_prior_average_order"]) * 2.05
        + float(order_features["contact_prior_crossing_fraction"]) * 0.24
    )
    breaker_turn = _rounded(
        float(order_features["breaker_density"]) * 1.65
        + float(order_features["breaker_boundary_alignment"]) * 0.42
        + float(order_features["hydrophobic_interruption_rate"]) * 0.22
    )
    charge_frustration = _rounded(
        float(order_features["charge_blockiness"]) * 0.28
        + float(order_features["same_charge_repulsion_pressure"]) * 1.75
        - float(order_features["opposite_charge_pairing_potential"]) * 0.55
    )
    mixed = _rounded(
        min(local_alpha, local_beta) * 0.95
        + compact_core * 0.24
        + long_range * 0.18
        - local_disorder * 0.18
        - domain * 0.10
    )
    evidence = MotifEvidenceVector(
        alpha_periodicity_evidence=_rounded(alpha_periodicity),
        beta_alternation_evidence=_rounded(beta_alternation),
        compact_core_evidence=compact_core,
        disorder_run_evidence=_rounded(disorder_run),
        domain_boundary_evidence=domain,
        long_range_closure_evidence=long_range,
        breaker_turn_evidence=breaker_turn,
        charge_frustration_evidence=charge_frustration,
        local_alpha_pressure_evidence=local_alpha,
        local_beta_pressure_evidence=local_beta,
        local_disorder_pressure_evidence=local_disorder,
        mixed_motif_evidence=mixed,
    )
    motif_features = dict(order_features)
    motif_features.update(
        {
            "alpha_periodicity_signal": alpha_periodicity,
            "beta_alternation_signal": beta_alternation,
            "disorder_run_signal": round(disorder_run, 6),
        }
    )
    return evidence, motif_features


def _class_scores(evidence: MotifEvidenceVector) -> dict[str, float]:
    return {
        "alpha_rich": _rounded(
            evidence.local_alpha_pressure_evidence * 0.48
            + evidence.alpha_periodicity_evidence * 0.30
            + evidence.compact_core_evidence * 0.18
            + evidence.long_range_closure_evidence * 0.08
            - evidence.local_disorder_pressure_evidence * 0.16
            - evidence.disorder_run_evidence * 0.18
            - evidence.domain_boundary_evidence * 0.08
        ),
        "beta_rich": _rounded(
            evidence.local_beta_pressure_evidence * 0.42
            + evidence.beta_alternation_evidence * 0.36
            + evidence.long_range_closure_evidence * 0.28
            + evidence.compact_core_evidence * 0.10
            - evidence.local_alpha_pressure_evidence * 0.12
            - evidence.local_disorder_pressure_evidence * 0.12
            - evidence.disorder_run_evidence * 0.14
        ),
        "alpha_beta_mixed": _rounded(
            evidence.mixed_motif_evidence * 0.60
            + min(
                evidence.local_alpha_pressure_evidence,
                evidence.local_beta_pressure_evidence,
            )
            * 0.22
            + evidence.compact_core_evidence * 0.10
            + evidence.long_range_closure_evidence * 0.06
            - evidence.local_disorder_pressure_evidence * 0.10
            - evidence.disorder_run_evidence * 0.14
        ),
        "disordered_flexible": _rounded(
            evidence.local_disorder_pressure_evidence * 0.46
            + evidence.disorder_run_evidence * 0.62
            + evidence.breaker_turn_evidence * 0.18
            + evidence.charge_frustration_evidence * 0.10
            - evidence.compact_core_evidence * 0.20
        ),
        "multidomain_boundary": _rounded(
            evidence.domain_boundary_evidence * 0.82
            + evidence.long_range_closure_evidence * 0.20
            + evidence.charge_frustration_evidence * 0.08
            - evidence.compact_core_evidence * 0.10
        ),
    }


def _topology_signature_from_motifs(
    evidence: MotifEvidenceVector,
    sequence_complexity: float,
    uncertainty_radius: float,
) -> FoldingTopologySignature:
    return FoldingTopologySignature(
        sequence_complexity=sequence_complexity,
        secondary_structure_balance=_rounded(
            1.0
            - abs(
                evidence.local_alpha_pressure_evidence
                - evidence.local_beta_pressure_evidence
            )
            * 0.70
            - evidence.local_disorder_pressure_evidence * 0.18
            - evidence.disorder_run_evidence * 0.10
        ),
        contact_map_closure=_rounded(
            evidence.compact_core_evidence * 0.46
            + evidence.long_range_closure_evidence * 0.34
            + evidence.mixed_motif_evidence * 0.12
            - evidence.local_disorder_pressure_evidence * 0.10
            - evidence.disorder_run_evidence * 0.08
        ),
        hydrophobic_core_closure=evidence.compact_core_evidence,
        loop_disorder_pressure=_rounded(
            evidence.local_disorder_pressure_evidence * 0.48
            + evidence.disorder_run_evidence * 0.32
            + evidence.breaker_turn_evidence * 0.28
            + evidence.charge_frustration_evidence * 0.10
        ),
        domain_boundary_stability=_rounded(
            1.0 - evidence.domain_boundary_evidence * 0.70
            - evidence.charge_frustration_evidence * 0.10
            + evidence.compact_core_evidence * 0.10
        ),
        long_range_contact_order=evidence.long_range_closure_evidence,
        conformational_flexibility=_rounded(
            evidence.local_disorder_pressure_evidence * 0.38
            + evidence.disorder_run_evidence * 0.24
            + evidence.breaker_turn_evidence * 0.30
            + evidence.charge_frustration_evidence * 0.14
            + evidence.domain_boundary_evidence * 0.08
        ),
        knot_or_entanglement_signature=_rounded(
            evidence.domain_boundary_evidence * 0.28
            + evidence.long_range_closure_evidence * 0.38
            + evidence.charge_frustration_evidence * 0.12
        ),
        uncertainty_radius=uncertainty_radius,
    )


def predict_motif_alignment(
    sequence: str,
    *,
    protein_id: str = "sequence",
) -> MotifPrediction:
    normalized = normalize_sequence(sequence)
    sequence_complexity = float(sequence_features(normalized)["sequence_complexity"])
    evidence, _motif_features = motif_evidence_from_sequence(normalized)
    scores = _class_scores(evidence)
    ranked = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    raw_class, top_score = ranked[0]
    second_class, second_score = ranked[1]
    conflict = round(second_score / top_score, 6) if top_score else 1.0
    confidence = _rounded((top_score - second_score) * 1.55 + top_score * 0.20)
    uncertainty = _rounded(0.22 + conflict * 0.42 + (1.0 - confidence) * 0.20)
    claim_strength = "weak"
    if confidence >= 0.58 and conflict <= 0.55:
        claim_strength = "strong"
    elif confidence >= 0.36 and conflict <= 0.72:
        claim_strength = "medium"
    abstained = conflict >= 0.82 or confidence < 0.28
    predicted_class = ABSTAINED_CLASS if abstained else raw_class
    dominant = tuple(
        key
        for key, value in sorted(
            evidence.to_dict().items(),
            key=lambda item: item[1],
            reverse=True,
        )[:3]
        if value >= 0.20
    )
    conflicting = tuple(
        key
        for key, value in sorted(
            evidence.to_dict().items(),
            key=lambda item: item[1],
            reverse=True,
        )[3:7]
        if value >= top_score * 0.72 or value >= 0.35
    )
    signature = _topology_signature_from_motifs(
        evidence,
        sequence_complexity=sequence_complexity,
        uncertainty_radius=uncertainty,
    )
    motif_signal_seen = (
        max(
            evidence.alpha_periodicity_evidence,
            evidence.beta_alternation_evidence,
            evidence.local_alpha_pressure_evidence,
            evidence.local_beta_pressure_evidence,
            evidence.local_disorder_pressure_evidence,
        )
        >= 0.30
        or evidence.compact_core_evidence >= 0.28
    )
    return MotifPrediction(
        protein_id=protein_id,
        sequence_length=len(normalized),
        motif_evidence=evidence,
        topology_signature=signature,
        predicted_fold_class=predicted_class,
        raw_predicted_fold_class=raw_class,
        dominant_features=dominant,
        conflicting_features=conflicting,
        evidence_conflict_score=conflict,
        confidence=confidence,
        uncertainty_radius=uncertainty,
        claim_strength=claim_strength,
        forced_prediction=not abstained,
        abstained=abstained,
        motif_signal_seen=motif_signal_seen,
    )


def likely_failure_reason(
    prediction: MotifPrediction,
    structure_class: str,
) -> str:
    evidence = prediction.motif_evidence
    if prediction.abstained:
        return "evidence_conflict_gated_prediction"
    if prediction.predicted_fold_class == structure_class:
        return ""
    if (
        prediction.predicted_fold_class == "multidomain_boundary"
        and evidence.domain_boundary_evidence >= evidence.compact_core_evidence
    ):
        return "multidomain_signal_confused_with_local_segmentation"
    if (
        prediction.predicted_fold_class == "alpha_beta_mixed"
        and structure_class == "beta_rich"
    ):
        return "beta_like_alternation_under_weighted"
    if (
        prediction.predicted_fold_class == "alpha_beta_mixed"
        and structure_class == "disordered_flexible"
    ):
        return "disorder_pressure_under_weighted"
    if (
        prediction.predicted_fold_class == "alpha_rich"
        and structure_class != "alpha_rich"
    ):
        return "alpha_periodicity_over_weighted"
    if evidence.compact_core_evidence > 0.46 and structure_class == "disordered_flexible":
        return "compact_core_signal_over_weighted"
    if evidence.charge_frustration_evidence > 0.30:
        return "charge_pattern_over_weighted"
    return "motif_evidence_conflict"


def motif_alignment_rows(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
) -> list[dict[str, object]]:
    evidence_by_id = {row.protein_id: row for row in evidence_rows}
    rows: list[dict[str, object]] = []
    for reference in references:
        structure = evidence_by_id[reference.protein_id]
        prediction = predict_motif_alignment(
            reference.sequence,
            protein_id=reference.protein_id,
        )
        label_class = normalize_fold_class(reference.reference_fold_class)
        structure_class = structure.structure_fold_class
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
        wrong_high_confidence = (
            prediction.forced_prediction
            and not class_match
            and prediction.confidence >= 0.58
        )
        evidence_dict = prediction.motif_evidence.to_dict()
        rows.append(
            {
                "protein_id": reference.protein_id,
                "sequence_length": prediction.sequence_length,
                "topology_evidence_vector_kind": MOTIF_SIGNATURE_KIND,
                "motif_signal_seen": prediction.motif_signal_seen,
                "alpha_periodicity_evidence": evidence_dict[
                    "alpha_periodicity_evidence"
                ],
                "beta_alternation_evidence": evidence_dict[
                    "beta_alternation_evidence"
                ],
                "compact_core_evidence": evidence_dict["compact_core_evidence"],
                "disorder_run_evidence": evidence_dict["disorder_run_evidence"],
                "domain_boundary_evidence": evidence_dict[
                    "domain_boundary_evidence"
                ],
                "long_range_closure_evidence": evidence_dict[
                    "long_range_closure_evidence"
                ],
                "breaker_turn_evidence": evidence_dict["breaker_turn_evidence"],
                "charge_frustration_evidence": evidence_dict[
                    "charge_frustration_evidence"
                ],
                "local_alpha_pressure_evidence": evidence_dict[
                    "local_alpha_pressure_evidence"
                ],
                "local_beta_pressure_evidence": evidence_dict[
                    "local_beta_pressure_evidence"
                ],
                "local_disorder_pressure_evidence": evidence_dict[
                    "local_disorder_pressure_evidence"
                ],
                "mixed_motif_evidence": evidence_dict["mixed_motif_evidence"],
                "motif_evidence": prediction.motif_evidence.to_dict(),
                "predicted_topology_signature": signature_to_dict(
                    prediction.topology_signature
                ),
                "dominant_features": ";".join(prediction.dominant_features),
                "conflicting_features": ";".join(prediction.conflicting_features),
                "evidence_conflict_score": prediction.evidence_conflict_score,
                "confidence": prediction.confidence,
                "uncertainty_radius": prediction.uncertainty_radius,
                "claim_strength": prediction.claim_strength,
                "forced_prediction": prediction.forced_prediction,
                "abstained": prediction.abstained,
                "predicted_fold_class": prediction.predicted_fold_class,
                "raw_predicted_fold_class": prediction.raw_predicted_fold_class,
                "structure_fold_class": structure_class,
                "label_fold_class": label_class,
                "prediction_vs_structure_score": similarity,
                "prediction_vs_label_score": label_similarity,
                "prediction_structure_class_match": class_match,
                "prediction_label_class_match": label_match,
                "high_confidence_wrong": wrong_high_confidence,
                "likely_failure_reason": likely_failure_reason(
                    prediction,
                    structure_class,
                ),
                "folding_problem_solved": False,
                "folding_solution_claim_created": False,
            }
        )
    return rows


def failure_diagnosis_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    diagnosis: list[dict[str, object]] = []
    for row in rows:
        if bool(row["prediction_structure_class_match"]):
            continue
        diagnosis.append(
            {
                "protein_id": row["protein_id"],
                "predicted_class": row["predicted_fold_class"],
                "raw_predicted_class": row["raw_predicted_fold_class"],
                "reference_label_class": row["label_fold_class"],
                "external_label_class": row["label_fold_class"],
                "structure_derived_class": row["structure_fold_class"],
                "dominant_prediction_features": row["dominant_features"],
                "conflicting_features": row["conflicting_features"],
                "likely_failure_reason": row["likely_failure_reason"],
                "confidence": row["confidence"],
                "uncertainty_radius": row["uncertainty_radius"],
                "claim_strength": row["claim_strength"],
                "high_confidence_wrong": row["high_confidence_wrong"],
                "abstained": row["abstained"],
            }
        )
    return diagnosis


def evidence_conflict_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    conflicts: list[dict[str, object]] = []
    for row in rows:
        conflicts.append(
            {
                "protein_id": row["protein_id"],
                "predicted_class": row["predicted_fold_class"],
                "raw_predicted_class": row["raw_predicted_fold_class"],
                "structure_derived_class": row["structure_fold_class"],
                "evidence_conflict_score": row["evidence_conflict_score"],
                "confidence": row["confidence"],
                "uncertainty_radius": row["uncertainty_radius"],
                "claim_strength": row["claim_strength"],
                "forced_prediction": row["forced_prediction"],
                "abstained": row["abstained"],
                "dominant_features": row["dominant_features"],
                "conflicting_features": row["conflicting_features"],
            }
        )
    return conflicts


def _mean_from_rows(rows: Sequence[Mapping[str, object]], key: str) -> float:
    return round(_mean(float(row.get(key, 0.0)) for row in rows), 6)


def build_motif_alignment_report(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
    *,
    source_benchmark_file: Path,
    structure_evidence_file: Path,
) -> dict[str, object]:
    rows = motif_alignment_rows(references, evidence_rows)
    controls = order_aware_control_separation_rows(references)
    forced_rows = [row for row in rows if bool(row["forced_prediction"])]
    structure_matches = sum(
        1 for row in rows if bool(row["prediction_structure_class_match"])
    )
    label_matches = sum(1 for row in rows if bool(row["prediction_label_class_match"]))
    high_confidence_wrong = sum(
        1 for row in rows if bool(row["high_confidence_wrong"])
    )
    motif_signal_seen = any(bool(row["motif_signal_seen"]) for row in rows)
    contact_prior_signal_seen = any(
        float(extract_order_aware_features(reference.sequence)[0]["contact_prior_density"])
        > 0
        for reference in references
    )
    conflict_mean = _mean_from_rows(rows, "evidence_conflict_score")
    revision_required = (
        high_confidence_wrong > 0
        or structure_matches < len(rows)
        or any(row["abstained"] for row in rows)
    )
    return {
        "benchmark_kind": MOTIF_ALIGNMENT_BENCHMARK_KIND,
        "topology_evidence_vector_kind": MOTIF_SIGNATURE_KIND,
        "source_order_aware_benchmark_kind": ORDER_AWARE_BENCHMARK_KIND,
        "source_label_benchmark_kind": LABEL_BENCHMARK_KIND_V0,
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
        "sequence_order_sensitivity_score": _mean_from_rows(
            controls,
            "separation_score",
        ),
        "real_vs_shuffled_separation_mean": _mean_from_rows(
            controls,
            "separation_score",
        ),
        "contact_prior_signal_seen": contact_prior_signal_seen,
        "motif_signal_seen": motif_signal_seen,
        "evidence_conflict_mean": conflict_mean,
        "uncertainty_gating_used": True,
        "forced_prediction_count": len(forced_rows),
        "abstained_prediction_count": sum(1 for row in rows if bool(row["abstained"])),
        "high_confidence_wrong_count": high_confidence_wrong,
        "failure_diagnosis_available": True,
        "revision_required": revision_required,
        "claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "boundary_statement": (
            "This motif alignment layer converts sequence-only order features "
            "into explicit motif evidence, then uses uncertainty gating before "
            "class interpretation. Reference labels and structure classes are "
            "used only after prediction for diagnosis."
        ),
        "rows": rows,
    }


def write_motif_alignment_outputs(
    *,
    report: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    failure_rows: Sequence[Mapping[str, object]],
    conflict_rows: Sequence[Mapping[str, object]],
    report_path: Path,
    rows_path: Path,
    failure_path: Path,
    conflicts_path: Path,
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
        nested_json_fields={"motif_evidence", "predicted_topology_signature"},
    )
    _write_csv_rows(failure_rows, failure_path, nested_json_fields=set())
    _write_csv_rows(conflict_rows, conflicts_path, nested_json_fields=set())
    dashboard_path.write_text(render_motif_alignment_dashboard(report), encoding="utf-8")
    return report_path, rows_path, failure_path, conflicts_path, dashboard_path


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


def _evidence_bars(rows: Sequence[Mapping[str, object]]) -> str:
    sections = []
    for row in rows:
        evidence = row.get("motif_evidence", {})
        if not isinstance(evidence, Mapping):
            continue
        bars = []
        for key, value in evidence.items():
            width = round(float(value) * 100, 3)
            bars.append(
                "<div class=\"bar-row\">"
                f"<span>{_escape(key.replace('_evidence', ''))}</span>"
                "<div class=\"bar-track\">"
                f"<div class=\"bar\" style=\"width:{width}%\"></div>"
                "</div>"
                f"<strong>{_escape(value)}</strong>"
                "</div>"
            )
        sections.append(
            "<section class=\"protein\">"
            f"<h3>{_escape(row['protein_id'])}</h3>"
            f"<p>predicted: {_escape(row['predicted_fold_class'])} | "
            f"structure: {_escape(row['structure_fold_class'])} | "
            f"claim: {_escape(row['claim_strength'])}</p>"
            + "".join(bars)
            + "</section>"
        )
    return "".join(sections)


def _failure_table(rows: Sequence[Mapping[str, object]]) -> str:
    failure_rows = failure_diagnosis_rows(rows)
    body = []
    for row in failure_rows:
        body.append(
            "<tr>"
            f"<td>{_escape(row['protein_id'])}</td>"
            f"<td>{_escape(row['predicted_class'])}</td>"
            f"<td>{_escape(row['structure_derived_class'])}</td>"
            f"<td>{_escape(row['likely_failure_reason'])}</td>"
            f"<td>{_escape(row['claim_strength'])}</td>"
            f"<td>{_escape(row['high_confidence_wrong'])}</td>"
            "</tr>"
        )
    return (
        "<table><thead><tr><th>protein_id</th><th>predicted</th>"
        "<th>structure</th><th>failure_reason</th><th>claim</th>"
        "<th>high_conf_wrong</th></tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table>"
    )


def render_motif_alignment_dashboard(report: Mapping[str, object]) -> str:
    rows = [row for row in report.get("rows", []) if isinstance(row, Mapping)]
    metrics = {
        "Motif Signal": report.get("motif_signal_seen", False),
        "Conflict Mean": report.get("evidence_conflict_mean", 0.0),
        "Forced": report.get("forced_prediction_count", 0),
        "Abstained": report.get("abstained_prediction_count", 0),
        "High-Confidence Wrong": report.get("high_confidence_wrong_count", 0),
        "Claim Allowed": report.get("claim_allowed", False),
    }
    metric_cards = "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(value)}</strong>"
        "</div>"
        for label, value in metrics.items()
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Motif-to-Structure Alignment</title>
<style>
body {{ margin: 0; font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #17211c; background: #fff; }}
header {{ padding: 28px clamp(18px, 4vw, 54px); background: #eef5f2; border-bottom: 1px solid #cad6d0; }}
main {{ padding: 24px clamp(18px, 4vw, 54px) 52px; }}
h1 {{ margin: 0 0 10px; font-size: clamp(28px, 4vw, 44px); letter-spacing: 0; }}
h2 {{ font-size: 22px; letter-spacing: 0; }}
h3 {{ margin: 18px 0 8px; font-size: 16px; }}
p {{ max-width: 980px; line-height: 1.55; }}
.warning {{ display: inline-block; padding: 8px 10px; border: 1px solid #9b3d2e; color: #9b3d2e; background: #fff8f6; font-weight: 700; }}
.metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 10px; margin-top: 18px; }}
.metric {{ border: 1px solid #cad6d0; background: #f7faf8; padding: 12px; min-height: 82px; }}
.metric span {{ display: block; color: #5d6762; font-size: 13px; }}
.metric strong {{ display: block; margin-top: 8px; font-size: 22px; }}
.protein {{ border-top: 1px solid #cad6d0; padding: 12px 0; }}
.bar-row {{ display: grid; grid-template-columns: minmax(180px, 260px) 1fr 70px; gap: 10px; align-items: center; margin: 8px 0; }}
.bar-track {{ height: 16px; border: 1px solid #cad6d0; background: #fff; }}
.bar {{ height: 100%; background: #0d766e; }}
table {{ width: 100%; border-collapse: collapse; font-size: 14px; }}
th, td {{ border: 1px solid #cad6d0; padding: 8px; text-align: left; vertical-align: top; }}
th {{ background: #eef5f2; }}
@media (max-width: 760px) {{ .bar-row {{ grid-template-columns: 1fr; }} }}
</style>
</head>
<body>
<header>
  <p class="warning">NOT A FOLDING SOLUTION / MOTIF ALIGNMENT REVIEW</p>
  <h1>Motif-to-Structure Alignment</h1>
  <p>This page shows motif evidence before fold-class interpretation, then
  surfaces conflicts, abstentions, and high-confidence wrong cases.</p>
  <div class="metrics">{metric_cards}</div>
</header>
<main>
  <section>
    <h2>Motif Evidence Per Protein</h2>
    {_evidence_bars(rows)}
  </section>
  <section>
    <h2>Failure Reasons</h2>
    {_failure_table(rows)}
  </section>
</main>
</body>
</html>
"""


def load_motif_alignment_inputs(
    benchmark_file: Path,
    structure_evidence_file: Path,
) -> tuple[tuple[FoldingReferenceExample, ...], tuple[StructureEvidenceRow, ...]]:
    return load_order_aware_inputs(benchmark_file, structure_evidence_file)
