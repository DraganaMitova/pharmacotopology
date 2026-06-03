from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_hierarchical_gates import (
    HIERARCHICAL_GATE_BENCHMARK_KIND,
    HierarchicalGatePrediction,
    hierarchical_gate_rows,
    load_hierarchical_gate_inputs,
    predict_hierarchical_gate,
)
from pharmacotopology.folding_motif_alignment import (
    ABSTAINED_CLASS,
    MOTIF_SIGNATURE_KIND,
    motif_evidence_from_sequence,
)
from pharmacotopology.folding_order_aware import extract_order_aware_features
from pharmacotopology.folding_structure_benchmark import (
    LABEL_BENCHMARK_KIND_V0,
    StructureEvidenceRow,
)
from pharmacotopology.folding_topology import (
    CONTACT_PROXY_DIMENSIONS,
    HYDROPHOBIC_AMINO_ACIDS,
    FoldingReferenceExample,
    contact_map_proxy_similarity,
    normalize_fold_class,
    normalize_sequence,
    sequence_features,
)


REGIME_ANALYSIS_BENCHMARK_KIND = "protein_regime_routing_failure_cohort_analysis"
REGIME_SIGNATURE_KIND = "sequence_only_protein_regime_tag"

PROTEIN_REGIMES: tuple[str, ...] = (
    "compact_single_domain",
    "multidomain_modular",
    "intrinsically_disordered",
    "repeat_like",
    "coiled_coil_or_fibrous",
    "membrane_like",
    "low_complexity_region_rich",
    "small_peptide_or_fragment",
    "ambiguous_regime",
)

RISK_ABSTAIN_REGIMES = frozenset(
    {
        "small_peptide_or_fragment",
        "repeat_like",
        "coiled_coil_or_fibrous",
        "membrane_like",
        "ambiguous_regime",
    }
)


@dataclass(frozen=True)
class ProteinRegimePrediction:
    protein_id: str
    sequence_length: int
    protein_regime: str
    regime_confidence: float
    regime_reason: str
    regime_allowed_gate_path: str
    regime_specific_warning: str
    regime_detection_used: bool
    sequence_complexity: float
    hydrophobic_fraction: float
    hydrophobic_run_fraction: float
    low_complexity_pressure: float
    membrane_pressure: float
    repeat_pressure: float
    compact_single_domain_pressure: float
    multidomain_modular_pressure: float
    intrinsically_disordered_pressure: float
    coiled_coil_or_fibrous_pressure: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class RegimeRoutedPrediction:
    protein_id: str
    sequence_length: int
    regime_prediction: ProteinRegimePrediction
    hierarchy_prediction: HierarchicalGatePrediction
    predicted_fold_class: str
    gate_path: tuple[str, ...]
    gate_decision_reason: str
    confidence: float
    uncertainty_radius: float
    claim_strength: str
    forced_prediction: bool
    abstained: bool
    regime_router_changed_prediction: bool


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


def _rounded(value: float) -> float:
    return round(_clamp(value), 6)


def _mean(values: Iterable[float]) -> float:
    values = tuple(values)
    if not values:
        return 0.0
    return sum(values) / len(values)


def _bool_mean(values: Iterable[bool]) -> float:
    values = tuple(values)
    if not values:
        return 0.0
    return round(sum(1 for value in values if value) / len(values), 6)


def _fraction(sequence: str, alphabet: frozenset[str]) -> float:
    if not sequence:
        return 0.0
    return sum(1 for residue in sequence if residue in alphabet) / len(sequence)


def _max_run_fraction(sequence: str, alphabet: frozenset[str]) -> float:
    current = 0
    best = 0
    for residue in sequence:
        if residue in alphabet:
            current += 1
            best = max(best, current)
        else:
            current = 0
    return round(best / max(len(sequence), 1), 6)


def _repeat_kmer_pressure(sequence: str, kmer_size: int = 3) -> float:
    if len(sequence) < kmer_size * 4:
        return 0.0
    kmers = [
        sequence[index : index + kmer_size]
        for index in range(0, len(sequence) - kmer_size + 1)
    ]
    if not kmers:
        return 0.0
    counts = Counter(kmers)
    repeated = sum(count for count in counts.values() if count > 1)
    return round(repeated / len(kmers), 6)


def _score_gap(scores: Mapping[str, float]) -> float:
    ranked = sorted(scores.values(), reverse=True)
    if len(ranked) < 2:
        return 1.0
    return round(ranked[0] - ranked[1], 6)


def _allowed_gate_path(regime: str) -> str:
    mapping = {
        "compact_single_domain": "secondary_structure_gate_or_abstain",
        "multidomain_modular": "segmentation_gate_or_abstain",
        "intrinsically_disordered": "disorder_gate_or_abstain",
        "repeat_like": "abstain_until_repeat_specific_gate_exists",
        "coiled_coil_or_fibrous": "abstain_until_fibrous_gate_exists",
        "membrane_like": "abstain_until_membrane_gate_exists",
        "low_complexity_region_rich": "disorder_or_low_complexity_abstain",
        "small_peptide_or_fragment": "abstain_until_fragment_context_exists",
        "ambiguous_regime": "abstain_until_regime_conflict_resolved",
    }
    return mapping.get(regime, "abstain_until_regime_conflict_resolved")


def detect_protein_regime(
    sequence: str,
    *,
    protein_id: str = "sequence",
) -> ProteinRegimePrediction:
    normalized = normalize_sequence(sequence)
    length = len(normalized)
    basic_features = sequence_features(normalized)
    motif_evidence, _motif_features = motif_evidence_from_sequence(normalized)
    order_features, _contact_edges = extract_order_aware_features(normalized)
    hydrophobic_fraction = float(basic_features["hydrophobic_fraction"])
    sequence_complexity = float(basic_features["sequence_complexity"])
    proline_glycine = float(basic_features["proline_glycine_fraction"])
    charged_fraction = float(basic_features["charged_fraction"])
    disorder_promoting = float(basic_features["disorder_promoting_fraction"])
    hydrophobic_run = _max_run_fraction(normalized, HYDROPHOBIC_AMINO_ACIDS)
    repeat_kmer = _repeat_kmer_pressure(normalized)
    length_long = _clamp((length - 180) / 260)
    length_fragment = _clamp((90 - length) / 45)

    low_complexity_pressure = _rounded(
        (1.0 - sequence_complexity) * 1.10
        + motif_evidence.disorder_run_evidence * 0.36
        + proline_glycine * 0.34
        + disorder_promoting * 0.16
        + charged_fraction * 0.16
        + float(order_features["charge_blockiness"]) * 0.34
    )
    membrane_pressure = _rounded(
        max(hydrophobic_fraction - 0.40, 0.0) * 2.30
        + float(order_features["hydrophobic_cluster_density"]) * 0.42
        + float(order_features["hydrophobic_cluster_max_fraction"]) * 0.80
        + hydrophobic_run * 1.20
        + (0.12 if length >= 180 else 0.0)
        - motif_evidence.local_disorder_pressure_evidence * 0.22
    )
    repeat_pressure = _rounded(
        repeat_kmer * 0.90
        + float(order_features["hydrophobic_cluster_periodicity"]) * 0.22
        + max(0.78 - sequence_complexity, 0.0) * 0.50
        + float(order_features["segment_boundary_contrast"]) * 0.12
    )
    intrinsically_disordered_pressure = _rounded(
        motif_evidence.local_disorder_pressure_evidence * 0.54
        + motif_evidence.disorder_run_evidence * 0.28
        + motif_evidence.breaker_turn_evidence * 0.18
        + proline_glycine * 0.28
        + disorder_promoting * 0.26
        + max(charged_fraction - 0.24, 0.0) * 0.32
        + max(0.35 - hydrophobic_fraction, 0.0) * 0.16
        + float(order_features["charge_blockiness"]) * 0.12
        + low_complexity_pressure * 0.12
        - motif_evidence.compact_core_evidence * 0.26
    )
    multidomain_modular_pressure = _rounded(
        motif_evidence.domain_boundary_evidence * 0.50
        + length_long * 0.26
        + motif_evidence.long_range_closure_evidence * 0.16
        + float(order_features["segment_boundary_contrast"]) * 0.18
        - motif_evidence.local_disorder_pressure_evidence * 0.14
    )
    coiled_pressure = _rounded(
        motif_evidence.local_alpha_pressure_evidence * 0.34
        + motif_evidence.alpha_periodicity_evidence * 0.30
        + float(order_features["hydrophobic_cluster_periodicity"]) * 0.18
        + (0.10 if length >= 120 else 0.0)
        - motif_evidence.local_beta_pressure_evidence * 0.12
        - motif_evidence.local_disorder_pressure_evidence * 0.22
    )
    compact_single_domain_pressure = _rounded(
        motif_evidence.compact_core_evidence * 0.34
        + motif_evidence.long_range_closure_evidence * 0.20
        + sequence_complexity * 0.18
        + (1.0 - motif_evidence.local_disorder_pressure_evidence) * 0.20
        - motif_evidence.domain_boundary_evidence * 0.16
        - length_long * 0.14
    )

    scores = {
        "compact_single_domain": compact_single_domain_pressure,
        "multidomain_modular": multidomain_modular_pressure,
        "intrinsically_disordered": intrinsically_disordered_pressure,
        "repeat_like": repeat_pressure,
        "coiled_coil_or_fibrous": coiled_pressure,
        "membrane_like": membrane_pressure,
        "low_complexity_region_rich": low_complexity_pressure,
        "small_peptide_or_fragment": length_fragment,
    }
    warning = ""
    reason = ""
    if length < 70:
        regime = "small_peptide_or_fragment"
        confidence = _rounded(0.72 + length_fragment * 0.20)
        reason = "sequence length is below the broad-fold context needed by the current gates"
    elif membrane_pressure >= 0.68 and hydrophobic_fraction >= 0.48:
        regime = "membrane_like"
        confidence = _rounded(0.52 + membrane_pressure * 0.38)
        reason = "hydrophobic enrichment and cluster topology resemble a membrane-like regime"
    elif (
        80 <= length <= 140
        and motif_evidence.local_disorder_pressure_evidence >= 0.56
        and proline_glycine >= 0.16
        and motif_evidence.local_beta_pressure_evidence >= 0.35
    ):
        regime = "ambiguous_regime"
        confidence = _rounded(
            0.58
            + motif_evidence.local_disorder_pressure_evidence * 0.20
            + proline_glycine * 0.16
        )
        warning = "folded_small_domain_may_mimic_disorder_from_sequence"
        reason = (
            "compact small-domain beta pressure coexists with high breaker/disorder pressure"
        )
    elif (
        intrinsically_disordered_pressure >= 0.50
        or (
            disorder_promoting >= 0.50
            and hydrophobic_fraction <= 0.34
            and (
                motif_evidence.breaker_turn_evidence >= 0.40
                or charged_fraction >= 0.28
            )
        )
        or (
            charged_fraction >= 0.36
            and hydrophobic_fraction <= 0.30
            and motif_evidence.domain_boundary_evidence >= 0.50
        )
    ):
        regime = "intrinsically_disordered"
        confidence = _rounded(0.50 + intrinsically_disordered_pressure * 0.40)
        reason = "local disorder, breaker, and disorder-run evidence dominate compact closure"
    elif low_complexity_pressure >= 0.62:
        regime = "low_complexity_region_rich"
        confidence = _rounded(0.52 + low_complexity_pressure * 0.36)
        reason = "sequence complexity and local disorder composition are low-complexity-rich"
    elif multidomain_modular_pressure >= 0.55 and length >= 180:
        regime = "multidomain_modular"
        confidence = _rounded(0.48 + multidomain_modular_pressure * 0.42)
        reason = "length and sequence-derived boundary/closure evidence support a modular regime"
        if (
            length >= 420
            and abs(
                motif_evidence.local_alpha_pressure_evidence
                - motif_evidence.local_beta_pressure_evidence
            )
            <= 0.05
            and motif_evidence.local_alpha_pressure_evidence >= 0.44
            and motif_evidence.local_beta_pressure_evidence >= 0.44
        ):
            warning = "long_modular_secondary_structure_ambiguous"
    elif repeat_pressure >= 0.56:
        regime = "repeat_like"
        confidence = _rounded(0.50 + repeat_pressure * 0.40)
        reason = "repeat and periodicity pressure are high enough to require a repeat-specific gate"
    elif coiled_pressure >= 0.58:
        regime = "coiled_coil_or_fibrous"
        confidence = _rounded(0.50 + coiled_pressure * 0.38)
        reason = "alpha periodicity and hydrophobic periodicity resemble a fibrous/coiled regime"
    else:
        top_regime, top_score = max(scores.items(), key=lambda item: item[1])
        if _score_gap(scores) <= 0.045 and top_score >= 0.48:
            regime = "ambiguous_regime"
            confidence = _rounded(0.44 + top_score * 0.28)
            warning = f"near_tie_between_sequence_regime_scores:{top_regime}"
            reason = "the top sequence-derived regime scores are too close to route safely"
        else:
            regime = "compact_single_domain"
            confidence = _rounded(0.46 + compact_single_domain_pressure * 0.42)
            reason = "compact closure is the least risky broad regime under current gates"

    return ProteinRegimePrediction(
        protein_id=protein_id,
        sequence_length=length,
        protein_regime=regime,
        regime_confidence=confidence,
        regime_reason=reason,
        regime_allowed_gate_path=_allowed_gate_path(regime),
        regime_specific_warning=warning,
        regime_detection_used=True,
        sequence_complexity=round(sequence_complexity, 6),
        hydrophobic_fraction=round(hydrophobic_fraction, 6),
        hydrophobic_run_fraction=hydrophobic_run,
        low_complexity_pressure=low_complexity_pressure,
        membrane_pressure=membrane_pressure,
        repeat_pressure=repeat_pressure,
        compact_single_domain_pressure=compact_single_domain_pressure,
        multidomain_modular_pressure=multidomain_modular_pressure,
        intrinsically_disordered_pressure=intrinsically_disordered_pressure,
        coiled_coil_or_fibrous_pressure=coiled_pressure,
    )


def _claim_strength(confidence: float, abstained: bool) -> str:
    if abstained:
        return "abstained"
    if confidence >= 0.72:
        return "strong"
    if confidence >= 0.58:
        return "medium"
    return "weak"


def _abstained_routed_prediction(
    *,
    hierarchy: HierarchicalGatePrediction,
    regime: ProteinRegimePrediction,
    reason: str,
    gate_tag: str,
) -> RegimeRoutedPrediction:
    return RegimeRoutedPrediction(
        protein_id=hierarchy.protein_id,
        sequence_length=hierarchy.sequence_length,
        regime_prediction=regime,
        hierarchy_prediction=hierarchy,
        predicted_fold_class=ABSTAINED_CLASS,
        gate_path=tuple(hierarchy.gate_path) + (gate_tag,),
        gate_decision_reason=reason,
        confidence=0.0,
        uncertainty_radius=_rounded(
            max(hierarchy.uncertainty_radius, 0.58)
            + (1.0 - regime.regime_confidence) * 0.10
        ),
        claim_strength="abstained",
        forced_prediction=False,
        abstained=True,
        regime_router_changed_prediction=not hierarchy.abstained,
    )


def _folded_domain_mimic_disorder_conflict(
    sequence: str,
    *,
    hierarchy: HierarchicalGatePrediction,
    regime: ProteinRegimePrediction,
) -> bool:
    features = sequence_features(sequence)
    evidence = hierarchy.motif_evidence
    return (
        80 <= hierarchy.sequence_length <= 180
        and float(features["sequence_complexity"]) >= 0.90
        and float(features["aromatic_fraction"]) >= 0.08
        and float(features["hydrophobic_fraction"]) >= 0.34
        and regime.compact_single_domain_pressure >= 0.38
        and evidence.compact_core_evidence >= 0.38
        and evidence.long_range_closure_evidence >= 0.34
        and evidence.disorder_run_evidence <= 0.07
        and evidence.local_disorder_pressure_evidence >= 0.50
    )


def predict_regime_routed_gate(
    sequence: str,
    *,
    protein_id: str = "sequence",
) -> RegimeRoutedPrediction:
    regime = detect_protein_regime(sequence, protein_id=protein_id)
    hierarchy = predict_hierarchical_gate(sequence, protein_id=protein_id)
    gate_path = tuple(hierarchy.gate_path)

    if hierarchy.abstained:
        return RegimeRoutedPrediction(
            protein_id=protein_id,
            sequence_length=hierarchy.sequence_length,
            regime_prediction=regime,
            hierarchy_prediction=hierarchy,
            predicted_fold_class=hierarchy.predicted_fold_class,
            gate_path=gate_path + ("regime_router:kept_hierarchy_abstention",),
            gate_decision_reason=(
                f"{hierarchy.gate_decision_reason}; regime router kept abstention "
                f"for {regime.protein_regime}"
            ),
            confidence=hierarchy.confidence,
            uncertainty_radius=hierarchy.uncertainty_radius,
            claim_strength=hierarchy.claim_strength,
            forced_prediction=False,
            abstained=True,
            regime_router_changed_prediction=False,
        )

    if regime.protein_regime in RISK_ABSTAIN_REGIMES:
        return _abstained_routed_prediction(
            hierarchy=hierarchy,
            regime=regime,
            reason=(
                f"{regime.regime_reason}; regime router abstained before accepting "
                "a broad fold-class claim"
            ),
            gate_tag=f"regime_router:abstained_{regime.protein_regime}",
        )

    if regime.protein_regime == "intrinsically_disordered":
        if hierarchy.predicted_fold_class == "disordered_flexible":
            if _folded_domain_mimic_disorder_conflict(
                sequence,
                hierarchy=hierarchy,
                regime=regime,
            ):
                return _abstained_routed_prediction(
                    hierarchy=hierarchy,
                    regime=regime,
                    reason=(
                        "folded-domain mimic disorder conflict: high disorder "
                        "pressure is present, but the sequence also has high "
                        "complexity, aromatic enrichment, compact-core closure, "
                        "long-range closure, and low true disorder-run evidence"
                    ),
                    gate_tag=(
                        "regime_router:"
                        "abstained_folded_domain_mimic_disorder_conflict"
                    ),
                )
            return RegimeRoutedPrediction(
                protein_id=protein_id,
                sequence_length=hierarchy.sequence_length,
                regime_prediction=regime,
                hierarchy_prediction=hierarchy,
                predicted_fold_class=hierarchy.predicted_fold_class,
                gate_path=gate_path + ("regime_router:accepted_disorder_gate",),
                gate_decision_reason=(
                    f"{hierarchy.gate_decision_reason}; sequence regime agrees "
                    "with disorder/flexibility routing"
                ),
                confidence=hierarchy.confidence,
                uncertainty_radius=hierarchy.uncertainty_radius,
                claim_strength=hierarchy.claim_strength,
                forced_prediction=True,
                abstained=False,
                regime_router_changed_prediction=False,
            )
        return _abstained_routed_prediction(
            hierarchy=hierarchy,
            regime=regime,
            reason=(
                "intrinsically disordered regime conflicts with a folded "
                "hierarchy class"
            ),
            gate_tag="regime_router:abstained_disorder_fold_conflict",
        )

    if regime.protein_regime == "low_complexity_region_rich":
        if hierarchy.predicted_fold_class == "disordered_flexible":
            return RegimeRoutedPrediction(
                protein_id=protein_id,
                sequence_length=hierarchy.sequence_length,
                regime_prediction=regime,
                hierarchy_prediction=hierarchy,
                predicted_fold_class=hierarchy.predicted_fold_class,
                gate_path=gate_path + ("regime_router:accepted_low_complexity_disorder",),
                gate_decision_reason=(
                    f"{hierarchy.gate_decision_reason}; low-complexity regime "
                    "does not contradict a disorder/flexibility class"
                ),
                confidence=hierarchy.confidence,
                uncertainty_radius=hierarchy.uncertainty_radius,
                claim_strength=hierarchy.claim_strength,
                forced_prediction=True,
                abstained=False,
                regime_router_changed_prediction=False,
            )
        return _abstained_routed_prediction(
            hierarchy=hierarchy,
            regime=regime,
            reason="low-complexity regime lacks a safe folded-class route",
            gate_tag="regime_router:abstained_low_complexity_fold_conflict",
        )

    if regime.protein_regime == "multidomain_modular":
        if regime.regime_specific_warning:
            return _abstained_routed_prediction(
                hierarchy=hierarchy,
                regime=regime,
                reason=(
                    f"{regime.regime_specific_warning}; long modular sequence "
                    "has ambiguous secondary-structure pressure"
                ),
                gate_tag="regime_router:abstained_long_modular_ambiguity",
            )
        if hierarchy.predicted_fold_class == "multidomain_boundary":
            return RegimeRoutedPrediction(
                protein_id=protein_id,
                sequence_length=hierarchy.sequence_length,
                regime_prediction=regime,
                hierarchy_prediction=hierarchy,
                predicted_fold_class=hierarchy.predicted_fold_class,
                gate_path=gate_path + ("regime_router:accepted_multidomain_gate",),
                gate_decision_reason=(
                    f"{hierarchy.gate_decision_reason}; sequence regime agrees "
                    "with modular/domain routing"
                ),
                confidence=hierarchy.confidence,
                uncertainty_radius=hierarchy.uncertainty_radius,
                claim_strength=hierarchy.claim_strength,
                forced_prediction=True,
                abstained=False,
                regime_router_changed_prediction=False,
            )
        return _abstained_routed_prediction(
            hierarchy=hierarchy,
            regime=regime,
            reason="multidomain regime conflicts with a non-segmentation class",
            gate_tag="regime_router:abstained_multidomain_conflict",
        )

    if regime.protein_regime == "compact_single_domain":
        if hierarchy.predicted_fold_class in {
            "alpha_rich",
            "beta_rich",
            "alpha_beta_mixed",
        }:
            return RegimeRoutedPrediction(
                protein_id=protein_id,
                sequence_length=hierarchy.sequence_length,
                regime_prediction=regime,
                hierarchy_prediction=hierarchy,
                predicted_fold_class=hierarchy.predicted_fold_class,
                gate_path=gate_path + ("regime_router:accepted_compact_single_domain",),
                gate_decision_reason=(
                    f"{hierarchy.gate_decision_reason}; compact-single-domain "
                    "regime accepts secondary-structure gate output"
                ),
                confidence=hierarchy.confidence,
                uncertainty_radius=hierarchy.uncertainty_radius,
                claim_strength=_claim_strength(hierarchy.confidence, False),
                forced_prediction=True,
                abstained=False,
                regime_router_changed_prediction=False,
            )
        return _abstained_routed_prediction(
            hierarchy=hierarchy,
            regime=regime,
            reason="compact-single-domain regime conflicts with boundary/disorder class",
            gate_tag="regime_router:abstained_compact_domain_conflict",
        )

    return _abstained_routed_prediction(
        hierarchy=hierarchy,
        regime=regime,
        reason="unhandled regime route abstained by default",
        gate_tag="regime_router:abstained_unhandled_regime",
    )


def _reference_regime(
    *,
    reference: FoldingReferenceExample,
    structure: StructureEvidenceRow,
) -> str:
    sequence = normalize_sequence(reference.sequence)
    features = sequence_features(sequence)
    hydrophobic = float(features["hydrophobic_fraction"])
    complexity = float(features["sequence_complexity"])
    proline_glycine = float(features["proline_glycine_fraction"])
    structure_features = structure.structure_features
    label_class = normalize_fold_class(reference.reference_fold_class)
    structure_class = structure.structure_fold_class
    if len(sequence) < 70:
        return "small_peptide_or_fragment"
    if hydrophobic >= 0.50 and len(sequence) >= 180:
        return "membrane_like"
    if structure.evidence_kind == "disorder_reference" or structure_class == "disordered_flexible":
        return "intrinsically_disordered"
    if complexity < 0.72 or proline_glycine >= 0.28:
        return "low_complexity_region_rich"
    if (
        label_class == "multidomain_boundary"
        or (
            float(structure_features.get("domain_boundary_signal", 0.0)) >= 0.55
            and float(structure_features.get("residue_count", len(sequence))) >= 180
        )
    ):
        return "multidomain_modular"
    return "compact_single_domain"


def _structure_source_quality(structure: StructureEvidenceRow) -> str:
    if structure.evidence_kind == "disorder_reference":
        return "curated_disorder_reference"
    residue_count = float(structure.structure_features.get("residue_count", 0.0))
    contact_density = float(structure.structure_features.get("contact_density", 0.0))
    if residue_count < 70:
        return "short_coordinate_fragment"
    if contact_density < 0.015:
        return "sparse_long_range_contact_graph"
    return "coordinate_contact_graph"


def _bucket(value: float, *, low: float, high: float) -> str:
    if value < low:
        return "low"
    if value >= high:
        return "high"
    return "mid"


def _length_bucket(length: int) -> str:
    if length < 70:
        return "fragment_lt70"
    if length < 140:
        return "small_70_139"
    if length < 260:
        return "medium_140_259"
    if length < 420:
        return "large_260_419"
    return "very_large_420_plus"


def _manual_review_reasons(row: Mapping[str, object]) -> tuple[str, ...]:
    reasons = []
    if bool(row["structure_label_disagreement"]):
        reasons.append("structure_label_disagreement")
    if bool(row["ambiguous_reference"]):
        reasons.append("ambiguous_reference")
    if bool(row["possible_bad_row"]):
        reasons.append("possible_bad_row")
    if bool(row["high_confidence_wrong"]):
        reasons.append("high_confidence_wrong")
    if str(row["regime_specific_warning"]):
        reasons.append(str(row["regime_specific_warning"]))
    return tuple(reasons)


def _failure_cohort(row: Mapping[str, object]) -> str:
    if bool(row["prediction_structure_class_match"]) and bool(row["forced_prediction"]):
        return "correct_forced"
    if (
        bool(row["hierarchical_high_confidence_wrong"])
        and bool(row["abstained"])
        and bool(row["regime_router_changed_prediction"])
    ):
        return "regime_prevented_high_confidence_wrong"
    if bool(row["high_confidence_wrong"]):
        return "high_confidence_wrong_after_regime"
    if bool(row["forced_prediction"]):
        return "wrong_forced_low_confidence"
    if bool(row["structure_label_disagreement"]):
        return "abstained_on_structure_label_disagreement"
    if str(row["protein_regime"]) in RISK_ABSTAIN_REGIMES:
        return "abstained_on_regime_risk"
    return "abstained_unresolved"


def regime_analysis_rows(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
) -> list[dict[str, object]]:
    evidence_by_id = {row.protein_id: row for row in evidence_rows}
    hierarchy_by_id = {
        str(row["protein_id"]): row
        for row in hierarchical_gate_rows(references, evidence_rows)
    }
    rows: list[dict[str, object]] = []
    for reference in references:
        structure = evidence_by_id[reference.protein_id]
        routed = predict_regime_routed_gate(
            reference.sequence,
            protein_id=reference.protein_id,
        )
        hierarchy_row = hierarchy_by_id[reference.protein_id]
        regime = routed.regime_prediction
        label_class = normalize_fold_class(reference.reference_fold_class)
        structure_class = structure.structure_fold_class
        class_match = routed.predicted_fold_class == structure_class
        label_match = routed.predicted_fold_class == label_class
        similarity = contact_map_proxy_similarity(
            routed.hierarchy_prediction.topology_signature,
            structure.structure_topology_signature,
            CONTACT_PROXY_DIMENSIONS,
        )
        label_similarity = contact_map_proxy_similarity(
            routed.hierarchy_prediction.topology_signature,
            reference.reference_topology_signature,
            CONTACT_PROXY_DIMENSIONS,
        )
        high_confidence_wrong = (
            routed.forced_prediction and not class_match and routed.confidence >= 0.58
        )
        structure_label_disagreement = structure_class != label_class
        reference_regime = _reference_regime(reference=reference, structure=structure)
        ambiguous_reference = structure_label_disagreement or (
            structure.evidence_kind == "coordinate_contact_graph"
            and float(structure.structure_features.get("residue_count", 0.0)) < 70
        )
        possible_bad_row = (
            structure_label_disagreement
            and routed.regime_prediction.protein_regime in {
                "small_peptide_or_fragment",
                "membrane_like",
                "multidomain_modular",
            }
        )
        structure_features = structure.structure_features
        row: dict[str, object] = {
            "protein_id": reference.protein_id,
            "sequence_length": routed.sequence_length,
            "topology_evidence_vector_kind": MOTIF_SIGNATURE_KIND,
            "protein_regime_signature_kind": REGIME_SIGNATURE_KIND,
            "regime_detection_used": regime.regime_detection_used,
            "protein_regime": regime.protein_regime,
            "reference_regime": reference_regime,
            "regime_accuracy_match": regime.protein_regime == reference_regime,
            "regime_confidence": regime.regime_confidence,
            "regime_reason": regime.regime_reason,
            "regime_allowed_gate_path": regime.regime_allowed_gate_path,
            "regime_specific_warning": regime.regime_specific_warning,
            "sequence_complexity": regime.sequence_complexity,
            "hydrophobic_fraction": regime.hydrophobic_fraction,
            "hydrophobic_run_fraction": regime.hydrophobic_run_fraction,
            "low_complexity_pressure": regime.low_complexity_pressure,
            "membrane_pressure": regime.membrane_pressure,
            "repeat_pressure": regime.repeat_pressure,
            "compact_single_domain_pressure": regime.compact_single_domain_pressure,
            "multidomain_modular_pressure": regime.multidomain_modular_pressure,
            "intrinsically_disordered_pressure": (
                regime.intrinsically_disordered_pressure
            ),
            "coiled_coil_or_fibrous_pressure": (
                regime.coiled_coil_or_fibrous_pressure
            ),
            "hierarchical_predicted_fold_class": hierarchy_row["predicted_fold_class"],
            "hierarchical_gate_path": hierarchy_row["gate_path"],
            "hierarchical_confidence": hierarchy_row["confidence"],
            "hierarchical_abstained": hierarchy_row["abstained"],
            "hierarchical_high_confidence_wrong": hierarchy_row[
                "high_confidence_wrong"
            ],
            "regime_router_changed_prediction": (
                routed.regime_router_changed_prediction
            ),
            "predicted_fold_class": routed.predicted_fold_class,
            "structure_fold_class": structure_class,
            "label_fold_class": label_class,
            "prediction_vs_structure_score": similarity,
            "prediction_vs_label_score": label_similarity,
            "prediction_structure_class_match": class_match,
            "prediction_label_class_match": label_match,
            "confidence": routed.confidence,
            "uncertainty_radius": routed.uncertainty_radius,
            "claim_strength": routed.claim_strength,
            "forced_prediction": routed.forced_prediction,
            "abstained": routed.abstained,
            "high_confidence_wrong": high_confidence_wrong,
            "gate_path": " | ".join(routed.gate_path),
            "gate_decision_reason": routed.gate_decision_reason,
            "structure_label_disagreement": structure_label_disagreement,
            "ambiguous_reference": ambiguous_reference,
            "possible_bad_row": possible_bad_row,
            "structure_source_quality": _structure_source_quality(structure),
            "length_bucket": _length_bucket(routed.sequence_length),
            "compactness_bucket": _bucket(
                float(hierarchy_row["compactness_gate_score"]),
                low=0.42,
                high=0.56,
            ),
            "disorder_bucket": _bucket(
                float(hierarchy_row["disorder_gate_score"]),
                low=0.24,
                high=0.42,
            ),
            "segmentation_bucket": _bucket(
                float(hierarchy_row["segmentation_gate_score"]),
                low=0.30,
                high=0.52,
            ),
            "contact_density_bucket": _bucket(
                float(structure_features.get("contact_density", 0.0)),
                low=0.02,
                high=0.07,
            ),
            "long_range_contact_bucket": _bucket(
                float(structure_features.get("long_range_contact_fraction", 0.0)),
                low=0.18,
                high=0.42,
            ),
            "folding_problem_solved": False,
            "folding_solution_claim_created": False,
        }
        manual_reasons = _manual_review_reasons(row)
        row["manual_review_required"] = bool(manual_reasons)
        row["manual_review_reasons"] = ";".join(manual_reasons)
        row["failure_cohort"] = _failure_cohort(row)
        rows.append(row)
    return rows


def _rate_by_key(
    rows: Sequence[Mapping[str, object]],
    *,
    key: str,
    predicate_key: str,
) -> dict[str, float]:
    totals = Counter(str(row[key]) for row in rows)
    positives = Counter(
        str(row[key]) for row in rows if bool(row.get(predicate_key, False))
    )
    return {
        value: round(positives[value] / total, 6)
        for value, total in sorted(totals.items())
    }


def _count_by_key(
    rows: Sequence[Mapping[str, object]],
    *,
    key: str,
    predicate_key: str,
) -> dict[str, int]:
    counts = Counter(
        str(row[key]) for row in rows if bool(row.get(predicate_key, False))
    )
    return {value: counts[value] for value in sorted({str(row[key]) for row in rows})}


def _gate_path_by_regime(
    rows: Sequence[Mapping[str, object]],
) -> dict[str, dict[str, int]]:
    distribution: dict[str, Counter[str]] = {}
    for row in rows:
        regime = str(row["protein_regime"])
        path = str(row["gate_path"]).split(" | ")[-1]
        distribution.setdefault(regime, Counter())[path] += 1
    return {
        regime: dict(sorted(counter.items()))
        for regime, counter in sorted(distribution.items())
    }


def failure_cohort_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, ...], dict[str, object]] = {}
    fields = (
        "failure_cohort",
        "protein_regime",
        "predicted_fold_class",
        "structure_fold_class",
        "label_fold_class",
        "length_bucket",
        "compactness_bucket",
        "disorder_bucket",
        "segmentation_bucket",
        "contact_density_bucket",
        "long_range_contact_bucket",
        "structure_source_quality",
    )
    for row in rows:
        key = tuple(str(row[field]) for field in fields)
        if key not in grouped:
            grouped[key] = {field: row[field] for field in fields}
            grouped[key].update(
                {
                    "count": 0,
                    "forced_prediction_count": 0,
                    "abstained_prediction_count": 0,
                    "correct_count": 0,
                    "high_confidence_wrong_count": 0,
                    "structure_label_disagreement_count": 0,
                    "possible_bad_rows_count": 0,
                    "manual_review_count": 0,
                    "example_proteins": [],
                }
            )
        item = grouped[key]
        item["count"] = int(item["count"]) + 1
        item["forced_prediction_count"] = int(item["forced_prediction_count"]) + int(
            bool(row["forced_prediction"])
        )
        item["abstained_prediction_count"] = int(
            item["abstained_prediction_count"]
        ) + int(bool(row["abstained"]))
        item["correct_count"] = int(item["correct_count"]) + int(
            bool(row["prediction_structure_class_match"])
        )
        item["high_confidence_wrong_count"] = int(
            item["high_confidence_wrong_count"]
        ) + int(bool(row["high_confidence_wrong"]))
        item["structure_label_disagreement_count"] = int(
            item["structure_label_disagreement_count"]
        ) + int(bool(row["structure_label_disagreement"]))
        item["possible_bad_rows_count"] = int(item["possible_bad_rows_count"]) + int(
            bool(row["possible_bad_row"])
        )
        item["manual_review_count"] = int(item["manual_review_count"]) + int(
            bool(row["manual_review_required"])
        )
        examples = item["example_proteins"]
        if isinstance(examples, list) and len(examples) < 4:
            examples.append(row["protein_id"])

    output = []
    for item in grouped.values():
        output.append(
            {
                **item,
                "accuracy": round(
                    int(item["correct_count"]) / max(int(item["count"]), 1),
                    6,
                ),
                "example_proteins": ";".join(str(value) for value in item["example_proteins"]),
            }
        )
    return sorted(
        output,
        key=lambda item: (
            -int(item["count"]),
            str(item["failure_cohort"]),
            str(item["protein_regime"]),
        ),
    )


def high_confidence_wrong_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    keys = (
        "protein_id",
        "protein_regime",
        "reference_regime",
        "predicted_fold_class",
        "structure_fold_class",
        "label_fold_class",
        "confidence",
        "gate_path",
        "failure_cohort",
        "manual_review_reasons",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if bool(row["high_confidence_wrong"])
    ]


def abstention_analysis_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    keys = (
        "protein_id",
        "protein_regime",
        "reference_regime",
        "hierarchical_predicted_fold_class",
        "predicted_fold_class",
        "structure_fold_class",
        "label_fold_class",
        "regime_router_changed_prediction",
        "hierarchical_high_confidence_wrong",
        "structure_label_disagreement",
        "possible_bad_row",
        "gate_decision_reason",
        "failure_cohort",
        "manual_review_reasons",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if bool(row["abstained"])
    ]


def _dominant_failure_cohort(rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    candidates = [row for row in rows if row.get("failure_cohort") != "correct_forced"]
    if not candidates:
        return {
            "failure_cohort": "none",
            "count": 0,
            "protein_regime": "",
            "interpretation": "no failure cohort dominated this run",
        }
    cohort_counts = Counter(str(row["failure_cohort"]) for row in candidates)
    dominant_name, dominant_count = max(
        cohort_counts.items(),
        key=lambda item: (item[1], item[0]),
    )
    dominant_rows = [
        row for row in candidates if str(row["failure_cohort"]) == dominant_name
    ]
    regimes = Counter(str(row["protein_regime"]) for row in dominant_rows)
    examples = [str(row["protein_id"]) for row in dominant_rows[:6]]
    return {
        "failure_cohort": dominant_name,
        "count": dominant_count,
        "regime_distribution": dict(sorted(regimes.items())),
        "manual_review_count": sum(
            1 for row in dominant_rows if bool(row["manual_review_required"])
        ),
        "example_proteins": ";".join(examples),
        "interpretation": (
            "largest non-correct cohort after sequence-only regime routing"
        ),
    }


def _new_failure_modes(rows: Sequence[Mapping[str, object]]) -> tuple[str, ...]:
    modes = []
    if any(bool(row["regime_router_changed_prediction"]) for row in rows):
        modes.append("sequence_regime_routing_changes_fold_claim_surface")
    if any(str(row["protein_regime"]) == "membrane_like" for row in rows):
        modes.append("membrane_like_broad_fold_aliasing_risk")
    if any(str(row["protein_regime"]) == "ambiguous_regime" for row in rows):
        modes.append("folded_small_domain_disorder_mimicry_risk")
    if any(bool(row["structure_label_disagreement"]) for row in rows):
        modes.append("structure_label_disagreement_requires_manual_review")
    if any(
        str(row["regime_specific_warning"]) == "long_modular_secondary_structure_ambiguous"
        for row in rows
    ):
        modes.append("long_modular_secondary_structure_ambiguity")
    return tuple(modes)


def _old_failure_modes_still_closed(
    rows: Sequence[Mapping[str, object]],
) -> dict[str, dict[str, object]]:
    counters = {
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
        "flexible_segmentation_false_multidomain_count": sum(
            1
            for row in rows
            if row["structure_fold_class"] == "disordered_flexible"
            and row["predicted_fold_class"] == "multidomain_boundary"
        ),
    }
    return {
        key: {
            "count": value,
            "status": "closed" if value == 0 else "reopened",
        }
        for key, value in counters.items()
    }


def _manual_review_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    keys = (
        "protein_id",
        "protein_regime",
        "reference_regime",
        "predicted_fold_class",
        "structure_fold_class",
        "label_fold_class",
        "failure_cohort",
        "manual_review_reasons",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if bool(row["manual_review_required"])
    ]


def build_regime_analysis_report(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
    *,
    source_benchmark_file: Path,
    structure_evidence_file: Path,
) -> dict[str, object]:
    rows = regime_analysis_rows(references, evidence_rows)
    cohorts = failure_cohort_rows(rows)
    forced_rows = [row for row in rows if bool(row["forced_prediction"])]
    structure_matches = sum(
        1 for row in rows if bool(row["prediction_structure_class_match"])
    )
    label_matches = sum(1 for row in rows if bool(row["prediction_label_class_match"]))
    high_confidence_wrong_count = sum(
        1 for row in rows if bool(row["high_confidence_wrong"])
    )
    structure_label_disagreement_count = sum(
        1 for row in rows if bool(row["structure_label_disagreement"])
    )
    ambiguous_reference_count = sum(
        1 for row in rows if bool(row["ambiguous_reference"])
    )
    possible_bad_rows_count = sum(1 for row in rows if bool(row["possible_bad_row"]))
    old_modes = _old_failure_modes_still_closed(rows)
    return {
        "benchmark_kind": REGIME_ANALYSIS_BENCHMARK_KIND,
        "source_hierarchical_gate_benchmark_kind": HIERARCHICAL_GATE_BENCHMARK_KIND,
        "source_label_benchmark_kind": LABEL_BENCHMARK_KIND_V0,
        "topology_evidence_vector_kind": MOTIF_SIGNATURE_KIND,
        "protein_regime_signature_kind": REGIME_SIGNATURE_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "structure_evidence_file": str(structure_evidence_file),
        "predictor_input_boundary": "sequence_only_no_labels_no_structure_answers",
        "regime_detection_boundary": "sequence_only_no_labels_no_pdb_no_cath_no_disprot_truth",
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
        "high_confidence_wrong_count": high_confidence_wrong_count,
        "regime_detection_used": all(bool(row["regime_detection_used"]) for row in rows),
        "regime_accuracy": _bool_mean(
            bool(row["regime_accuracy_match"]) for row in rows
        ),
        "regime_confidence_mean": round(
            _mean(float(row["regime_confidence"]) for row in rows),
            6,
        ),
        "accuracy_by_regime": _rate_by_key(
            rows,
            key="protein_regime",
            predicate_key="prediction_structure_class_match",
        ),
        "abstention_by_regime": _rate_by_key(
            rows,
            key="protein_regime",
            predicate_key="abstained",
        ),
        "high_confidence_wrong_by_regime": _count_by_key(
            rows,
            key="protein_regime",
            predicate_key="high_confidence_wrong",
        ),
        "gate_path_by_regime": _gate_path_by_regime(rows),
        "dominant_failure_cohort": _dominant_failure_cohort(rows),
        "new_failure_modes_detected": list(_new_failure_modes(rows)),
        "old_failure_modes_still_closed": old_modes,
        "structure_label_disagreement_count": structure_label_disagreement_count,
        "ambiguous_reference_count": ambiguous_reference_count,
        "possible_bad_rows_count": possible_bad_rows_count,
        "rows_requiring_manual_review": _manual_review_rows(rows),
        "manual_review_row_count": sum(
            1 for row in rows if bool(row["manual_review_required"])
        ),
        "regime_router_changed_prediction_count": sum(
            1 for row in rows if bool(row["regime_router_changed_prediction"])
        ),
        "hierarchical_high_confidence_wrong_count_before_regime_routing": sum(
            1 for row in rows if bool(row["hierarchical_high_confidence_wrong"])
        ),
        "hierarchical_high_confidence_wrong_prevented_by_regime_routing": sum(
            1
            for row in rows
            if bool(row["hierarchical_high_confidence_wrong"])
            and bool(row["regime_router_changed_prediction"])
            and bool(row["abstained"])
        ),
        "failure_cohort_distribution": dict(
            sorted(Counter(str(row["failure_cohort"]) for row in rows).items())
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
            "This layer detects sequence-derived protein regimes before accepting "
            "hierarchical fold gates. Structure classes, labels, accessions, PDB "
            "coordinates, CATH/DisProt truth channels, and reference regimes are "
            "used only after prediction for scoring and failure-cohort analysis."
        ),
        "rows": rows,
    }


def write_regime_analysis_outputs(
    *,
    report: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    cohort_rows: Sequence[Mapping[str, object]],
    high_confidence_wrong: Sequence[Mapping[str, object]],
    abstention_rows: Sequence[Mapping[str, object]],
    report_path: Path,
    rows_path: Path,
    failure_cohorts_path: Path,
    high_confidence_wrong_path: Path,
    abstention_analysis_path: Path,
    dashboard_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(rows, rows_path)
    _write_csv_rows(cohort_rows, failure_cohorts_path)
    _write_csv_rows(
        high_confidence_wrong,
        high_confidence_wrong_path,
        fieldnames=[
            "protein_id",
            "protein_regime",
            "reference_regime",
            "predicted_fold_class",
            "structure_fold_class",
            "label_fold_class",
            "confidence",
            "gate_path",
            "failure_cohort",
            "manual_review_reasons",
        ],
    )
    _write_csv_rows(abstention_rows, abstention_analysis_path)
    dashboard_path.write_text(
        render_regime_dashboard(report),
        encoding="utf-8",
    )
    return (
        report_path,
        rows_path,
        failure_cohorts_path,
        high_confidence_wrong_path,
        abstention_analysis_path,
        dashboard_path,
    )


def _write_csv_rows(
    rows: Sequence[Mapping[str, object]],
    path: Path,
    *,
    fieldnames: Sequence[str] | None = None,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as file:
        names = list(fieldnames or (list(rows[0]) if rows else []))
        if not names:
            return path
        writer = csv.DictWriter(file, fieldnames=names, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in names})
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
        "prediction_vs_label_accuracy",
        "forced_prediction_count",
        "abstained_prediction_count",
        "high_confidence_wrong_count",
        "regime_accuracy",
        "regime_confidence_mean",
        "manual_review_row_count",
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
            rendered = ", ".join(
                f"{_escape(nested_key)}: {_escape(nested_value)}"
                for nested_key, nested_value in value.items()
            )
        else:
            rendered = _escape(value)
        rows.append(
            "<tr>"
            f"<td>{_escape(key)}</td>"
            f"<td>{rendered}</td>"
            "</tr>"
        )
    return (
        f"<section><h2>{_escape(title)}</h2>"
        "<table><thead><tr><th>key</th><th>value</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></section>"
    )


def _heatmap_table(report: Mapping[str, object]) -> str:
    rows = report.get("rows", [])
    if not isinstance(rows, Sequence):
        return ""
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    regimes = sorted({str(row["protein_regime"]) for row in rows if isinstance(row, Mapping)})
    cohorts = sorted({str(row["failure_cohort"]) for row in rows if isinstance(row, Mapping)})
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        counts[str(row["failure_cohort"])][str(row["protein_regime"])] += 1
    header = "".join(f"<th>{_escape(regime)}</th>" for regime in regimes)
    body = []
    for cohort in cohorts:
        cells = []
        for regime in regimes:
            count = counts[cohort][regime]
            shade = min(count * 18, 90)
            cells.append(
                f"<td style=\"background:rgba(42, 111, 151, {shade / 100});\">"
                f"{_escape(count)}</td>"
            )
        body.append(f"<tr><td>{_escape(cohort)}</td>{''.join(cells)}</tr>")
    return (
        "<section><h2>Failure Cohort Heatmap</h2>"
        "<table><thead><tr><th>cohort</th>"
        + header
        + "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></section>"
    )


def _manual_review_table(report: Mapping[str, object]) -> str:
    rows = report.get("rows_requiring_manual_review", [])
    if not isinstance(rows, Sequence) or not rows:
        return (
            "<section><h2>Manual Review Rows</h2>"
            "<p>No manual review rows recorded.</p></section>"
        )
    body = []
    for row in rows:
        if not isinstance(row, Mapping):
            continue
        body.append(
            "<tr>"
            f"<td>{_escape(row.get('protein_id', ''))}</td>"
            f"<td>{_escape(row.get('protein_regime', ''))}</td>"
            f"<td>{_escape(row.get('predicted_fold_class', ''))}</td>"
            f"<td>{_escape(row.get('structure_fold_class', ''))}</td>"
            f"<td>{_escape(row.get('label_fold_class', ''))}</td>"
            f"<td>{_escape(row.get('failure_cohort', ''))}</td>"
            f"<td>{_escape(row.get('manual_review_reasons', ''))}</td>"
            "</tr>"
        )
    return (
        "<section><h2>Manual Review Rows</h2>"
        "<table><thead><tr><th>protein</th><th>regime</th><th>predicted</th>"
        "<th>structure</th><th>label</th><th>cohort</th><th>review reason</th>"
        "</tr></thead><tbody>"
        + "".join(body)
        + "</tbody></table></section>"
    )


def render_regime_dashboard(report: Mapping[str, object]) -> str:
    title = "Protein Regime Routing Dashboard"
    return (
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
        f"<title>{_escape(title)}</title>"
        "<style>"
        "body{font-family:Arial,sans-serif;margin:0;background:#f8faf7;color:#1f2933;}"
        "header{background:#20313d;color:white;padding:28px 32px;}"
        "main{padding:24px 32px;}"
        ".metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:18px 0;}"
        ".metric{background:white;border:1px solid #d7ddd8;border-radius:6px;padding:14px;}"
        ".metric span{display:block;font-size:12px;color:#52616b;margin-bottom:6px;}"
        ".metric strong{font-size:20px;}"
        "section{margin:24px 0;}"
        "table{border-collapse:collapse;width:100%;background:white;border:1px solid #d7ddd8;}"
        "th,td{border:1px solid #d7ddd8;padding:8px;text-align:left;vertical-align:top;font-size:13px;}"
        "th{background:#e9efec;}"
        "code{background:#eef3f0;padding:2px 4px;border-radius:4px;}"
        "</style></head><body>"
        f"<header><h1>{_escape(title)}</h1>"
        "<p>Sequence-only protein regimes route broad fold gates before "
        "structure/label scoring.</p></header><main>"
        f"<div class=\"metrics\">{_metric_cards(report)}</div>"
        + _mapping_table("Accuracy By Regime", report.get("accuracy_by_regime"))
        + _mapping_table("Abstention By Regime", report.get("abstention_by_regime"))
        + _mapping_table(
            "High-Confidence Wrong By Regime",
            report.get("high_confidence_wrong_by_regime"),
        )
        + _mapping_table("Gate Path By Regime", report.get("gate_path_by_regime"))
        + _heatmap_table(report)
        + _manual_review_table(report)
        + _mapping_table(
            "Old Failure Counters Stayed Closed",
            report.get("old_failure_modes_still_closed"),
        )
        + _mapping_table(
            "Dominant Failure Cohort",
            report.get("dominant_failure_cohort"),
        )
        + "</main></body></html>"
    )
