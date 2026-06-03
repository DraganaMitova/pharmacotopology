from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

from pharmacotopology.folding_axis_adjudication import UNKNOWN_BY_AXIS
from pharmacotopology.folding_axis_profile import (
    AXIS_PROFILE_SIGNATURE_KIND,
    axis_profile_rows,
)
from pharmacotopology.folding_hierarchical_gates import predict_hierarchical_gate
from pharmacotopology.folding_order_aware import extract_order_aware_features
from pharmacotopology.folding_regime_analysis import (
    REGIME_ANALYSIS_BENCHMARK_KIND,
    detect_protein_regime,
)
from pharmacotopology.folding_structure_benchmark import StructureEvidenceRow
from pharmacotopology.folding_topology import (
    FoldingReferenceExample,
    HYDROPHOBIC_AMINO_ACIDS,
    normalize_sequence,
    sequence_features,
)


ARCHITECTURE_AXIS_BENCHMARK_KIND = "architecture_axis_evidence_adjudication"
ARCHITECTURE_AXIS_SIGNATURE_KIND = "sequence_only_architecture_axis_evidence_packet"
ARCHITECTURE_AXIS_VALUES = (
    "compact_single_domain",
    "multidomain_or_segmented",
    "repeat_like",
    "fragment_scope",
    "unknown",
)
ARCHITECTURE_UNKNOWN = UNKNOWN_BY_AXIS["architecture_axis"]


@dataclass(frozen=True)
class ArchitectureEvidencePacket:
    row_id: str
    sequence_hash: str
    sequence_length: int
    length_band: str
    segmentation_pressure: float
    repeat_pressure: float
    compact_domain_pressure: float
    fragment_scope_pressure: float
    multidomain_pressure: float
    membrane_segmentation_pressure: float
    disorder_segmentation_pressure: float
    evidence_strength: float
    architecture_axis_prediction: str
    architecture_axis_confidence: float
    architecture_axis_claim_allowed: bool
    architecture_axis_abstention_reason: str
    architecture_axis_decision_reason: str
    protein_regime: str
    regime_confidence: float
    composition_regime_contrast: float
    linker_separator_pressure: float
    compact_repeat_support: float
    repeat_window_count: int
    repeat_recurrence_support: float
    repeat_unit_consistency: float
    repeat_vs_compact_margin: float
    hydrophobic_periodicity_only_risk: float
    repeat_compact_single_domain_ambiguity_guard: bool = False
    architecture_secondary_leakage_used: bool = False
    architecture_label_leakage_used: bool = False

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


def _rounded(value: float) -> float:
    return round(_clamp(value), 6)


def _bool_mean(values: Sequence[bool]) -> float:
    if not values:
        return 0.0
    return round(sum(1 for value in values if value) / len(values), 6)


def _known(value: object) -> bool:
    return str(value) != ARCHITECTURE_UNKNOWN


def _length_band(length: int) -> str:
    if length < 70:
        return "short_fragment_under_70"
    if length < 140:
        return "small_70_139"
    if length < 260:
        return "medium_140_259"
    if length < 420:
        return "large_260_419"
    return "very_large_420_plus"


def _sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:16]


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


def _window_profiles(sequence: str, window_size: int = 36) -> tuple[dict[str, float], ...]:
    profiles = []
    for start in range(0, len(sequence), window_size):
        window = sequence[start : start + window_size]
        if len(window) < 18:
            continue
        features = sequence_features(window)
        profiles.append(
            {
                "hydrophobic": float(features["hydrophobic_fraction"]),
                "charged": float(features["charged_fraction"]),
                "aromatic": float(features["aromatic_fraction"]),
                "disorder": float(features["disorder_promoting_fraction"]),
                "breaker": _fraction(window, frozenset("PG")),
            }
        )
    return tuple(profiles)


def _repeat_window_count(sequence: str) -> int:
    profiles = _window_profiles(sequence, window_size=24)
    return sum(
        1
        for profile in profiles
        if profile["hydrophobic"] >= 0.34
        and profile["breaker"] <= 0.18
        and profile["disorder"] <= 0.46
    )


def _repeat_unit_consistency(sequence: str) -> float:
    profiles = _window_profiles(sequence, window_size=24)
    if len(profiles) < 3:
        return 0.0
    hydrophobic_values = [profile["hydrophobic"] for profile in profiles]
    mean_hydrophobic = sum(hydrophobic_values) / len(hydrophobic_values)
    mean_abs_deviation = sum(
        abs(value - mean_hydrophobic) for value in hydrophobic_values
    ) / len(hydrophobic_values)
    return _rounded(1.0 - mean_abs_deviation / 0.28)


def _repeat_recurrence_support(
    sequence: str,
    *,
    hydrophobic_periodicity: float,
    composition_contrast: float,
    linker_pressure: float,
) -> float:
    repeat_windows = _repeat_window_count(sequence)
    unit_consistency = _repeat_unit_consistency(sequence)
    multi_window_support = _clamp((repeat_windows - 2) / 4)
    return _rounded(
        multi_window_support * 0.36
        + unit_consistency * 0.28
        + hydrophobic_periodicity * 0.22
        + composition_contrast * 0.10
        - linker_pressure * 0.18
    )


def _hydrophobic_periodicity_only_risk(
    *,
    hydrophobic_periodicity: float,
    repeat_recurrence_support: float,
    composition_contrast: float,
    linker_pressure: float,
    long_range_closure_evidence: float,
) -> float:
    return _rounded(
        hydrophobic_periodicity * 0.42
        + max(0.12 - composition_contrast, 0.0) * 2.30
        + max(0.46 - repeat_recurrence_support, 0.0) * 0.46
        + long_range_closure_evidence * 0.10
        - linker_pressure * 0.22
    )


def _composition_regime_contrast(sequence: str) -> float:
    profiles = _window_profiles(sequence)
    if len(profiles) < 2:
        return 0.0
    contrasts = []
    for left, right in zip(profiles, profiles[1:]):
        contrasts.append(
            (
                abs(left["hydrophobic"] - right["hydrophobic"])
                + abs(left["charged"] - right["charged"])
                + abs(left["aromatic"] - right["aromatic"])
                + abs(left["disorder"] - right["disorder"])
                + abs(left["breaker"] - right["breaker"])
            )
            / 5
        )
    return round(sum(contrasts) / len(contrasts), 6)


def _linker_separator_pressure(sequence: str) -> float:
    profiles = _window_profiles(sequence, window_size=24)
    if not profiles:
        return 0.0
    linker_like = [
        profile
        for profile in profiles
        if profile["breaker"] >= 0.13
        or (profile["charged"] >= 0.32 and profile["hydrophobic"] <= 0.35)
        or profile["disorder"] >= 0.50
    ]
    return round(len(linker_like) / len(profiles), 6)


def architecture_evidence_packet_from_sequence(
    sequence: str,
    *,
    protein_id: str = "sequence",
    external_safe_quarantine: bool = False,
) -> ArchitectureEvidencePacket:
    normalized = normalize_sequence(sequence)
    length = len(normalized)
    features = sequence_features(normalized)
    order_features, _edges = extract_order_aware_features(
        normalized,
        protein_id=protein_id,
    )
    regime = detect_protein_regime(normalized, protein_id=protein_id)
    hierarchy = predict_hierarchical_gate(normalized, protein_id=protein_id)
    evidence = hierarchy.motif_evidence
    scores = hierarchy.gate_scores

    hydrophobic_fraction = float(features["hydrophobic_fraction"])
    hydrophobic_run = _max_run_fraction(normalized, HYDROPHOBIC_AMINO_ACIDS)
    composition_contrast = _composition_regime_contrast(normalized)
    linker_pressure = _linker_separator_pressure(normalized)
    hydrophobic_periodicity = float(order_features["hydrophobic_cluster_periodicity"])
    repeat_window_count = _repeat_window_count(normalized)
    repeat_unit_consistency = _repeat_unit_consistency(normalized)
    repeat_recurrence_support = _repeat_recurrence_support(
        normalized,
        hydrophobic_periodicity=hydrophobic_periodicity,
        composition_contrast=composition_contrast,
        linker_pressure=linker_pressure,
    )
    compact_repeat_support = _rounded(
        regime.repeat_pressure * 0.44
        + hydrophobic_periodicity * 0.30
        + float(order_features["contact_prior_density"]) * 0.18
        - regime.intrinsically_disordered_pressure * 0.20
    )
    segmentation_pressure = _rounded(
        evidence.domain_boundary_evidence * 0.34
        + scores.segmentation_gate_score * 0.30
        + float(order_features["segment_boundary_contrast"]) * 0.75
        + composition_contrast * 0.38
        + linker_pressure * 0.18
        - regime.intrinsically_disordered_pressure * 0.12
    )
    repeat_pressure = _rounded(
        regime.repeat_pressure * 0.58
        + hydrophobic_periodicity * 0.24
        + compact_repeat_support * 0.20
        - linker_pressure * 0.12
    )
    compact_domain_pressure = _rounded(
        regime.compact_single_domain_pressure * 0.48
        + scores.compactness_gate_score * 0.32
        + evidence.long_range_closure_evidence * 0.18
        + float(order_features["contact_prior_mean_weight"]) * 0.34
        - segmentation_pressure * 0.16
        - regime.intrinsically_disordered_pressure * 0.10
    )
    fragment_scope_pressure = _rounded(
        (1.0 if length < 70 else 0.0) * 0.78
        + _clamp((95 - length) / 40) * 0.22
        + (0.12 if scores.compactness_gate_score >= 0.42 else 0.0)
        - segmentation_pressure * 0.08
    )
    multidomain_pressure = _rounded(
        regime.multidomain_modular_pressure * 0.42
        + evidence.domain_boundary_evidence * 0.27
        + segmentation_pressure * 0.24
        + linker_pressure * 0.16
        + _clamp((length - 260) / 360) * 0.18
        - repeat_pressure * 0.16
        - regime.intrinsically_disordered_pressure * 0.10
    )
    membrane_segmentation_pressure = _rounded(
        (0.42 if regime.protein_regime == "membrane_like" else 0.0)
        + max(hydrophobic_fraction - 0.44, 0.0) * 1.55
        + hydrophobic_run * 1.75
        + float(order_features["hydrophobic_cluster_max_fraction"]) * 1.20
    )
    disorder_segmentation_pressure = _rounded(
        regime.intrinsically_disordered_pressure * 0.54
        + evidence.disorder_run_evidence * 0.28
        + float(order_features["breaker_density"]) * 0.58
        + linker_pressure * 0.16
        - compact_domain_pressure * 0.10
    )

    prediction = ARCHITECTURE_UNKNOWN
    confidence = 0.0
    claim_allowed = False
    abstention_reason = "architecture evidence is insufficient or mixed"
    decision_reason = abstention_reason
    evidence_strength = max(
        compact_domain_pressure,
        fragment_scope_pressure,
        multidomain_pressure,
        repeat_pressure,
    )
    repeat_vs_compact_margin = round(repeat_pressure - compact_domain_pressure, 6)
    hydrophobic_periodicity_only_risk = _hydrophobic_periodicity_only_risk(
        hydrophobic_periodicity=hydrophobic_periodicity,
        repeat_recurrence_support=repeat_recurrence_support,
        composition_contrast=composition_contrast,
        linker_pressure=linker_pressure,
        long_range_closure_evidence=evidence.long_range_closure_evidence,
    )
    repeat_compact_guard = (
        external_safe_quarantine
        and 75 <= length <= 190
        and repeat_pressure >= 0.60
        and repeat_vs_compact_margin <= 0.34
        and evidence.long_range_closure_evidence >= 0.34
        and composition_contrast <= 0.075
        and hydrophobic_periodicity_only_risk >= 0.40
        and repeat_recurrence_support < 0.70
    )

    if length < 70 and fragment_scope_pressure >= 0.70:
        prediction = "fragment_scope"
        confidence = _rounded(fragment_scope_pressure)
        claim_allowed = True
        abstention_reason = ""
        decision_reason = "short sequence scope is treated as fragment context"
    elif (
        75 <= length <= 170
        and regime.protein_regime == "compact_single_domain"
        and compact_domain_pressure >= 0.43
        and multidomain_pressure <= 0.36
        and segmentation_pressure <= 0.43
        and repeat_pressure <= 0.36
        and disorder_segmentation_pressure <= 0.50
        and membrane_segmentation_pressure <= 0.28
        and evidence.long_range_closure_evidence >= 0.43
    ):
        prediction = "compact_single_domain"
        confidence = _rounded(
            compact_domain_pressure * 0.72
            + (1.0 - segmentation_pressure) * 0.16
            + (1.0 - disorder_segmentation_pressure) * 0.12
        )
        claim_allowed = True
        abstention_reason = ""
        decision_reason = (
            "small-to-medium compact evidence has low segmentation, repeat, "
            "membrane, and disorder pressure"
        )
    elif (
        length >= 400
        and regime.protein_regime == "multidomain_modular"
        and multidomain_pressure >= 0.43
        and evidence.domain_boundary_evidence >= 0.52
        and segmentation_pressure >= 0.36
        and repeat_pressure <= 0.35
        and membrane_segmentation_pressure <= 0.36
        and disorder_segmentation_pressure <= 0.42
    ):
        prediction = "multidomain_or_segmented"
        confidence = _rounded(
            multidomain_pressure * 0.70
            + segmentation_pressure * 0.16
            + (1.0 - repeat_pressure) * 0.08
            + (1.0 - disorder_segmentation_pressure) * 0.06
        )
        claim_allowed = True
        abstention_reason = ""
        decision_reason = (
            "large modular evidence has boundary pressure and is not better "
            "explained as repeat, membrane, or disorder segmentation"
        )
    elif (
        repeat_pressure >= 0.62
        and compact_repeat_support >= 0.48
        and disorder_segmentation_pressure <= 0.40
        and membrane_segmentation_pressure <= 0.36
    ):
        if repeat_compact_guard:
            abstention_reason = "repeat_compact_single_domain_ambiguity_quarantine"
            decision_reason = (
                "repeat-like architecture was not claimed because the signal is "
                "mostly hydrophobic periodicity with compact-domain closure and "
                "weak recurrence support"
            )
        else:
            prediction = "repeat_like"
            confidence = _rounded(
                repeat_pressure * 0.80
                + compact_repeat_support * 0.14
                + (1.0 - disorder_segmentation_pressure) * 0.06
            )
            claim_allowed = True
            abstention_reason = ""
            decision_reason = (
                "repeat evidence is stable and not better explained as disorder"
            )
    elif length < 75:
        abstention_reason = "short sequence did not clear fragment-scope confidence"
        decision_reason = abstention_reason
    elif regime.protein_regime == "membrane_like":
        abstention_reason = (
            "membrane segmentation can reflect helices rather than global domains"
        )
        decision_reason = abstention_reason
    elif repeat_pressure > 0.35:
        abstention_reason = "repeat/segmentation signals are mixed"
        decision_reason = abstention_reason
    elif disorder_segmentation_pressure > 0.50:
        abstention_reason = "disorder-like segmentation pressure is too high"
        decision_reason = abstention_reason
    elif length > 170 and regime.protein_regime == "compact_single_domain":
        abstention_reason = (
            "compact regime is present but length/segmentation evidence needs "
            "architecture-specific support"
        )
        decision_reason = abstention_reason
    elif regime.protein_regime == "multidomain_modular":
        abstention_reason = (
            "modular regime lacks enough non-repeat architecture evidence"
        )
        decision_reason = abstention_reason

    return ArchitectureEvidencePacket(
        row_id=protein_id,
        sequence_hash=_sequence_hash(normalized),
        sequence_length=length,
        length_band=_length_band(length),
        segmentation_pressure=segmentation_pressure,
        repeat_pressure=repeat_pressure,
        compact_domain_pressure=compact_domain_pressure,
        fragment_scope_pressure=fragment_scope_pressure,
        multidomain_pressure=multidomain_pressure,
        membrane_segmentation_pressure=membrane_segmentation_pressure,
        disorder_segmentation_pressure=disorder_segmentation_pressure,
        evidence_strength=_rounded(evidence_strength),
        architecture_axis_prediction=prediction,
        architecture_axis_confidence=confidence,
        architecture_axis_claim_allowed=claim_allowed,
        architecture_axis_abstention_reason=abstention_reason,
        architecture_axis_decision_reason=decision_reason,
        protein_regime=regime.protein_regime,
        regime_confidence=regime.regime_confidence,
        composition_regime_contrast=composition_contrast,
        linker_separator_pressure=linker_pressure,
        compact_repeat_support=compact_repeat_support,
        repeat_window_count=repeat_window_count,
        repeat_recurrence_support=repeat_recurrence_support,
        repeat_unit_consistency=repeat_unit_consistency,
        repeat_vs_compact_margin=repeat_vs_compact_margin,
        hydrophobic_periodicity_only_risk=hydrophobic_periodicity_only_risk,
        repeat_compact_single_domain_ambiguity_guard=repeat_compact_guard,
    )


def architecture_axis_rows(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
) -> list[dict[str, object]]:
    profile_rows_by_id = {
        str(row["protein_id"]): row
        for row in axis_profile_rows(references, evidence_rows)
    }
    rows: list[dict[str, object]] = []
    for reference in references:
        packet = architecture_evidence_packet_from_sequence(
            reference.sequence,
            protein_id=reference.protein_id,
        )
        profile_row = profile_rows_by_id[reference.protein_id]
        truth = str(profile_row["adjudicated_truth_architecture_axis"])
        prediction = packet.architecture_axis_prediction
        scorable = _known(prediction) and _known(truth)
        match = scorable and prediction == truth
        conflict = scorable and not match
        recovered_from_profile = (
            _known(prediction)
            and prediction == str(profile_row["profile_architecture_axis"])
        )
        row = packet.to_dict()
        row.update(
            {
                "architecture_axis_signature_kind": (
                    ARCHITECTURE_AXIS_SIGNATURE_KIND
                ),
                "source_axis_profile_signature_kind": (
                    AXIS_PROFILE_SIGNATURE_KIND
                ),
                "profile_architecture_axis": profile_row[
                    "profile_architecture_axis"
                ],
                "adjudicated_truth_architecture_axis": truth,
                "architecture_axis_scorable": scorable,
                "architecture_axis_match": match,
                "architecture_axis_same_axis_conflict": conflict,
                "architecture_axis_recovered_from_profile": recovered_from_profile,
                "global_fold_class_claim_allowed": False,
                "axis_profile_claim_allowed": True,
                "folding_problem_solved": False,
                "folding_solution_claim_created": False,
                "claim_allowed": False,
                "drug_design_created": False,
                "molecule_generated": False,
                "protein_sequence_design_created": False,
            }
        )
        rows.append(row)
    return rows


def architecture_axis_conflict_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    keys = (
        "row_id",
        "architecture_axis_prediction",
        "adjudicated_truth_architecture_axis",
        "architecture_axis_confidence",
        "architecture_axis_decision_reason",
        "protein_regime",
        "length_band",
        "segmentation_pressure",
        "repeat_pressure",
        "repeat_recurrence_support",
        "repeat_unit_consistency",
        "repeat_vs_compact_margin",
        "hydrophobic_periodicity_only_risk",
        "compact_domain_pressure",
        "multidomain_pressure",
        "membrane_segmentation_pressure",
        "disorder_segmentation_pressure",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if bool(row["architecture_axis_same_axis_conflict"])
    ]


def architecture_axis_abstention_rows(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    keys = (
        "row_id",
        "protein_regime",
        "length_band",
        "architecture_axis_prediction",
        "architecture_axis_abstention_reason",
        "evidence_strength",
        "segmentation_pressure",
        "repeat_pressure",
        "repeat_recurrence_support",
        "repeat_unit_consistency",
        "repeat_vs_compact_margin",
        "hydrophobic_periodicity_only_risk",
        "compact_domain_pressure",
        "fragment_scope_pressure",
        "multidomain_pressure",
        "membrane_segmentation_pressure",
        "disorder_segmentation_pressure",
        "adjudicated_truth_architecture_axis",
    )
    return [
        {key: row[key] for key in keys}
        for row in rows
        if not bool(row["architecture_axis_claim_allowed"])
    ]


def _coverage(rows: Sequence[Mapping[str, object]]) -> float:
    return _bool_mean(
        [bool(row["architecture_axis_claim_allowed"]) for row in rows]
    )


def _claim_count(rows: Sequence[Mapping[str, object]], axis_value: str) -> int:
    return sum(
        1
        for row in rows
        if bool(row["architecture_axis_claim_allowed"])
        and row["architecture_axis_prediction"] == axis_value
    )


def _claim_distribution(rows: Sequence[Mapping[str, object]]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for row in rows:
        if bool(row["architecture_axis_claim_allowed"]):
            counts[str(row["architecture_axis_prediction"])] += 1
    return {value: counts[value] for value in ARCHITECTURE_AXIS_VALUES}


def build_architecture_axis_report(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
    *,
    source_benchmark_file: Path,
    structure_evidence_file: Path,
) -> dict[str, object]:
    rows = architecture_axis_rows(references, evidence_rows)
    claim_allowed_rows = [
        row for row in rows if bool(row["architecture_axis_claim_allowed"])
    ]
    conflicts = [
        row for row in rows if bool(row["architecture_axis_same_axis_conflict"])
    ]
    safe_claims = [
        row
        for row in claim_allowed_rows
        if not bool(row["architecture_axis_same_axis_conflict"])
    ]
    multidomain_truth_rows = [
        row
        for row in rows
        if row["adjudicated_truth_architecture_axis"] == "multidomain_or_segmented"
    ]
    architecture_axis_coverage = _coverage(rows)
    previous_profile_architecture_axis_coverage = _bool_mean(
        [_known(row["profile_architecture_axis"]) for row in rows]
    )
    return {
        "benchmark_kind": ARCHITECTURE_AXIS_BENCHMARK_KIND,
        "source_regime_analysis_benchmark_kind": REGIME_ANALYSIS_BENCHMARK_KIND,
        "architecture_axis_signature_kind": ARCHITECTURE_AXIS_SIGNATURE_KIND,
        "source_axis_profile_signature_kind": AXIS_PROFILE_SIGNATURE_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "structure_evidence_file": str(structure_evidence_file),
        "predictor_input_boundary": "sequence_only_no_labels_no_structure_answers",
        "truth_scoring_boundary": (
            "labels_structure_sources_and_reference_axes_used_only_after_architecture_prediction"
        ),
        "benchmark_size": len(rows),
        "architecture_axis_coverage": architecture_axis_coverage,
        "previous_profile_architecture_axis_coverage": (
            previous_profile_architecture_axis_coverage
        ),
        "architecture_axis_claim_allowed_count": len(claim_allowed_rows),
        "architecture_axis_abstained_count": len(rows) - len(claim_allowed_rows),
        "architecture_axis_conflict_count": len(conflicts),
        "architecture_axis_safe_claim_count": len(safe_claims),
        "architecture_axis_same_axis_conflict_count": len(conflicts),
        "architecture_axis_accuracy": _bool_mean(
            [
                bool(row["architecture_axis_match"])
                for row in claim_allowed_rows
                if bool(row["architecture_axis_scorable"])
            ]
        ),
        "fragment_scope_detected_count": _claim_count(rows, "fragment_scope"),
        "multidomain_claim_count": _claim_count(rows, "multidomain_or_segmented"),
        "multidomain_abstained_count": sum(
            1
            for row in multidomain_truth_rows
            if row["architecture_axis_prediction"] != "multidomain_or_segmented"
        ),
        "repeat_like_claim_count": _claim_count(rows, "repeat_like"),
        "compact_single_domain_claim_count": _claim_count(
            rows,
            "compact_single_domain",
        ),
        "architecture_claim_distribution": _claim_distribution(rows),
        "architecture_claim_without_secondary_leakage_count": sum(
            1
            for row in claim_allowed_rows
            if bool(row["architecture_secondary_leakage_used"])
        ),
        "architecture_claim_without_label_leakage_count": sum(
            1
            for row in claim_allowed_rows
            if bool(row["architecture_label_leakage_used"])
        ),
        "architecture_axis_recovered_from_profile_count": sum(
            1
            for row in claim_allowed_rows
            if bool(row["architecture_axis_recovered_from_profile"])
        ),
        "global_fold_class_claim_allowed": False,
        "axis_profile_claim_allowed": True,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "claim_allowed": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "artifact_reproducible": True,
        "boundary_statement": (
            "This layer creates sequence-only architecture evidence packets and "
            "scores them against architecture truth only after prediction. It "
            "does not use secondary labels as architecture evidence, recover "
            "global fold-class coverage, override safety guards, export raw "
            "sequences, or claim that folding is solved."
        ),
        "rows": rows,
    }


def build_architecture_axis_certificate(
    report: Mapping[str, object],
    *,
    output_names: Sequence[str],
) -> dict[str, object]:
    return {
        "certificate_kind": "architecture_axis_evidence_safety_certificate",
        "benchmark_kind": report["benchmark_kind"],
        "architecture_axis_signature_kind": report[
            "architecture_axis_signature_kind"
        ],
        "source_benchmark_file": report["source_benchmark_file"],
        "structure_evidence_file": report["structure_evidence_file"],
        "architecture_axis_coverage": report["architecture_axis_coverage"],
        "architecture_axis_claim_allowed_count": report[
            "architecture_axis_claim_allowed_count"
        ],
        "architecture_axis_same_axis_conflict_count": report[
            "architecture_axis_same_axis_conflict_count"
        ],
        "architecture_claim_without_secondary_leakage_count": report[
            "architecture_claim_without_secondary_leakage_count"
        ],
        "architecture_claim_without_label_leakage_count": report[
            "architecture_claim_without_label_leakage_count"
        ],
        "global_fold_class_claim_allowed": report[
            "global_fold_class_claim_allowed"
        ],
        "axis_profile_claim_allowed": report["axis_profile_claim_allowed"],
        "folding_problem_solved": report["folding_problem_solved"],
        "raw_sequences_exported": False,
        "output_artifacts": tuple(output_names),
    }


def write_architecture_axis_outputs(
    *,
    report: Mapping[str, object],
    rows: Sequence[Mapping[str, object]],
    conflicts: Sequence[Mapping[str, object]],
    abstentions: Sequence[Mapping[str, object]],
    report_path: Path,
    rows_path: Path,
    conflicts_path: Path,
    abstentions_path: Path,
    dashboard_path: Path,
    certificate_path: Path,
) -> tuple[Path, Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    output_names = (
        report_path.name,
        rows_path.name,
        conflicts_path.name,
        abstentions_path.name,
        dashboard_path.name,
        certificate_path.name,
    )
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv_rows(rows, rows_path)
    _write_csv_rows(conflicts, conflicts_path)
    _write_csv_rows(abstentions, abstentions_path)
    dashboard_path.write_text(render_architecture_axis_dashboard(report), encoding="utf-8")
    certificate = build_architecture_axis_certificate(
        report,
        output_names=output_names,
    )
    certificate_path.write_text(
        json.dumps(certificate, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return (
        report_path,
        rows_path,
        conflicts_path,
        abstentions_path,
        dashboard_path,
        certificate_path,
    )


def _write_csv_rows(
    rows: Sequence[Mapping[str, object]],
    path: Path,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as file:
        if not rows:
            return path
        fieldnames = list(rows[0])
        writer = csv.DictWriter(file, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
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
        "architecture_axis_coverage",
        "previous_profile_architecture_axis_coverage",
        "architecture_axis_claim_allowed_count",
        "architecture_axis_same_axis_conflict_count",
        "fragment_scope_detected_count",
        "compact_single_domain_claim_count",
        "multidomain_claim_count",
        "repeat_like_claim_count",
        "claim_allowed",
        "folding_problem_solved",
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


def _claim_preview(report: Mapping[str, object]) -> str:
    rows = report.get("rows", [])
    if not isinstance(rows, Sequence):
        return ""
    preview_rows = [
        row
        for row in rows
        if isinstance(row, Mapping) and bool(row["architecture_axis_claim_allowed"])
    ][:32]
    body = "".join(
        "<tr>"
        f"<td>{_escape(row['row_id'])}</td>"
        f"<td>{_escape(row['length_band'])}</td>"
        f"<td>{_escape(row['protein_regime'])}</td>"
        f"<td>{_escape(row['architecture_axis_prediction'])}</td>"
        f"<td>{_escape(row['architecture_axis_confidence'])}</td>"
        f"<td>{_escape(row['architecture_axis_decision_reason'])}</td>"
        "</tr>"
        for row in preview_rows
    )
    return (
        "<section><h2>Architecture Axis Claims</h2>"
        "<table><thead><tr>"
        "<th>row_id</th><th>length_band</th><th>regime</th><th>axis</th>"
        "<th>confidence</th><th>decision reason</th>"
        "</tr></thead><tbody>"
        + body
        + "</tbody></table></section>"
    )


def render_architecture_axis_dashboard(report: Mapping[str, object]) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Architecture Axis Evidence Adjudication</title>
  <style>
    body {{
      margin: 0;
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f6f7f2;
      color: #1f2523;
    }}
    header {{
      padding: 32px;
      background: #22302d;
      color: #f6f7f2;
    }}
    main {{
      max-width: 1180px;
      margin: 0 auto;
      padding: 24px;
    }}
    h1, h2 {{
      margin: 0 0 12px;
      letter-spacing: 0;
    }}
    section {{
      margin: 24px 0;
    }}
    .metrics {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
      margin-top: 20px;
    }}
    .metric {{
      background: #ffffff;
      border: 1px solid #d5ddd5;
      border-radius: 6px;
      padding: 14px;
    }}
    .metric span {{
      display: block;
      color: #58635e;
      font-size: 12px;
      overflow-wrap: anywhere;
    }}
    .metric strong {{
      display: block;
      margin-top: 8px;
      font-size: 24px;
    }}
    .rule-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
      gap: 12px;
    }}
    .rule {{
      background: #fff;
      border: 1px solid #d5ddd5;
      border-radius: 6px;
      padding: 14px;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: #fff;
      border: 1px solid #d5ddd5;
      border-radius: 6px;
      overflow: hidden;
    }}
    th, td {{
      padding: 9px 10px;
      border-bottom: 1px solid #e8eee8;
      text-align: left;
      font-size: 13px;
      vertical-align: top;
      overflow-wrap: anywhere;
    }}
    th {{
      background: #e8eee8;
      color: #34423d;
    }}
  </style>
</head>
<body>
  <header>
    <h1>Architecture Axis Evidence Adjudication</h1>
    <p>Sequence-side architecture claims are scored after prediction, separate from secondary structure and global fold class.</p>
    <div class="metrics">{_metric_cards(report)}</div>
  </header>
  <main>
    <section>
      <h2>Boundary Rules</h2>
      <div class="rule-grid">
        <div class="rule"><strong>Architecture Is Not Secondary Structure</strong><br>Alpha, beta, and mixed labels are not used as architecture evidence.</div>
        <div class="rule"><strong>Length Is Not Enough</strong><br>Long sequences need boundary evidence and must not be repeat, membrane, or disorder explanations.</div>
        <div class="rule"><strong>Fragments Stay Scoped</strong><br>Short solved chains are treated as fragment scope, not global fold truth.</div>
        <div class="rule"><strong>No Global Class Recovery</strong><br>The collapsed fold class remains refused.</div>
      </div>
    </section>
    {_mapping_table("Claim Distribution", report.get("architecture_claim_distribution"))}
    {_claim_preview(report)}
  </main>
</body>
</html>
"""
