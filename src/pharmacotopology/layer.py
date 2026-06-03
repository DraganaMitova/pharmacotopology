from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from pharmacotopology.field_alphabet import (
    FIELD_CLOSED_CHANNELS,
    FieldAction,
    FieldEvent,
    FieldKey,
    FieldResponse,
    FieldStatus,
)
from pharmacotopology.input_quarantine import InputQuarantine
from pharmacotopology.session_record import (
    ManualSessionRecord,
    boundary_status_dict,
    now_timestamp,
)
from pharmacotopology.state_certificate import create_evidence_sealed_certificate
from pharmacotopology.trace_io import append_jsonl, event_group


PATHOLOGY_DIMENSIONS: tuple[str, ...] = (
    "salience_amplification",
    "recurrence_overbinding",
    "symbolic_closure_pressure",
    "threat_propagation",
    "falsification_weakness",
    "boundary_instability",
    "agency_confusion",
    "sensory_intrusion",
    "cognitive_fragmentation",
    "negative_shutdown",
    "sleep_instability",
)

TOPOLOGY_PRESSURE_MIN = 0.0
TOPOLOGY_PRESSURE_MAX = 1.25

DEFAULT_PHARMACOTOPOLOGY_SURFACE = "ι.clean.pharmacotopology.simulation.request"

PHARMACOTOPOLOGY_GENERATED_FILES: tuple[str, ...] = (
    "input_stream.jsonl",
    "memory.jsonl",
    "audit.jsonl",
    "session_records.jsonl",
    "clean_pharmacotopology_layer_report.json",
    "calibration_readiness_report.json",
    "pharmacotopology_rankings.csv",
    "pharmacotopology_deltas.csv",
    "sensitivity_analysis_report.json",
    "sensitivity_rankings.csv",
    "pharmacotopology_dashboard.html",
    "multi_profile_dashboard.html",
    "sensitivity_explorer_report.json",
    "sensitivity_explorer_samples.csv",
    "folding_topology_benchmark_report.json",
    "folding_topology_benchmark.csv",
    "real_folding_500_report.json",
    "real_folding_500_rows.csv",
    "real_folding_500_dashboard.html",
    "real_folding_500_certificate.json",
    "real_folding_500_failures.csv",
    "real_folding_500_confusion_matrix.csv",
    "real_folding_10_report.json",
    "real_folding_10_rows.csv",
    "real_folding_10_dashboard.html",
    "real_folding_10_certificate.json",
    "real_folding_10_failures.csv",
    "real_folding_10_confusion_matrix.csv",
    "real_folding_10_structure_report.json",
    "real_folding_10_structure_rows.csv",
    "real_folding_10_structure_dashboard.html",
    "real_folding_10_order_controls.csv",
    "real_folding_10_falsification_report.json",
    "real_folding_10_order_aware_report.json",
    "real_folding_10_order_aware_rows.csv",
    "real_folding_10_contact_prior.csv",
    "real_folding_10_control_separation.csv",
    "real_folding_10_order_aware_dashboard.html",
    "real_folding_10_motif_alignment_report.json",
    "real_folding_10_motif_alignment_rows.csv",
    "real_folding_10_failure_diagnosis.csv",
    "real_folding_10_evidence_conflicts.csv",
    "real_folding_10_motif_alignment_dashboard.html",
    "field_validation.json",
    "field_metrics.json",
)

EVIDENCE_STAGE_WEIGHTS: dict[str, float] = {
    "hypothesis_only": 0.0,
    "mechanism_proxy": 0.25,
    "external_dataset_proxy": 0.50,
    "externally_calibrated": 0.75,
    "replicated_external": 0.90,
}


@dataclass(frozen=True)
class TopologyProfile:
    profile_id: str
    description: str
    dimensions: dict[str, float]


@dataclass(frozen=True)
class MechanismVector:
    mechanism_id: str
    mechanism_family: str
    receptor_profile: tuple[str, ...]
    abstract_compound_class: str
    protein_family: str
    protein_mechanism_class: str
    protein_state_shift: str
    pathway_network_perturbation: str
    deltas: dict[str, float]
    collapse_cost: float
    evidence_stage: str = "hypothesis_only"
    evidence_weight: float = 0.0
    uncertainty_radius: float = 0.35
    primary_evidence_sources: tuple[str, ...] = ()
    evidence_refs: tuple[str, ...] = ()
    calibration_blockers: tuple[str, ...] = (
        "no_primary_evidence_source_attached",
        "no_empirical_delta_interval",
    )
    confidence_interval_kind: str = "model_uncertainty_interval"
    assumption_notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class PerturbationResult:
    mechanism_id: str
    mechanism_family: str
    receptor_profile: tuple[str, ...]
    abstract_compound_class: str
    protein_family: str
    protein_mechanism_class: str
    protein_state_shift: str
    pathway_network_perturbation: str
    evidence_stage: str
    evidence_weight: float
    uncertainty_radius: float
    primary_evidence_sources: tuple[str, ...]
    evidence_refs: tuple[str, ...]
    calibration_blockers: tuple[str, ...]
    confidence_interval_kind: str
    assumption_notes: tuple[str, ...]
    resulting_state: dict[str, float]
    topology_delta: dict[str, float]
    scored_dimensions: tuple[str, ...]
    improved_dimensions: tuple[str, ...]
    worsened_dimensions: tuple[str, ...]
    pathology_reduction_score: float
    pathology_reduction_interval: dict[str, float]
    collapse_cost_score: float
    collapse_cost_interval: dict[str, float]
    net_topology_health_score: float
    net_topology_health_interval: dict[str, float]
    fit_label: str
    evidence_readiness_label: str


@dataclass(frozen=True)
class CalibrationReadinessReport:
    calibration_status: str
    practical_use: str
    mechanism_vectors_reviewed: int
    evidence_backed_vectors: int
    uncalibrated_vectors: int
    mean_evidence_weight: float
    mean_uncertainty_radius: float
    clinical_use_allowed: bool
    external_validation_required: bool
    blockers: tuple[str, ...]
    next_steps: tuple[str, ...]


@dataclass(frozen=True)
class CleanPharmacotopologyLayerReport:
    run_kind: str
    pharmacotopology_review_valid: bool
    simulation_only: bool
    hypothesis_numbers_only: bool
    calibration_status: str
    practical_use: str
    evidence_backed_vectors: int
    uncalibrated_vectors: int
    mean_evidence_weight: float
    mean_uncertainty_radius: float
    clinical_use_allowed: bool
    mechanism_vectors_reviewed: int
    topology_dimensions_reviewed: int
    top_mechanism_id: str
    top_net_topology_health_score: float
    top_pathology_reduction_score: float
    top_collapse_cost_score: float
    destabilizing_mechanism_count: int
    clinical_advice_created: bool
    medication_recommendation_created: bool
    real_patient_inference_created: bool
    brand_name_mapping_created: bool
    treatment_claim_created: bool
    voice_opened: bool
    schedule_created: bool
    command_queue_created: bool
    autonomous_loop_created: bool
    surface_leakage: float
    stop_integrity: float


DEFAULT_SCHIZOPHRENIA_LIKE_PROFILE = TopologyProfile(
    profile_id="schizophrenia_like_topology_profile",
    description=(
        "Simulated unstable topology profile only; not a diagnosis, patient model, "
        "or medication target."
    ),
    dimensions={
        "salience_amplification": 0.91,
        "recurrence_overbinding": 0.86,
        "symbolic_closure_pressure": 0.80,
        "threat_propagation": 0.76,
        "falsification_weakness": 0.72,
        "boundary_instability": 0.68,
        "agency_confusion": 0.70,
        "sensory_intrusion": 0.74,
        "cognitive_fragmentation": 0.66,
        "negative_shutdown": 0.58,
        "sleep_instability": 0.64,
    },
)

DEFAULT_DEPRESSION_LIKE_PROFILE = TopologyProfile(
    profile_id="depression_like_topology_profile",
    description=(
        "Synthetic low-agency shutdown pressure profile only; not a diagnosis, "
        "patient model, or medication target."
    ),
    dimensions={
        "salience_amplification": 0.42,
        "recurrence_overbinding": 0.64,
        "symbolic_closure_pressure": 0.54,
        "threat_propagation": 0.50,
        "falsification_weakness": 0.48,
        "boundary_instability": 0.38,
        "agency_confusion": 0.46,
        "sensory_intrusion": 0.34,
        "cognitive_fragmentation": 0.58,
        "negative_shutdown": 0.88,
        "sleep_instability": 0.72,
    },
)

DEFAULT_MANIA_LIKE_PROFILE = TopologyProfile(
    profile_id="mania_like_topology_profile",
    description=(
        "Synthetic high-salience activation pressure profile only; not a diagnosis, "
        "patient model, or medication target."
    ),
    dimensions={
        "salience_amplification": 0.92,
        "recurrence_overbinding": 0.74,
        "symbolic_closure_pressure": 0.66,
        "threat_propagation": 0.52,
        "falsification_weakness": 0.62,
        "boundary_instability": 0.70,
        "agency_confusion": 0.68,
        "sensory_intrusion": 0.56,
        "cognitive_fragmentation": 0.60,
        "negative_shutdown": 0.28,
        "sleep_instability": 0.90,
    },
)

DEFAULT_ANXIETY_LIKE_PROFILE = TopologyProfile(
    profile_id="anxiety_like_topology_profile",
    description=(
        "Synthetic threat-propagation pressure profile only; not a diagnosis, "
        "patient model, or medication target."
    ),
    dimensions={
        "salience_amplification": 0.62,
        "recurrence_overbinding": 0.66,
        "symbolic_closure_pressure": 0.58,
        "threat_propagation": 0.90,
        "falsification_weakness": 0.54,
        "boundary_instability": 0.60,
        "agency_confusion": 0.42,
        "sensory_intrusion": 0.68,
        "cognitive_fragmentation": 0.50,
        "negative_shutdown": 0.46,
        "sleep_instability": 0.78,
    },
)

DEFAULT_MIXED_STATE_LIKE_PROFILE = TopologyProfile(
    profile_id="mixed_state_like_topology_profile",
    description=(
        "Synthetic simultaneous activation-and-shutdown pressure profile only; "
        "not a diagnosis, patient model, or medication target."
    ),
    dimensions={
        "salience_amplification": 0.82,
        "recurrence_overbinding": 0.78,
        "symbolic_closure_pressure": 0.70,
        "threat_propagation": 0.82,
        "falsification_weakness": 0.60,
        "boundary_instability": 0.76,
        "agency_confusion": 0.64,
        "sensory_intrusion": 0.58,
        "cognitive_fragmentation": 0.68,
        "negative_shutdown": 0.72,
        "sleep_instability": 0.88,
    },
)

DEFAULT_NORMAL_BOUNDED_PROFILE = TopologyProfile(
    profile_id="normal_bounded_topology_profile",
    description=(
        "Reference profile for bounded pattern review; lower pressure means the "
        "dimension is closer to stable review."
    ),
    dimensions={
        "salience_amplification": 0.32,
        "recurrence_overbinding": 0.28,
        "symbolic_closure_pressure": 0.25,
        "threat_propagation": 0.24,
        "falsification_weakness": 0.20,
        "boundary_instability": 0.22,
        "agency_confusion": 0.18,
        "sensory_intrusion": 0.20,
        "cognitive_fragmentation": 0.24,
        "negative_shutdown": 0.26,
        "sleep_instability": 0.20,
    },
)

DEFAULT_TOPOLOGY_PROFILES: dict[str, TopologyProfile] = {
    "schizophrenia_like": DEFAULT_SCHIZOPHRENIA_LIKE_PROFILE,
    "depression_like": DEFAULT_DEPRESSION_LIKE_PROFILE,
    "mania_like": DEFAULT_MANIA_LIKE_PROFILE,
    "anxiety_like": DEFAULT_ANXIETY_LIKE_PROFILE,
    "mixed_state_like": DEFAULT_MIXED_STATE_LIKE_PROFILE,
}


def get_topology_profile(profile_key: str) -> TopologyProfile:
    try:
        return DEFAULT_TOPOLOGY_PROFILES[profile_key]
    except KeyError as exc:
        available = ", ".join(sorted(DEFAULT_TOPOLOGY_PROFILES))
        raise ValueError(
            f"Unknown topology profile {profile_key!r}. Available: {available}"
        ) from exc

DEFAULT_MECHANISM_VECTORS: tuple[MechanismVector, ...] = (
    MechanismVector(
        mechanism_id="d2_antagonist_like",
        mechanism_family="dopamine_salience_dampening",
        receptor_profile=("D2 antagonism",),
        abstract_compound_class="abstract_dopamine_d2_dampener",
        protein_family="dopamine_receptor_gpcr_family",
        protein_mechanism_class="receptor_state_dampening",
        protein_state_shift="reduced_d2_like_signaling_tone",
        pathway_network_perturbation="dopamine_salience_network_dampening",
        deltas={
            "salience_amplification": -0.35,
            "recurrence_overbinding": -0.22,
            "threat_propagation": -0.18,
            "falsification_weakness": -0.10,
            "negative_shutdown": 0.18,
        },
        collapse_cost=0.28,
    ),
    MechanismVector(
        mechanism_id="d2_partial_agonist_like",
        mechanism_family="dopamine_salience_modulation",
        receptor_profile=("D2 partial agonism",),
        abstract_compound_class="abstract_dopamine_d2_modulator",
        protein_family="dopamine_receptor_gpcr_family",
        protein_mechanism_class="receptor_state_modulation",
        protein_state_shift="stabilized_d2_like_signaling_tone",
        pathway_network_perturbation="dopamine_salience_network_modulation",
        deltas={
            "salience_amplification": -0.20,
            "recurrence_overbinding": -0.14,
            "symbolic_closure_pressure": -0.08,
            "negative_shutdown": 0.04,
        },
        collapse_cost=0.12,
    ),
    MechanismVector(
        mechanism_id="5ht2a_antagonist_like",
        mechanism_family="serotonin_sensory_closure_modulation",
        receptor_profile=("5-HT2A antagonism",),
        abstract_compound_class="abstract_serotonin_5ht2a_dampener",
        protein_family="serotonin_receptor_gpcr_family",
        protein_mechanism_class="receptor_state_dampening",
        protein_state_shift="reduced_5ht2a_like_signaling_tone",
        pathway_network_perturbation="serotonin_sensory_closure_modulation",
        deltas={
            "sensory_intrusion": -0.25,
            "symbolic_closure_pressure": -0.20,
            "sleep_instability": -0.18,
        },
        collapse_cost=0.12,
    ),
    MechanismVector(
        mechanism_id="nmda_support_like",
        mechanism_family="glutamate_falsification_support",
        receptor_profile=("NMDA support",),
        abstract_compound_class="abstract_nmda_support_modulator",
        protein_family="glutamate_receptor_ion_channel_family",
        protein_mechanism_class="channel_support_modulation",
        protein_state_shift="supported_nmda_like_gating_function",
        pathway_network_perturbation="glutamate_falsification_network_support",
        deltas={
            "cognitive_fragmentation": -0.20,
            "falsification_weakness": -0.18,
            "boundary_instability": -0.15,
        },
        collapse_cost=0.05,
    ),
    MechanismVector(
        mechanism_id="gaba_stabilizing_like",
        mechanism_family="inhibition_threat_dampening",
        receptor_profile=("GABA modulation",),
        abstract_compound_class="abstract_gaba_stabilizing_modulator",
        protein_family="gaba_receptor_inhibitory_channel_family",
        protein_mechanism_class="inhibitory_channel_modulation",
        protein_state_shift="increased_gaba_like_inhibitory_stability",
        pathway_network_perturbation="inhibition_threat_network_dampening",
        deltas={
            "threat_propagation": -0.20,
            "sensory_intrusion": -0.10,
            "cognitive_fragmentation": 0.07,
        },
        collapse_cost=0.18,
    ),
    MechanismVector(
        mechanism_id="muscarinic_modulation_like",
        mechanism_family="cholinergic_boundary_review",
        receptor_profile=("muscarinic modulation",),
        abstract_compound_class="abstract_muscarinic_modulator",
        protein_family="muscarinic_acetylcholine_receptor_gpcr_family",
        protein_mechanism_class="receptor_state_modulation",
        protein_state_shift="modulated_muscarinic_like_boundary_signaling",
        pathway_network_perturbation="cholinergic_boundary_review_modulation",
        deltas={
            "salience_amplification": -0.18,
            "sensory_intrusion": -0.14,
            "boundary_instability": -0.16,
            "cognitive_fragmentation": -0.08,
        },
        collapse_cost=0.10,
    ),
    MechanismVector(
        mechanism_id="adrenergic_dampening_like",
        mechanism_family="arousal_threat_dampening",
        receptor_profile=("adrenergic dampening",),
        abstract_compound_class="abstract_adrenergic_dampener",
        protein_family="adrenergic_receptor_gpcr_family",
        protein_mechanism_class="arousal_receptor_dampening",
        protein_state_shift="reduced_adrenergic_like_arousal_tone",
        pathway_network_perturbation="arousal_threat_network_dampening",
        deltas={
            "threat_propagation": -0.22,
            "sleep_instability": -0.12,
            "salience_amplification": -0.08,
            "negative_shutdown": 0.05,
        },
        collapse_cost=0.10,
    ),
    MechanismVector(
        mechanism_id="histamine_sedation_load_like",
        mechanism_family="sedation_heavy_arousal_reduction",
        receptor_profile=("histamine blockade", "sedation load"),
        abstract_compound_class="abstract_histamine_sedation_load",
        protein_family="histamine_receptor_gpcr_family",
        protein_mechanism_class="wakefulness_receptor_dampening",
        protein_state_shift="reduced_histamine_like_wakefulness_tone",
        pathway_network_perturbation="sedation_heavy_arousal_reduction",
        deltas={
            "sleep_instability": -0.25,
            "threat_propagation": -0.08,
            "salience_amplification": -0.05,
            "cognitive_fragmentation": 0.12,
            "negative_shutdown": 0.20,
        },
        collapse_cost=0.40,
    ),
    MechanismVector(
        mechanism_id="anticholinergic_burden_like",
        mechanism_family="collapse_cost_control",
        receptor_profile=("anticholinergic burden",),
        abstract_compound_class="abstract_anticholinergic_burden",
        protein_family="acetylcholine_receptor_family",
        protein_mechanism_class="cognitive_burden_receptor_dampening",
        protein_state_shift="reduced_cholinergic_like_support_tone",
        pathway_network_perturbation="cognitive_boundary_collapse_cost_increase",
        deltas={
            "cognitive_fragmentation": 0.24,
            "boundary_instability": 0.08,
            "agency_confusion": 0.08,
            "negative_shutdown": 0.18,
        },
        collapse_cost=0.34,
    ),
    MechanismVector(
        mechanism_id="glutamate_amplification_stressor_like",
        mechanism_family="destabilizing_control",
        receptor_profile=("excitatory amplification", "stress load"),
        abstract_compound_class="abstract_excitatory_amplification_stressor",
        protein_family="glutamate_receptor_signaling_family",
        protein_mechanism_class="excitatory_state_amplification",
        protein_state_shift="increased_glutamate_like_excitatory_pressure",
        pathway_network_perturbation="salience_threat_recurrence_amplification",
        deltas={
            "salience_amplification": 0.21,
            "recurrence_overbinding": 0.19,
            "symbolic_closure_pressure": 0.18,
            "threat_propagation": 0.20,
            "sensory_intrusion": 0.14,
        },
        collapse_cost=0.08,
    ),
)


def reset_pharmacotopology_run_dir(run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    for name in PHARMACOTOPOLOGY_GENERATED_FILES:
        path = run_dir / name
        if path.exists():
            path.unlink()


def _clamp(value: float, lower: float, upper: float) -> float:
    return max(lower, min(upper, value))


def _rounded(value: float) -> float:
    return round(value, 6)


def _dimension_distance(
    state: Mapping[str, float],
    target: Mapping[str, float],
    dimensions: Sequence[str] = PATHOLOGY_DIMENSIONS,
) -> float:
    distances = [abs(float(state[dim]) - float(target[dim])) for dim in dimensions]
    return sum(distances) / len(distances)


def _collapse_cost_score(
    vector: MechanismVector,
    topology_delta: Mapping[str, float],
) -> float:
    return _collapse_cost_score_from_base(vector.collapse_cost, topology_delta)


def _collapse_cost_score_from_base(
    collapse_cost: float,
    topology_delta: Mapping[str, float],
) -> float:
    penalty = (
        max(float(topology_delta.get("negative_shutdown", 0.0)), 0.0) * 0.50
        + max(float(topology_delta.get("cognitive_fragmentation", 0.0)), 0.0)
        * 0.35
        + max(float(topology_delta.get("boundary_instability", 0.0)), 0.0)
        * 0.20
    )
    return _rounded(_clamp(collapse_cost + penalty, 0.0, 1.0))


def _interval(values: Sequence[float]) -> dict[str, float]:
    return {
        "lower": _rounded(min(values)),
        "upper": _rounded(max(values)),
    }


def _effective_evidence_weight(vector: MechanismVector) -> float:
    stage_weight = EVIDENCE_STAGE_WEIGHTS.get(vector.evidence_stage, 0.0)
    return _rounded(_clamp(max(vector.evidence_weight, stage_weight), 0.0, 1.0))


def _effective_uncertainty_radius(vector: MechanismVector) -> float:
    evidence_weight = _effective_evidence_weight(vector)
    return _rounded(
        _clamp(vector.uncertainty_radius * (1.0 - evidence_weight), 0.0, 1.0)
    )


def _evidence_readiness_label(evidence_weight: float) -> str:
    if evidence_weight <= 0.0:
        return "uncalibrated_hypothesis"
    if evidence_weight < 0.40:
        return "mechanism_proxy_only"
    if evidence_weight < 0.75:
        return "partially_calibrated"
    return "externally_calibrated_research"


def _score_scaled_vector(
    source: TopologyProfile,
    target: TopologyProfile,
    vector: MechanismVector,
    *,
    delta_scale: float,
    collapse_scale: float,
) -> tuple[
    dict[str, float],
    dict[str, float],
    tuple[str, ...],
    tuple[str, ...],
    tuple[str, ...],
    float,
    float,
    float,
]:
    scored_dimensions = tuple(
        dimension
        for dimension in PATHOLOGY_DIMENSIONS
        if abs(float(vector.deltas.get(dimension, 0.0))) > 0.0
    ) or PATHOLOGY_DIMENSIONS
    baseline_distance = _dimension_distance(
        source.dimensions,
        target.dimensions,
        scored_dimensions,
    )
    resulting_state: dict[str, float] = {}
    topology_delta: dict[str, float] = {}

    for dimension in PATHOLOGY_DIMENSIONS:
        original = float(source.dimensions[dimension])
        shifted = _clamp(
            original + (float(vector.deltas.get(dimension, 0.0)) * delta_scale),
            TOPOLOGY_PRESSURE_MIN,
            TOPOLOGY_PRESSURE_MAX,
        )
        resulting_state[dimension] = _rounded(shifted)
        topology_delta[dimension] = _rounded(shifted - original)

    resulting_distance = _dimension_distance(
        resulting_state,
        target.dimensions,
        scored_dimensions,
    )
    pathology_reduction_score = _rounded(
        _clamp((baseline_distance - resulting_distance) / baseline_distance, -1.0, 1.0)
    )
    collapse_cost_score = _collapse_cost_score_from_base(
        vector.collapse_cost * collapse_scale,
        topology_delta,
    )
    net_topology_health_score = _rounded(
        _clamp(pathology_reduction_score - (collapse_cost_score * 0.55), -1.0, 1.0)
    )

    improved_dimensions = tuple(
        dimension
        for dimension in PATHOLOGY_DIMENSIONS
        if abs(resulting_state[dimension] - target.dimensions[dimension])
        < abs(source.dimensions[dimension] - target.dimensions[dimension])
    )
    worsened_dimensions = tuple(
        dimension
        for dimension in PATHOLOGY_DIMENSIONS
        if abs(resulting_state[dimension] - target.dimensions[dimension])
        > abs(source.dimensions[dimension] - target.dimensions[dimension])
    )

    return (
        resulting_state,
        topology_delta,
        scored_dimensions,
        improved_dimensions,
        worsened_dimensions,
        pathology_reduction_score,
        collapse_cost_score,
        net_topology_health_score,
    )


def classify_fit(pathology_reduction_score: float, collapse_cost_score: float) -> str:
    if pathology_reduction_score < 0.0:
        return "destabilizing"
    if pathology_reduction_score >= 0.35 and collapse_cost_score <= 0.25:
        return "balanced_normalization_candidate"
    if pathology_reduction_score >= 0.35:
        return "partial_correction_high_collapse_cost"
    if pathology_reduction_score >= 0.15:
        return "weak_partial_correction"
    return "weak_or_unclear_correction"


def apply_mechanism_vector(
    source: TopologyProfile,
    target: TopologyProfile,
    vector: MechanismVector,
) -> PerturbationResult:
    (
        resulting_state,
        topology_delta,
        scored_dimensions,
        improved_dimensions,
        worsened_dimensions,
        pathology_reduction_score,
        collapse_cost_score,
        net_topology_health_score,
    ) = _score_scaled_vector(
        source,
        target,
        vector,
        delta_scale=1.0,
        collapse_scale=1.0,
    )

    evidence_weight = _effective_evidence_weight(vector)
    uncertainty_radius = _effective_uncertainty_radius(vector)
    uncertainty_scales = (
        max(0.0, 1.0 - uncertainty_radius),
        1.0,
        1.0 + uncertainty_radius,
    )
    interval_scores = tuple(
        _score_scaled_vector(
            source,
            target,
            vector,
            delta_scale=scale,
            collapse_scale=scale,
        )
        for scale in uncertainty_scales
    )

    return PerturbationResult(
        mechanism_id=vector.mechanism_id,
        mechanism_family=vector.mechanism_family,
        receptor_profile=vector.receptor_profile,
        abstract_compound_class=vector.abstract_compound_class,
        protein_family=vector.protein_family,
        protein_mechanism_class=vector.protein_mechanism_class,
        protein_state_shift=vector.protein_state_shift,
        pathway_network_perturbation=vector.pathway_network_perturbation,
        evidence_stage=vector.evidence_stage,
        evidence_weight=evidence_weight,
        uncertainty_radius=uncertainty_radius,
        primary_evidence_sources=vector.primary_evidence_sources,
        evidence_refs=vector.evidence_refs,
        calibration_blockers=vector.calibration_blockers,
        confidence_interval_kind=vector.confidence_interval_kind,
        assumption_notes=vector.assumption_notes,
        resulting_state=resulting_state,
        topology_delta=topology_delta,
        scored_dimensions=scored_dimensions,
        improved_dimensions=improved_dimensions,
        worsened_dimensions=worsened_dimensions,
        pathology_reduction_score=pathology_reduction_score,
        pathology_reduction_interval=_interval(
            tuple(score[5] for score in interval_scores)
        ),
        collapse_cost_score=collapse_cost_score,
        collapse_cost_interval=_interval(tuple(score[6] for score in interval_scores)),
        net_topology_health_score=net_topology_health_score,
        net_topology_health_interval=_interval(
            tuple(score[7] for score in interval_scores)
        ),
        fit_label=classify_fit(pathology_reduction_score, collapse_cost_score),
        evidence_readiness_label=_evidence_readiness_label(evidence_weight),
    )


def rank_perturbation_results(
    results: Sequence[PerturbationResult],
) -> tuple[PerturbationResult, ...]:
    return tuple(
        sorted(
            results,
            key=lambda item: (
                item.net_topology_health_score,
                item.pathology_reduction_score,
                -item.collapse_cost_score,
                item.mechanism_id,
            ),
            reverse=True,
        )
    )


def build_calibration_readiness_report(
    results: Sequence[PerturbationResult],
) -> CalibrationReadinessReport:
    vector_count = len(results)
    evidence_backed_vectors = sum(
        1 for result in results if result.evidence_weight > 0.0 and result.evidence_refs
    )
    uncalibrated_vectors = sum(
        1
        for result in results
        if result.evidence_readiness_label == "uncalibrated_hypothesis"
    )
    mean_evidence_weight = _rounded(
        sum(result.evidence_weight for result in results) / vector_count
        if vector_count
        else 0.0
    )
    mean_uncertainty_radius = _rounded(
        sum(result.uncertainty_radius for result in results) / vector_count
        if vector_count
        else 0.0
    )
    if evidence_backed_vectors == 0:
        calibration_status = "uncalibrated_hypothesis_workbench"
    elif evidence_backed_vectors < vector_count:
        calibration_status = "partially_evidence_backed_workbench"
    else:
        calibration_status = "evidence_backed_research_workbench"

    return CalibrationReadinessReport(
        calibration_status=calibration_status,
        practical_use="bounded_hypothesis_comparison_and_falsification",
        mechanism_vectors_reviewed=vector_count,
        evidence_backed_vectors=evidence_backed_vectors,
        uncalibrated_vectors=uncalibrated_vectors,
        mean_evidence_weight=mean_evidence_weight,
        mean_uncertainty_radius=mean_uncertainty_radius,
        clinical_use_allowed=False,
        external_validation_required=True,
        blockers=(
            "no_external_calibration_dataset",
            "no_empirical_delta_estimates",
            "no_patient_specific_inference_allowed",
            "no_brand_name_or_treatment_mapping",
        ),
        next_steps=(
            "attach_source_references_to_each_mechanism_vector",
            "replace_point_deltas_with_evidence_derived_ranges",
            "compare_rankings_against_external_non_patient_datasets",
            "track_assumption_changes_across_calibration_runs",
        ),
    )


def build_pharmacotopology_review(
    *,
    source: TopologyProfile = DEFAULT_SCHIZOPHRENIA_LIKE_PROFILE,
    target: TopologyProfile = DEFAULT_NORMAL_BOUNDED_PROFILE,
    mechanisms: Sequence[MechanismVector] = DEFAULT_MECHANISM_VECTORS,
) -> dict[str, Any]:
    results = rank_perturbation_results(
        tuple(apply_mechanism_vector(source, target, vector) for vector in mechanisms)
    )
    top_result = results[0]
    calibration_readiness = build_calibration_readiness_report(results)

    return {
        "Φ.kind": "Φ.clean.pharmacotopology_layer.v0",
        "Φ.scope": {
            "simulation_only": True,
            "hypothesis_numbers_only": True,
            "clinical_advice_created": False,
            "medication_recommendation_created": False,
            "real_patient_inference_created": False,
            "brand_name_mapping_created": False,
            "dimension_direction": "lower_pressure_is_closer_to_bounded_review",
        },
        "Φ.practical_use": {
            "research_use_label": "bounded_hypothesis_workbench",
            "allowed_use": "compare_and_falsify_mechanism_assumptions",
            "clinical_use_allowed": False,
            "requires_external_calibration_for_claims": True,
        },
        "Φ.source_profile": asdict(source),
        "Φ.target_profile": asdict(target),
        "Φ.mechanism_vectors": [asdict(vector) for vector in mechanisms],
        "Φ.results": [asdict(result) for result in results],
        "Φ.ranking": [
            {
                "rank": index,
                "mechanism_id": result.mechanism_id,
                "pathology_reduction_score": result.pathology_reduction_score,
                "collapse_cost_score": result.collapse_cost_score,
                "net_topology_health_score": result.net_topology_health_score,
                "net_topology_health_interval": result.net_topology_health_interval,
                "evidence_weight": result.evidence_weight,
                "uncertainty_radius": result.uncertainty_radius,
                "evidence_readiness_label": result.evidence_readiness_label,
                "fit_label": result.fit_label,
            }
            for index, result in enumerate(results, start=1)
        ],
        "Φ.review": {
            "valid": True,
            "source_profile_seen": True,
            "target_profile_seen": True,
            "mechanism_vectors_reviewed": len(mechanisms),
            "topology_dimensions_reviewed": len(PATHOLOGY_DIMENSIONS),
            "pathology_reduction_score_seen": True,
            "collapse_cost_score_seen": True,
            "net_topology_health_score_seen": True,
            "top_mechanism_id": top_result.mechanism_id,
            "top_net_topology_health_score": top_result.net_topology_health_score,
            "destabilizing_mechanism_count": sum(
                1 for result in results if result.fit_label == "destabilizing"
            ),
        },
        "Φ.calibration_readiness": asdict(calibration_readiness),
        "Φ.calibration": {
            "stage": "hypothesis_only",
            "future_sources": [
                "clinical_trial_outcomes",
                "side_effect_profiles",
                "receptor_binding_data",
                "patient_symptom_scales",
                "relapse_rates",
                "negative_symptom_data",
                "cognition_data",
                "real_world_discontinuation_rates",
            ],
        },
        "Φ.claim": {
            "clinical_advice_created": False,
            "medication_recommendation_created": False,
            "real_patient_inference_created": False,
            "brand_name_mapping_created": False,
            "treatment_claim_created": False,
            "cure_claim_created": False,
        },
        "Φ.neg": {
            "voice_opened": False,
            "schedule_created": False,
            "command_queue_created": False,
            "self_practice_created": False,
            "autonomous_loop_created": False,
            "repeat_entitlement_created": False,
        },
        "Φ.doctrine": {
            "mechanism_vector_is_medication_advice": False,
            "topology_delta_is_treatment_claim": False,
            "pathology_reduction_is_cure_claim": False,
            "collapse_cost_is_side_effect_prediction": False,
            "simulated_profile_is_patient_profile": False,
            "ranking_is_prescribing_guidance": False,
        },
    }


def build_pharmacotopology_packet(
    input_atom: dict[str, object],
    *,
    review: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    return {
        FieldKey.SCHEMA: "κ.clean.pharmacotopology_layer.1",
        FieldKey.CYCLE: 1,
        FieldKey.CONTEXT: "χ.clean.pharmacotopology_layer",
        FieldKey.INPUT: input_atom,
        FieldKey.PRIOR_COUNT: 0,
        FieldKey.PRIOR_TRACE: None,
        FieldKey.ACTION: FieldAction.WRITE_MEMORY,
        FieldKey.MEMORY: "μ.clean.pharmacotopology_layer.Φ.τ1.ρ0",
        FieldKey.BOUNDARY: list(FIELD_CLOSED_CHANNELS),
        FieldKey.RESPONSE: FieldResponse.READOUT,
        FieldKey.PHARMACOTOPOLOGY_REVIEW: review or build_pharmacotopology_review(),
        FieldKey.NEXT: "ν.claim_boundary",
        FieldKey.STOP: "ω.clean",
    }


def report_from_review(
    review: dict[str, Any],
) -> CleanPharmacotopologyLayerReport:
    review_flags = review["Φ.review"]
    claim = review["Φ.claim"]
    neg = review["Φ.neg"]
    ranking = review["Φ.ranking"]
    calibration = review["Φ.calibration_readiness"]
    top = ranking[0] if ranking else {}

    return CleanPharmacotopologyLayerReport(
        run_kind="clean_pharmacotopology_layer",
        pharmacotopology_review_valid=bool(review_flags["valid"]),
        simulation_only=bool(review["Φ.scope"]["simulation_only"]),
        hypothesis_numbers_only=bool(review["Φ.scope"]["hypothesis_numbers_only"]),
        calibration_status=str(calibration["calibration_status"]),
        practical_use=str(calibration["practical_use"]),
        evidence_backed_vectors=int(calibration["evidence_backed_vectors"]),
        uncalibrated_vectors=int(calibration["uncalibrated_vectors"]),
        mean_evidence_weight=float(calibration["mean_evidence_weight"]),
        mean_uncertainty_radius=float(calibration["mean_uncertainty_radius"]),
        clinical_use_allowed=bool(calibration["clinical_use_allowed"]),
        mechanism_vectors_reviewed=int(review_flags["mechanism_vectors_reviewed"]),
        topology_dimensions_reviewed=int(review_flags["topology_dimensions_reviewed"]),
        top_mechanism_id=str(review_flags["top_mechanism_id"]),
        top_net_topology_health_score=float(top["net_topology_health_score"]),
        top_pathology_reduction_score=float(top["pathology_reduction_score"]),
        top_collapse_cost_score=float(top["collapse_cost_score"]),
        destabilizing_mechanism_count=int(
            review_flags["destabilizing_mechanism_count"]
        ),
        clinical_advice_created=bool(claim["clinical_advice_created"]),
        medication_recommendation_created=bool(
            claim["medication_recommendation_created"]
        ),
        real_patient_inference_created=bool(claim["real_patient_inference_created"]),
        brand_name_mapping_created=bool(claim["brand_name_mapping_created"]),
        treatment_claim_created=bool(claim["treatment_claim_created"]),
        voice_opened=bool(neg["voice_opened"]),
        schedule_created=bool(neg["schedule_created"]),
        command_queue_created=bool(neg["command_queue_created"]),
        autonomous_loop_created=bool(neg["autonomous_loop_created"]),
        surface_leakage=0.0,
        stop_integrity=1.0,
    )


def refusal_report() -> CleanPharmacotopologyLayerReport:
    return CleanPharmacotopologyLayerReport(
        run_kind="clean_pharmacotopology_layer",
        pharmacotopology_review_valid=False,
        simulation_only=True,
        hypothesis_numbers_only=True,
        calibration_status="refused",
        practical_use="none",
        evidence_backed_vectors=0,
        uncalibrated_vectors=0,
        mean_evidence_weight=0.0,
        mean_uncertainty_radius=0.0,
        clinical_use_allowed=False,
        mechanism_vectors_reviewed=0,
        topology_dimensions_reviewed=0,
        top_mechanism_id="",
        top_net_topology_health_score=0.0,
        top_pathology_reduction_score=0.0,
        top_collapse_cost_score=0.0,
        destabilizing_mechanism_count=0,
        clinical_advice_created=False,
        medication_recommendation_created=False,
        real_patient_inference_created=False,
        brand_name_mapping_created=False,
        treatment_claim_created=False,
        voice_opened=False,
        schedule_created=False,
        command_queue_created=False,
        autonomous_loop_created=False,
        surface_leakage=0.0,
        stop_integrity=1.0,
    )


def write_report(
    run_dir: Path,
    report: CleanPharmacotopologyLayerReport,
) -> Path:
    path = run_dir / "clean_pharmacotopology_layer_report.json"
    path.write_text(
        json.dumps(asdict(report), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_calibration_readiness_report(
    run_dir: Path,
    review: dict[str, Any],
) -> Path:
    path = run_dir / "calibration_readiness_report.json"
    path.write_text(
        json.dumps(
            review["Φ.calibration_readiness"],
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def append_packet(run_dir: Path, packet: dict[str, Any]) -> None:
    event_group(run_dir, cycle=1)
    append_jsonl(
        run_dir / "memory.jsonl",
        {
            "cycle": 1,
            "content": json.dumps(packet, ensure_ascii=False, sort_keys=True),
            "created_at": 0.0,
        },
    )
    append_jsonl(
        run_dir / "session_records.jsonl",
        asdict(
            ManualSessionRecord(
                session_id="clean.pharmacotopology_layer.0",
                turn_index=1,
                iota=packet[FieldKey.INPUT],
                psi=FieldResponse.READOUT,
                mu=str(packet[FieldKey.MEMORY]),
                rho_n=0,
                boundary_status=boundary_status_dict(),
                omega="ω.clean",
                created_at=now_timestamp(),
            )
        ),
    )


def run_clean_pharmacotopology_layer(
    run_dir: Path,
    surface: str = DEFAULT_PHARMACOTOPOLOGY_SURFACE,
    source: TopologyProfile = DEFAULT_SCHIZOPHRENIA_LIKE_PROFILE,
    target: TopologyProfile = DEFAULT_NORMAL_BOUNDED_PROFILE,
    mechanisms: Sequence[MechanismVector] = DEFAULT_MECHANISM_VECTORS,
) -> CleanPharmacotopologyLayerReport:
    reset_pharmacotopology_run_dir(run_dir)
    quarantine = InputQuarantine(run_dir / "input_stream.jsonl")
    atom = quarantine.ingest(surface, kind="operator")

    certificate = create_evidence_sealed_certificate()
    if not certificate.clean_transfer_ready():
        append_jsonl(
            run_dir / "audit.jsonl",
            {
                "cycle": 0,
                "kind": FieldEvent.DENIED,
                "message": FieldStatus.CONTRACT_REFUSED,
            },
        )
        append_jsonl(
            run_dir / "audit.jsonl",
            {
                "cycle": 0,
                "kind": FieldEvent.STOPPED,
                "message": FieldStatus.EXTERNAL_STOP,
            },
        )
        report = refusal_report()
        write_report(run_dir, report)
        return report

    review = build_pharmacotopology_review(
        source=source,
        target=target,
        mechanisms=mechanisms,
    )
    packet = build_pharmacotopology_packet(atom.field_packet(), review=review)
    append_packet(run_dir, packet)
    report = report_from_review(packet[FieldKey.PHARMACOTOPOLOGY_REVIEW])
    write_calibration_readiness_report(
        run_dir,
        packet[FieldKey.PHARMACOTOPOLOGY_REVIEW],
    )
    write_report(run_dir, report)
    return report
