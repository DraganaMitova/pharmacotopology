from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

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
    "pharmacotopology_dashboard.html",
    "field_validation.json",
    "field_metrics.json",
)


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
    deltas: dict[str, float]
    collapse_cost: float
    evidence_stage: str = "hypothesis_only"


@dataclass(frozen=True)
class PerturbationResult:
    mechanism_id: str
    mechanism_family: str
    receptor_profile: tuple[str, ...]
    evidence_stage: str
    resulting_state: dict[str, float]
    topology_delta: dict[str, float]
    scored_dimensions: tuple[str, ...]
    improved_dimensions: tuple[str, ...]
    worsened_dimensions: tuple[str, ...]
    pathology_reduction_score: float
    collapse_cost_score: float
    net_topology_health_score: float
    fit_label: str


@dataclass(frozen=True)
class CleanPharmacotopologyLayerReport:
    run_kind: str
    pharmacotopology_review_valid: bool
    simulation_only: bool
    hypothesis_numbers_only: bool
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

DEFAULT_MECHANISM_VECTORS: tuple[MechanismVector, ...] = (
    MechanismVector(
        mechanism_id="d2_antagonist_like",
        mechanism_family="dopamine_salience_dampening",
        receptor_profile=("D2 antagonism",),
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
    penalty = (
        max(float(topology_delta.get("negative_shutdown", 0.0)), 0.0) * 0.50
        + max(float(topology_delta.get("cognitive_fragmentation", 0.0)), 0.0)
        * 0.35
        + max(float(topology_delta.get("boundary_instability", 0.0)), 0.0)
        * 0.20
    )
    return _rounded(_clamp(vector.collapse_cost + penalty, 0.0, 1.0))


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
            original + float(vector.deltas.get(dimension, 0.0)),
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
    collapse_cost_score = _collapse_cost_score(vector, topology_delta)
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

    return PerturbationResult(
        mechanism_id=vector.mechanism_id,
        mechanism_family=vector.mechanism_family,
        receptor_profile=vector.receptor_profile,
        evidence_stage=vector.evidence_stage,
        resulting_state=resulting_state,
        topology_delta=topology_delta,
        scored_dimensions=scored_dimensions,
        improved_dimensions=improved_dimensions,
        worsened_dimensions=worsened_dimensions,
        pathology_reduction_score=pathology_reduction_score,
        collapse_cost_score=collapse_cost_score,
        net_topology_health_score=net_topology_health_score,
        fit_label=classify_fit(pathology_reduction_score, collapse_cost_score),
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


def build_pharmacotopology_packet(input_atom: dict[str, object]) -> dict[str, Any]:
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
        FieldKey.PHARMACOTOPOLOGY_REVIEW: build_pharmacotopology_review(),
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
    top = ranking[0] if ranking else {}

    return CleanPharmacotopologyLayerReport(
        run_kind="clean_pharmacotopology_layer",
        pharmacotopology_review_valid=bool(review_flags["valid"]),
        simulation_only=bool(review["Φ.scope"]["simulation_only"]),
        hypothesis_numbers_only=bool(review["Φ.scope"]["hypothesis_numbers_only"]),
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

    packet = build_pharmacotopology_packet(atom.field_packet())
    append_packet(run_dir, packet)
    report = report_from_review(packet[FieldKey.PHARMACOTOPOLOGY_REVIEW])
    write_report(run_dir, report)
    return report
