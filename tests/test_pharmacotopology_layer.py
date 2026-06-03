from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.field_alphabet import FieldKey
from pharmacotopology.layer import (
    DEFAULT_MECHANISM_VECTORS,
    DEFAULT_NORMAL_BOUNDED_PROFILE,
    DEFAULT_SCHIZOPHRENIA_LIKE_PROFILE,
    DEFAULT_TOPOLOGY_PROFILES,
    PATHOLOGY_DIMENSIONS,
    apply_mechanism_vector,
    build_pharmacotopology_review,
    get_topology_profile,
    rank_perturbation_results,
    run_clean_pharmacotopology_layer,
)


def test_pharmacotopology_ranking_separates_helpful_and_destabilizing() -> None:
    results = rank_perturbation_results(
        tuple(
            apply_mechanism_vector(
                DEFAULT_SCHIZOPHRENIA_LIKE_PROFILE,
                DEFAULT_NORMAL_BOUNDED_PROFILE,
                vector,
            )
            for vector in DEFAULT_MECHANISM_VECTORS
        )
    )

    top = results[0]
    destabilizing = next(
        result
        for result in results
        if result.mechanism_id == "glutamate_amplification_stressor_like"
    )

    assert top.net_topology_health_score > 0.0
    assert destabilizing.fit_label == "destabilizing"
    assert destabilizing.pathology_reduction_score < 0.0
    assert top.net_topology_health_score > destabilizing.net_topology_health_score


def test_topology_profile_catalog_exposes_synthetic_profiles() -> None:
    assert set(DEFAULT_TOPOLOGY_PROFILES) == {
        "schizophrenia_like",
        "depression_like",
        "mania_like",
        "anxiety_like",
        "mixed_state_like",
    }

    for profile in DEFAULT_TOPOLOGY_PROFILES.values():
        assert tuple(profile.dimensions) == PATHOLOGY_DIMENSIONS
        assert "not a diagnosis" in profile.description

    assert get_topology_profile("mania_like").profile_id == (
        "mania_like_topology_profile"
    )


def test_pharmacotopology_review_denies_clinical_claims() -> None:
    review = build_pharmacotopology_review()

    assert review["Φ.review"]["valid"] is True
    assert review["Φ.scope"]["simulation_only"] is True
    assert review["Φ.scope"]["hypothesis_numbers_only"] is True
    assert review["Φ.claim"]["clinical_advice_created"] is False
    assert review["Φ.claim"]["medication_recommendation_created"] is False
    assert review["Φ.claim"]["real_patient_inference_created"] is False
    assert review["Φ.claim"]["brand_name_mapping_created"] is False
    assert review["Φ.doctrine"]["ranking_is_prescribing_guidance"] is False


def test_pharmacotopology_review_exposes_calibration_readiness() -> None:
    review = build_pharmacotopology_review()
    calibration = review["Φ.calibration_readiness"]
    top_result = review["Φ.results"][0]

    assert review["Φ.practical_use"]["research_use_label"] == (
        "bounded_hypothesis_workbench"
    )
    assert review["Φ.practical_use"]["clinical_use_allowed"] is False
    assert calibration["practical_use"] == (
        "bounded_hypothesis_comparison_and_falsification"
    )
    assert calibration["clinical_use_allowed"] is False
    assert calibration["external_validation_required"] is True
    assert calibration["calibration_status"] == "uncalibrated_hypothesis_workbench"
    assert top_result["evidence_readiness_label"] == "uncalibrated_hypothesis"
    assert top_result["net_topology_health_interval"]["lower"] <= (
        top_result["net_topology_health_score"]
    )
    assert top_result["net_topology_health_interval"]["upper"] >= (
        top_result["net_topology_health_score"]
    )


def test_clean_pharmacotopology_run_writes_bounded_artifacts(tmp_path: Path) -> None:
    report = run_clean_pharmacotopology_layer(tmp_path)

    assert report.pharmacotopology_review_valid is True
    assert report.mechanism_vectors_reviewed == len(DEFAULT_MECHANISM_VECTORS)
    assert report.topology_dimensions_reviewed > 0
    assert report.clinical_advice_created is False
    assert report.medication_recommendation_created is False
    assert report.voice_opened is False
    assert report.stop_integrity == 1.0
    assert report.calibration_status == "uncalibrated_hypothesis_workbench"
    assert report.practical_use == "bounded_hypothesis_comparison_and_falsification"
    assert report.clinical_use_allowed is False

    calibration_report = json.loads(
        (tmp_path / "calibration_readiness_report.json").read_text(encoding="utf-8")
    )
    assert calibration_report["clinical_use_allowed"] is False
    assert calibration_report["external_validation_required"] is True

    rows = [
        json.loads(line)
        for line in (tmp_path / "memory.jsonl").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    packet = json.loads(rows[0]["content"])
    assert FieldKey.PHARMACOTOPOLOGY_REVIEW in packet
    assert packet[FieldKey.RESPONSE] == "ψ.readout"
    assert packet[FieldKey.PHARMACOTOPOLOGY_REVIEW]["Φ.scope"][
        "clinical_advice_created"
    ] is False


def test_clean_pharmacotopology_run_accepts_synthetic_profile(
    tmp_path: Path,
) -> None:
    run_clean_pharmacotopology_layer(
        tmp_path,
        source=get_topology_profile("anxiety_like"),
    )

    rows = [
        json.loads(line)
        for line in (tmp_path / "memory.jsonl").read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    packet = json.loads(rows[0]["content"])
    review = packet[FieldKey.PHARMACOTOPOLOGY_REVIEW]

    assert review["Φ.source_profile"]["profile_id"] == "anxiety_like_topology_profile"
    assert review["Φ.claim"]["real_patient_inference_created"] is False
