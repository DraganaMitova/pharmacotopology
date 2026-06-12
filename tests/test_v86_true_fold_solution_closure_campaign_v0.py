from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V86_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V86"
E80_ROOT = ROOT / "data" / "protein_esperanto_engine" / "E80"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v86_certificate_blocks_without_real_closure_inputs() -> None:
    cert = _read(V86_ROOT / "v86_true_fold_solution_closure_campaign_certificate.json")

    assert cert["kind"] == "V86_TRUE_FOLD_SOLUTION_CLOSURE_CAMPAIGN_CERTIFICATE_v0"
    assert cert["status"] == "V86_TRUE_FOLD_SOLUTION_CLOSURE_BLOCKED"
    assert cert["campaign_controls_passed"] is True
    assert cert["engine_version_used"] == "E80"
    assert cert["baseline_engine_version"] == "E79"
    assert cert["source_batch_id"] == "V85_TRUE_UNIVERSAL_FOLDING_UNLOCK_CAMPAIGN"
    assert cert["highest_claim_tier_unlocked"] == "CLAIM_4_COARSE_PHYSICAL_SUPPORTED"
    assert cert["previous_claim_ceiling"]["source_claim_ceiling_inherited_from_v85"] is True
    assert cert["fresh_target_shortage"] is True
    assert cert["fresh_target_count"] == 0
    assert cert["coordinate_emission_target_count"] == 0
    assert cert["fold_constraint_packet_count"] == 0
    assert cert["distance_map_packet_count"] == 0
    assert cert["topology_constraint_packet_count"] == 0
    assert cert["real_fold_holdout_count"] == 0
    assert cert["real_or_validated_physical_execution_count"] == 0
    assert cert["proxy_physical_execution_used_for_claim"] is False
    assert cert["target_fold_claim_count"] == 0
    assert cert["every_required_hard_family_has_supported_target_fold_claim"] is False
    assert cert["unsupported_fold_claims"] == 0
    assert cert["unsupported_physical_claims"] == 0
    assert cert["coordinate_native_leakage"] is False
    assert cert["external_fold_model_prediction_input_used"] is False
    assert cert["token_only_acceptance_count"] == 0
    assert cert["external_blind_benchmark_exported"] is True
    assert cert["external_blind_benchmark_passed"] is False
    assert cert["general_solution_candidate_claim_allowed"] is False
    assert cert["universal_folding_solution_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["failed_controls"] == []
    assert "fresh_target_shortage" in cert["blocked_reasons"]
    assert "real_fold_holdout_count_zero" in cert["blocked_reasons"]
    assert "real_or_validated_physical_execution_count_zero" in cert["blocked_reasons"]
    assert "target_fold_claim_count_zero" in cert["blocked_reasons"]
    assert "hard_family_support_incomplete" in cert["blocked_reasons"]
    assert "external_blind_benchmark_not_passed" in cert["blocked_reasons"]


def test_v86_report_and_e80_certificate_expose_geometry_bridge() -> None:
    report = _read(V86_ROOT / "v86_true_fold_solution_closure_campaign_report.json")
    export = _read(V86_ROOT / "v86_external_blind_benchmark_export.json")
    e80_cert = _read(E80_ROOT / "e80_real_fold_geometry_and_physical_calibration_engine_certificate.json")

    assert report["kind"] == "V86_TRUE_FOLD_SOLUTION_CLOSURE_CAMPAIGN_REPORT_v0"
    assert set(report["components"]) == {
        "manifest",
        "fresh_target_resolver",
        "fold_geometry_compiler",
        "constraint_to_distance_map_compiler",
        "topology_constraint_compiler",
        "backbone_coordinate_emitter",
        "real_holdout_coordinate_loader",
        "physical_relaxation_executor",
        "fold_quality_evaluator",
        "external_blind_benchmark_export",
        "universal_solution_unlock_firewall",
    }
    assert report["components"]["manifest"]["engine_revision"] == "E80_REAL_FOLD_GEOMETRY_AND_PHYSICAL_CALIBRATION_ENGINE"
    assert report["components"]["manifest"]["external_fold_model_prediction_input_allowed"] is False
    assert report["components"]["manifest"]["proxy_physical_execution_for_target_fold_claim_allowed"] is False
    assert report["components"]["fresh_target_resolver"]["fresh_target_shortage"] is True
    assert report["components"]["fold_geometry_compiler"] == []
    assert report["components"]["backbone_coordinate_emitter"] == []
    assert report["components"]["real_holdout_coordinate_loader"]["real_fold_holdout_count"] == 0
    assert report["components"]["physical_relaxation_executor"]["real_or_validated_physical_execution_count"] == 0
    assert report["components"]["universal_solution_unlock_firewall"]["protein_folding_solved"] is False
    assert report["real_closure_input_files"]["fresh_targets_manifest"]["exists"] is False
    assert report["real_closure_input_files"]["physical_execution_rows"]["exists"] is False

    assert export["kind"] == "V86_EXTERNAL_BLIND_BENCHMARK_EXPORT_v0"
    assert export["rows"] == []
    assert export["empty_export_reason"] == "no_target_fold_claims_allowed_for_external_export"

    assert e80_cert["kind"] == "E80_REAL_FOLD_GEOMETRY_AND_PHYSICAL_CALIBRATION_ENGINE_CERTIFICATE_v0"
    assert e80_cert["engine_revision"] == "E80"
    assert e80_cert["baseline_engine_revision"] == "E79"
    assert e80_cert["protein_folding_solved_by_field_setting"] is False
    assert e80_cert["universal_folding_solution_claim_allowed"] is False
    assert e80_cert["protein_folding_solved"] is False
