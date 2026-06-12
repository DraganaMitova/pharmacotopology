from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V85_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V85"
E79_ROOT = ROOT / "data" / "protein_esperanto_engine" / "E79"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v85_certificate_blocks_universal_claim_without_fresh_real_unlock_evidence() -> None:
    cert = _read(V85_ROOT / "v85_true_universal_folding_unlock_campaign_certificate.json")

    assert cert["kind"] == "V85_TRUE_UNIVERSAL_FOLDING_UNLOCK_CAMPAIGN_CERTIFICATE_v0"
    assert cert["status"] == "V85_TRUE_UNIVERSAL_FOLDING_UNLOCK_BLOCKED"
    assert cert["campaign_controls_passed"] is True
    assert cert["engine_version_used"] == "E79"
    assert cert["baseline_engine_version"] == "E78"
    assert cert["source_batch_id"] == "V84_FINAL_CLOSED_LOOP_PROTEIN_FOLDING_CAMPAIGN"
    assert cert["highest_claim_tier_unlocked"] == "CLAIM_4_COARSE_PHYSICAL_SUPPORTED"
    assert cert["previous_claim_ceiling"]["source_claim_ceiling_inherited_from_v84"] is True
    assert cert["fresh_target_shortage"] is True
    assert cert["fresh_target_count"] == 0
    assert cert["real_fold_holdout_count"] == 0
    assert cert["real_or_validated_physical_execution_count"] == 0
    assert cert["proxy_physical_execution_used_for_claim"] is False
    assert cert["target_fold_claim_count"] == 0
    assert cert["general_solution_candidate_claim_allowed"] is False
    assert cert["universal_folding_solution_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["unsupported_fold_claims"] == 0
    assert cert["unsupported_physical_claims"] == 0
    assert cert["unsupported_claims"] == 0
    assert cert["coordinate_native_leakage"] is False
    assert cert["external_blind_benchmark_exported"] is True
    assert cert["external_blind_benchmark_passed"] is False
    assert cert["prior_calibration_inputs_count"] >= 1
    assert cert["prior_calibration_inputs_used_for_universal_unlock"] is False
    assert cert["deterministic_variants_for_universal_claim_allowed"] is False
    assert cert["failed_accepted_count"] == 0
    assert cert["failed_controls"] == []
    assert "fresh_target_shortage" in cert["blocked_reasons"]
    assert "no_real_or_validated_physical_execution" in cert["blocked_reasons"]
    assert "target_fold_claim_count_zero" in cert["blocked_reasons"]
    assert "family_generalization_not_sufficient" in cert["blocked_reasons"]
    assert "external_blind_benchmark_not_passed" in cert["blocked_reasons"]


def test_v85_report_and_e79_certificate_expose_unlock_components() -> None:
    report = _read(V85_ROOT / "v85_true_universal_folding_unlock_campaign_report.json")
    export = _read(V85_ROOT / "v85_external_blind_benchmark_export.json")
    e79_cert = _read(E79_ROOT / "e79_universal_folding_claim_unlock_engine_certificate.json")

    assert report["kind"] == "V85_TRUE_UNIVERSAL_FOLDING_UNLOCK_CAMPAIGN_REPORT_v0"
    assert set(report["components"]) == {
        "manifest",
        "fresh_target_resolver",
        "real_fold_holdout_loader",
        "atomistic_or_validated_physical_executor",
        "target_fold_evaluator",
        "family_generalization_evaluator",
        "external_blind_benchmark_export",
        "universal_claim_firewall",
    }
    assert report["components"]["manifest"]["fresh_targets_required_for_universal_claim"] is True
    assert report["components"]["manifest"]["proxy_physical_execution_for_claim_allowed"] is False
    assert report["components"]["fresh_target_resolver"]["fresh_target_shortage"] is True
    assert report["components"]["atomistic_or_validated_physical_executor"]["real_or_validated_execution_count"] == 0
    assert report["components"]["target_fold_evaluator"]["target_fold_claim_count"] == 0
    assert report["components"]["universal_claim_firewall"]["protein_folding_solved"] is False
    assert report["real_unlock_input_files"]["fresh_targets_manifest"]["exists"] is False
    assert report["real_unlock_input_files"]["real_physical_executions"]["exists"] is False
    assert report["prior_calibration_inputs"]

    assert export["kind"] == "V85_EXTERNAL_BLIND_BENCHMARK_EXPORT_v0"
    assert export["rows"] == []
    assert export["empty_export_reason"] == "no_target_fold_claims_allowed_for_external_export"

    assert e79_cert["kind"] == "E79_UNIVERSAL_FOLDING_CLAIM_UNLOCK_ENGINE_CERTIFICATE_v0"
    assert e79_cert["engine_revision"] == "E79"
    assert e79_cert["baseline_engine_revision"] == "E78"
    assert e79_cert["universal_folding_solution_claim_allowed_by_field_setting"] is False
    assert e79_cert["protein_folding_solved_by_field_setting"] is False
    assert e79_cert["universal_folding_solution_claim_allowed"] is False
    assert e79_cert["protein_folding_solved"] is False
