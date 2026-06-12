from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V87_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V87"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v87_fills_real_input_and_execution_zeroes_without_manual_solved_flag() -> None:
    cert = _read(V87_ROOT / "v87_real_inputs_and_execution_fill_campaign_certificate.json")

    assert cert["kind"] == "V87_REAL_INPUTS_AND_EXECUTION_FILL_CAMPAIGN_CERTIFICATE_v0"
    assert cert["status"] == "V87_REAL_INPUTS_AND_EXECUTION_FILL_PARTIAL_UNIVERSAL_BLOCKED"
    assert cert["campaign_controls_passed"] is True
    assert cert["engine_version_used"] == "E80"
    assert cert["source_batch_id"] == "V86_TRUE_FOLD_SOLUTION_CLOSURE_CAMPAIGN"
    assert cert["fresh_target_count"] > 0
    assert cert["coordinate_emission_target_count"] > 0
    assert cert["real_fold_holdout_count"] > 0
    assert cert["real_or_validated_physical_execution_count"] > 0
    assert cert["real_openmm_execution_count"] > 0
    assert cert["proxy_physical_execution_used_for_claim"] is False
    assert cert["target_fold_claim_count"] > 0
    assert cert["unsupported_fold_claims"] == 0
    assert cert["unsupported_physical_claims"] == 0
    assert cert["coordinate_native_leakage"] is False
    assert cert["token_only_acceptance_count"] == 0
    assert cert["external_blind_benchmark_exported"] is True
    assert cert["external_blind_benchmark_passed"] is True
    assert cert["fresh_target_shortage"] is True
    assert cert["every_required_hard_family_has_supported_target_fold_claim"] is False
    assert cert["supported_hard_families"] == ["globular_soluble"]
    assert cert["universal_folding_solution_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["target_fold_claim_count_nonzero_full_suite_required"] is True
    assert cert["failed_controls"] == []
    assert cert["blocked_reasons"] == [
        "fresh_target_shortage",
        "hard_family_support_incomplete",
    ]


def test_v87_report_records_sequence_only_hydration_postseal_holdouts_and_openmm() -> None:
    report = _read(V87_ROOT / "v87_real_inputs_and_execution_fill_campaign_report.json")
    export = _read(V87_ROOT / "v87_external_blind_benchmark_export.json")
    sequence_manifest = _read(V87_ROOT / "v87_sequence_only_target_hydration_manifest.json")

    assert report["kind"] == "V87_REAL_INPUTS_AND_EXECUTION_FILL_CAMPAIGN_REPORT_v0"
    assert sequence_manifest["kind"] == "V87_SEQUENCE_ONLY_TARGET_HYDRATION_MANIFEST_v0"
    assert sequence_manifest["coordinates_opened_during_sequence_hydration"] is False
    assert sequence_manifest["native_contacts_opened_during_sequence_hydration"] is False
    assert sequence_manifest["target_count"] > 0
    assert all(target["sequence_hydration_source"] for target in sequence_manifest["targets"])
    assert len(report["predicted_fold_packets"]) == sequence_manifest["target_count"]
    assert all(packet["prediction_hash"] for packet in report["predicted_fold_packets"])
    assert len(report["real_holdout_rows"]) == sequence_manifest["target_count"]
    assert all(row["opened_after_prediction_hash"] is True for row in report["real_holdout_rows"])
    assert all(row["coordinate_truth_used_before_prediction"] is False for row in report["real_holdout_rows"])
    assert all(row["native_contacts_used_before_prediction"] is False for row in report["real_holdout_rows"])
    assert len(report["physical_execution_rows"]) == sequence_manifest["target_count"]
    assert all(row["execution_backend"] == "openmm" for row in report["physical_execution_rows"])
    assert all(row["real_openmm_execution"] is True for row in report["physical_execution_rows"])
    assert all(row["proxy_only"] is False for row in report["physical_execution_rows"])
    assert report["universal_solution_unlock_firewall"]["protein_folding_solved"] is False
    assert export["kind"] == "V87_EXTERNAL_BLIND_BENCHMARK_EXPORT_v0"
    assert export["rows"]
