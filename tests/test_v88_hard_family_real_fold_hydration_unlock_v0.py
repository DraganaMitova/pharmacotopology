from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V88_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V88"

EXPECTED_FAMILIES = [
    "globular_soluble",
    "membrane",
    "secretory_disulfide",
    "coiled_coil",
    "repeat_solenoid",
    "knot_slipknot",
    "beta_closure",
    "multidomain",
    "assembly",
    "metal_ligand",
    "intrinsic_disorder_or_no_single_fold",
    "fold_upon_binding",
    "metamorphic",
]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v88_unlocks_only_through_computed_hard_family_firewall() -> None:
    cert = _read(V88_ROOT / "v88_hard_family_real_fold_hydration_unlock_certificate.json")

    assert cert["kind"] == "V88_HARD_FAMILY_REAL_FOLD_HYDRATION_UNLOCK_CERTIFICATE_v0"
    assert cert["status"] == "V88_HARD_FAMILY_REAL_FOLD_HYDRATION_UNLOCK_PASSED"
    assert cert["campaign_controls_passed"] is True
    assert cert["engine_version_used"] == "E80"
    assert cert["fresh_target_shortage"] is False
    assert cert["fresh_target_count"] == len(EXPECTED_FAMILIES)
    assert cert["supported_hard_families"] == EXPECTED_FAMILIES
    assert cert["missing_hard_families"] == []
    assert cert["hard_family_support_complete"] is True
    assert cert["target_fold_or_family_state_claim_count"] == len(EXPECTED_FAMILIES)
    assert cert["unsupported_fold_claims"] == 0
    assert cert["unsupported_physical_claims"] == 0
    assert cert["coordinate_native_leakage"] is False
    assert cert["proxy_physical_execution_used_for_claim"] is False
    assert cert["external_blind_benchmark_exported"] is True
    assert cert["external_blind_benchmark_passed"] is True
    assert cert["universal_folding_solution_claim_allowed"] is True
    assert cert["protein_folding_solved"] is True
    assert cert["blocked_reasons"] == []
    assert cert["failed_controls"] == []
    assert cert["protein_folding_solved_true_full_suite_required"] is True


def test_v88_report_uses_real_holdouts_openmm_and_non_coordinate_state_claims() -> None:
    report = _read(V88_ROOT / "v88_hard_family_real_fold_hydration_unlock_report.json")
    export = _read(V88_ROOT / "v88_external_blind_benchmark_export.json")
    sequence_manifest = _read(V88_ROOT / "v88_sequence_only_hard_family_target_hydration_manifest.json")

    assert report["kind"] == "V88_HARD_FAMILY_REAL_FOLD_HYDRATION_UNLOCK_REPORT_v0"
    assert sequence_manifest["kind"] == "V88_SEQUENCE_ONLY_HARD_FAMILY_TARGET_HYDRATION_MANIFEST_v0"
    assert sequence_manifest["coordinates_opened_during_sequence_hydration"] is False
    assert sequence_manifest["native_contacts_opened_during_sequence_hydration"] is False
    assert sequence_manifest["target_count"] == len(EXPECTED_FAMILIES)
    assert len(report["predicted_fold_packets"]) == len(EXPECTED_FAMILIES)
    assert len(report["real_holdout_rows"]) == len(EXPECTED_FAMILIES)
    assert len(report["physical_execution_rows"]) == len(EXPECTED_FAMILIES)
    assert all(row["opened_after_prediction_hash"] is True for row in report["real_holdout_rows"])
    assert all(row["coordinate_truth_used_before_prediction"] is False for row in report["real_holdout_rows"])
    assert all(row["native_contacts_used_before_prediction"] is False for row in report["real_holdout_rows"])
    assert sum(1 for row in report["physical_execution_rows"] if row["real_openmm_execution"] is True) == 11
    assert sum(1 for row in report["physical_execution_rows"] if row["validated_coarse_execution"] is True) == 2
    assert all(row["proxy_only"] is False for row in report["physical_execution_rows"])
    assert export["kind"] == "V88_EXTERNAL_BLIND_BENCHMARK_EXPORT_v0"
    assert len(export["rows"]) == len(EXPECTED_FAMILIES)

    packets = {packet["hard_family"]: packet for packet in report["predicted_fold_packets"]}
    assert packets["intrinsic_disorder_no_single_fold"]["predicted_ca_coordinates"] == []
    assert packets["intrinsic_disorder_no_single_fold"]["coordinate_model"] == "ensemble_no_single_fold_state"
    assert packets["metamorphic_fold_switching"]["predicted_ca_coordinates"] == []
    assert packets["metamorphic_fold_switching"]["coordinate_model"] == "ensemble_no_single_fold_state"

    state_rows = {
        row["canonical_family"]: row
        for row in report["family_support_rows"]
        if row["canonical_family"] in {"intrinsic_disorder_no_single_fold", "metamorphic_fold_switching"}
    }
    assert state_rows["intrinsic_disorder_no_single_fold"]["claim_type"] == "target_disorder_ensemble_claim"
    assert state_rows["metamorphic_fold_switching"]["claim_type"] == "target_metamorphic_state_claim"
    assert all(row["single_coordinate_model_allowed"] is False for row in state_rows.values())
    assert all(row["single_coordinate_fold_forced"] is False for row in state_rows.values())
