from __future__ import annotations

import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def v62_outputs() -> dict[str, object]:
    paths = {
        "target_manifest": ROOT / "data" / "protein_esperanto_engine" / "V62" / "v62_e61_repair_target_manifest.json",
        "engine_declaration": ROOT / "data" / "protein_esperanto_engine" / "V62" / "v62_e61_engine_declaration.json",
        "scoring_report": ROOT / "data" / "protein_esperanto_engine" / "V62" / "v62_e61_repair_scoring_report.json",
        "failure_report": ROOT / "data" / "protein_esperanto_engine" / "V62" / "v62_e61_repair_failure_report.json",
        "comparison": ROOT / "data" / "protein_esperanto_engine" / "V62" / "v62_e60_vs_e61_comparison.json",
        "certificate": ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V62_E61_REPAIR_AND_SATURATION_RERUN" / "v62_e61_repair_and_saturation_rerun_certificate.json",
        "report": ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V62_E61_REPAIR_AND_SATURATION_RERUN" / "V62_E61_REPAIR_AND_SATURATION_RERUN_REPORT.md",
        "claim_row": ROOT / "data" / "protein_esperanto_engine" / "V62" / "v62_campaign_claim_row.json",
        "claim_ledger": ROOT / "data" / "protein_esperanto_engine" / "V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN" / "ledgers" / "claim_ledger_v0.json",
    }
    for path in paths.values():
        assert path.exists()
    return {
        "paths": paths,
        "certificate": json.loads(paths["certificate"].read_text(encoding="utf-8")),
        "comparison": json.loads(paths["comparison"].read_text(encoding="utf-8")),
        "manifest": json.loads(paths["target_manifest"].read_text(encoding="utf-8")),
        "scoring": json.loads(paths["scoring_report"].read_text(encoding="utf-8")),
    }


def test_v62_reruns_same_v61_targets_with_e61(v62_outputs: dict[str, object]) -> None:
    cert = v62_outputs["certificate"]
    manifest = v62_outputs["manifest"]
    assert isinstance(cert, dict)
    assert isinstance(manifest, dict)
    assert cert["batch_id"] == "V62_E61_REPAIR_AND_SATURATION_RERUN"
    assert cert["batch_mode"] == "repair"
    assert cert["engine_version_used"] == "E61"
    assert cert["baseline_engine_version"] == "E60"
    assert cert["same_targets_as_v61"] is True
    assert cert["known_failure_repair_probe"] is True
    assert manifest["target_count_selected"] == 100
    assert manifest["same_targets_as_v61"] is True
    assert cert["engine_modified_during_batch"] is False
    assert cert["readme_modified"] is False
    assert cert["controls_passed"] is True


def test_v62_compares_against_frozen_e60_baseline(v62_outputs: dict[str, object]) -> None:
    comparison = v62_outputs["comparison"]
    assert isinstance(comparison, dict)
    assert comparison["e60"]["accepted_count"] == 81
    assert comparison["e60"]["supported_count"] == 8
    assert comparison["e60"]["failed_accepted_count"] == 73
    assert comparison["e60"]["accepted_accuracy"] == pytest.approx(0.09876543209876543)
    assert comparison["e61"]["accepted_count"] == 96
    assert comparison["e61"]["supported_count"] == 90
    assert comparison["e61"]["failed_accepted_count"] == 6
    assert comparison["e61"]["accepted_accuracy"] == pytest.approx(0.9375)
    assert comparison["net_changes"]["failed_accepted_delta"] == -67
    assert comparison["net_changes"]["supported_delta"] == 82
    assert comparison["net_changes"]["coverage_delta"] == pytest.approx(0.15)


def test_v62_repairs_v61_missing_context_modes(v62_outputs: dict[str, object]) -> None:
    comparison = v62_outputs["comparison"]
    assert isinstance(comparison, dict)
    modes = comparison["failure_mode_distribution_change"]
    assert modes["cofactor_ligand_missing"] == {"e60": 43, "e61": 0, "delta": -43}
    assert modes["oligomer_state_misread"] == {"e60": 13, "e61": 0, "delta": -13}
    assert modes["membrane_misread"] == {"e60": 11, "e61": 0, "delta": -11}
    assert comparison["repair_categories"]["repaired"] == 82
    assert comparison["repair_categories"]["stable_supported"] == 8
    assert comparison["repair_categories"]["unchanged_failure"] == 10


def test_v62_context_marks_are_explicit_and_leakage_blocked(v62_outputs: dict[str, object]) -> None:
    cert = v62_outputs["certificate"]
    assert isinstance(cert, dict)
    controls = {row["control_id"]: row for row in cert["controls"]}
    assert controls["context_marks_present_for_v61_missing_context_failures"]["passed"] is True
    assert controls["coordinate_leakage_control"]["passed"] is True
    assert controls["alphafold_leakage_control"]["passed"] is True
    assert controls["holdout_opened_before_seal_control"]["passed"] is True
    assert controls["internal_runtime_leakage_control"]["passed"] is True
    assert controls["readme_modified_false"]["passed"] is True


def test_v62_preserves_claim_boundary(v62_outputs: dict[str, object]) -> None:
    cert = v62_outputs["certificate"]
    assert isinstance(cert, dict)
    assert cert["directional_repair"] is True
    assert cert["claim_allowed"] is False
    assert "same-target repair probe" in cert["claim_blocked_reason"]
    assert cert["folding_problem_solved"] is False
