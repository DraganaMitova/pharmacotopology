from __future__ import annotations

import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v66_builds_requested_adaptive_200_repair_set() -> None:
    manifest = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V66"
        / "v66_e63_fast_membrane_repair_target_manifest.json"
    )
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V66_E63_FAST_MEMBRANE_REPAIR_REPLAY_200"
        / "v66_e63_fast_membrane_repair_replay_200_certificate.json"
    )
    assert manifest["target_count_selected"] == 200
    assert manifest["composition_rule"] == {
        "E61_E62_STABLE_SENTINEL": 47,
        "V64_E62_OLIGOMER_REGRESSION": 3,
        "V64_V63_MEMBRANE_MISREAD_FAILURE": 80,
        "V65_MEMBRANE_TOPOLOGY_PANEL": 70,
    }
    assert cert["batch_id"] == "V66_E63_FAST_MEMBRANE_REPAIR_REPLAY_200"
    assert cert["batch_mode"] == "adaptive_200_repair_replay"
    assert cert["engine_version_used"] == "E63"
    assert cert["baseline_engine_version"] == "E62"
    assert cert["targets_total"] == 200
    assert cert["controls_passed"] is True


def test_v66_e63_repair_metrics_improve_over_e62_baseline() -> None:
    comparison = _read(ROOT / "data" / "protein_esperanto_engine" / "V66" / "v66_e62_vs_e63_comparison.json")
    assert comparison["baseline"]["supported_count"] == 96
    assert comparison["baseline"]["failed_accepted_count"] == 104
    assert comparison["baseline"]["abstain_count"] == 7
    assert comparison["baseline"]["accepted_accuracy"] == pytest.approx(0.46113989637305697)
    assert comparison["e63"]["supported_count"] == 130
    assert comparison["e63"]["failed_accepted_count"] == 70
    assert comparison["e63"]["abstain_count"] == 28
    assert comparison["e63"]["accepted_accuracy"] == pytest.approx(0.5930232558139535)
    assert comparison["net_changes"]["supported_delta"] == 34
    assert comparison["net_changes"]["failed_accepted_delta"] == -34
    assert comparison["net_changes"]["abstain_delta"] == 21
    assert comparison["repair_passed"] is True


def test_v66_collapses_false_membrane_without_regressing_sentinels() -> None:
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V66_E63_FAST_MEMBRANE_REPAIR_REPLAY_200"
        / "v66_e63_fast_membrane_repair_replay_200_certificate.json"
    )
    assert cert["false_membrane_repair"] == {
        "e62_false_membrane_calls": 24,
        "e62_peripheral_as_tm": 7,
        "e63_false_membrane_calls": 0,
        "e63_peripheral_as_tm": 0,
        "oligomer_regressions_remaining_under_e63": 0,
        "sentinel_regressions_under_e63": 0,
        "true_tm_missed_under_e63": 0,
    }
    assert cert["failure_modes"] == {"membrane_misread": 70}
    assert cert["status"] == "V66_E63_FAST_MEMBRANE_REPAIR_PASSED_REVIEW_REQUIRED"
    assert cert["repair_passed"] is True


def test_v66_preserves_claim_boundary_and_points_to_v67() -> None:
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V66_E63_FAST_MEMBRANE_REPAIR_REPLAY_200"
        / "v66_e63_fast_membrane_repair_replay_200_certificate.json"
    )
    assert cert["claim_allowed"] is False
    assert cert["folding_problem_solved"] is False
    assert cert["next_required_batch"] == "V67_RCSB_NONREDUNDANT_200_DISCOVERY_E63"
