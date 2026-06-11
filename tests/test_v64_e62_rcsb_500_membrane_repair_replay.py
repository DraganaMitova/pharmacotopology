from __future__ import annotations

import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v64_reruns_exact_v63_500_with_e62() -> None:
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V64_E62_RCSB_500_MEMBRANE_REPAIR_REPLAY"
        / "v64_e62_rcsb_500_membrane_repair_replay_certificate.json"
    )
    manifest = _read(ROOT / "data" / "protein_esperanto_engine" / "V64" / "v64_e62_replay_target_manifest.json")
    assert cert["batch_id"] == "V64_E62_RCSB_500_MEMBRANE_REPAIR_REPLAY"
    assert cert["batch_mode"] == "paired_repair_replay"
    assert cert["engine_version_used"] == "E62"
    assert cert["baseline_engine_version"] == "E61"
    assert cert["same_targets_as_v63"] is True
    assert cert["targets_total"] == 500
    assert manifest["target_count_selected"] == 500
    assert manifest["same_targets_as_v63"] is True
    assert cert["controls_passed"] is True
    assert cert["engine_modified_during_batch"] is False


def test_v64_measures_directional_e62_membrane_repair() -> None:
    comparison = _read(ROOT / "data" / "protein_esperanto_engine" / "V64" / "v64_e61_vs_e62_comparison.json")
    assert comparison["e61"]["supported_count"] == 238
    assert comparison["e61"]["failed_accepted_count"] == 262
    assert comparison["e61"]["failure_modes"]["membrane_misread"] == 216
    assert comparison["e62"]["supported_count"] == 328
    assert comparison["e62"]["failed_accepted_count"] == 172
    assert comparison["e62"]["failure_modes"]["membrane_misread"] == 123
    assert comparison["e62"]["abstain_count"] == 0
    assert comparison["e62"]["accepted_accuracy"] == pytest.approx(0.656)
    assert comparison["net_changes"]["supported_delta"] == 90
    assert comparison["net_changes"]["failed_accepted_delta"] == -90
    assert comparison["net_changes"]["membrane_misread_delta"] == -93
    assert comparison["directional_repair"] is True
    assert comparison["membrane_still_fails"] is True


def test_v64_reports_repaired_persistent_and_new_failures() -> None:
    comparison = _read(ROOT / "data" / "protein_esperanto_engine" / "V64" / "v64_e61_vs_e62_comparison.json")
    assert comparison["repair_categories"] == {
        "new_failure": 3,
        "persistent_failure": 169,
        "repaired": 93,
        "stable_supported": 235,
    }
    assert comparison["membrane_repair"] == {
        "cofactor_or_oligomer_regression_count": 3,
        "e61_membrane_misread": 216,
        "e62_membrane_misread": 123,
        "membrane_misread_persistent": 123,
        "membrane_misread_repaired": 93,
        "soluble_false_membrane_call_count": 3,
    }
    assert comparison["new_failure_modes_from_regressions"] == {"oligomer_state_misread": 3}


def test_v64_preserves_claim_boundary() -> None:
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V64_E62_RCSB_500_MEMBRANE_REPAIR_REPLAY"
        / "v64_e62_rcsb_500_membrane_repair_replay_certificate.json"
    )
    assert cert["claim_allowed"] is False
    assert cert["directional_repair"] is True
    assert cert["membrane_still_fails"] is True
    assert cert["next_required_batch"] == "V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL"
    assert cert["folding_problem_solved"] is False
