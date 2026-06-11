from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v70_builds_requested_e66_repair_panel() -> None:
    manifest = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V70"
        / "v70_e66_metal_ligand_repair_target_manifest.json"
    )
    expected = {
        "V69_METAL_CLUSTER_FAILURE_REPLAY": 49,
        "V69_LIGAND_LOCKED_FAILURE_REPLAY": 14,
        "METAL_CLUSTER_POSITIVE_EXPANSION": 50,
        "LIGAND_LOCKED_POSITIVE_EXPANSION": 30,
        "TRUE_TM_METAL_CONFLICT_SENTINEL": 25,
        "ASSEMBLY_REQUIRED_METAL_CONFLICT_SENTINEL": 20,
        "V69_NON_METAL_FAILURE_TRACKING": 12,
    }
    assert manifest["batch_id"] == "V70_E66_METAL_CLUSTER_LIGAND_REPAIR_PANEL_200"
    assert manifest["engine_version_used"] == "E66"
    assert manifest["baseline_engine_version"] == "E65"
    assert manifest["target_count_selected"] == 200
    assert manifest["composition_rule"] == expected
    assert Counter(row["panel_group"] for row in manifest["selected_targets"]) == expected


def test_v70_repairs_v69_metal_and_ligand_failures_without_priority_regression() -> None:
    cert = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V70"
        / "v70_e66_metal_ligand_repair_certificate.json"
    )
    assert cert["status"] == "V70_E66_METAL_CLUSTER_LIGAND_REPAIR_PANEL_PASSED_REVIEW_REQUIRED"
    assert cert["controls_passed"] is True
    assert cert["v69_metal_cluster_failures_repaired"] == 49
    assert cert["v69_ligand_locked_failures_repaired"] == 14
    assert cert["targeted_failed_accepted_count"] == 0
    assert cert["true_TM_preserved"] == 25
    assert cert["assembly_required_preserved"] == 20
    assert cert["sentinel_regression_count"] == 0
    assert cert["non_metal_tracking_remaining"] == 7
    assert cert["next_required_batch"] == "V71_RCSB_NONREDUNDANT_200_DISCOVERY_E66"


def test_v70_reports_acceptance_split_and_tracking_failures() -> None:
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V70_E66_METAL_CLUSTER_LIGAND_REPAIR_PANEL_200"
        / "v70_e66_metal_ligand_repair_panel_200_certificate.json"
    )
    assert cert["accepted_supported"] == 193
    assert cert["clean_abstain_supported"] == 0
    assert cert["failed_accepted_count"] == 2
    assert cert["accepted_accuracy"] == pytest.approx(0.9897435897435898)
    assert cert["coverage"] == pytest.approx(0.975)
    assert cert["failure_modes"] == {"non_metal_tracking_remaining": 7}
