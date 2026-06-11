from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v72_uses_expected_e67_repair_panel_composition() -> None:
    manifest = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V72"
        / "v72_e67_disorder_boundary_repair_target_manifest.json"
    )
    expected = {
        "V71_DISORDER_FAILURE_REPLAY": 31,
        "DISORDER_BOUNDARY_POSITIVE_EXPANSION": 60,
        "GENERIC_OLIGOMER_CONTROL": 25,
        "TRUE_TM_DISORDER_CONFLICT_SENTINEL": 20,
        "ASSEMBLY_REQUIRED_DISORDER_PRIORITY_SENTINEL": 14,
        "METAL_LIGAND_DISORDER_CONFLICT_SENTINEL": 20,
        "V71_CLOSED_BETA_TRACKING": 30,
    }
    assert manifest["batch_id"] == "V72_E67_DISORDER_BOUNDARY_REPAIR_PANEL_200"
    assert manifest["engine_version_used"] == "E67"
    assert manifest["target_count_selected"] == 200
    assert manifest["composition_rule"] == expected
    assert Counter(row["panel_group"] for row in manifest["selected_targets"]) == expected


def test_v72_repairs_disorder_failures_and_preserves_sentinels() -> None:
    cert = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V72"
        / "v72_e67_disorder_boundary_repair_certificate.json"
    )
    assert cert["status"] == "V72_E67_DISORDER_BOUNDARY_REPAIR_PANEL_PASSED_REVIEW_REQUIRED"
    assert cert["controls_passed"] is True
    assert cert["targets_total"] == 200
    assert cert["accepted_count"] == 200
    assert cert["accepted_supported"] == 170
    assert cert["failed_accepted_count"] == 30
    assert cert["targeted_failed_accepted_count"] == 0
    assert cert["accepted_accuracy"] == pytest.approx(0.85)
    assert cert["coverage"] == pytest.approx(1.0)
    assert cert["v71_disorder_failures_repaired"] == 31
    assert cert["disorder_expansion_supported"] == 60
    assert cert["generic_oligomer_controls_preserved"] == 25
    assert cert["true_TM_preserved"] == 20
    assert cert["assembly_required_preserved"] == 14
    assert cert["metal_ligand_preserved"] == 20
    assert cert["sentinel_regression_count"] == 0
    assert cert["closed_beta_tracking_remaining"] == 30
    assert cert["failed_accepted_by_failure_mode"] == {"closed_beta_tracking_remaining": 30}
    assert cert["next_required_batch"] == "V73_BETA_CLOSURE_TOPOLOGY_REPAIR_PANEL_200"


def test_v72_scoring_keeps_beta_tracking_separate_from_targeted_failures() -> None:
    scoring = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V72"
        / "v72_e67_disorder_boundary_repair_scoring_report.json"
    )
    rows = scoring["rows"]
    assert sum(row["score_label"] == "supported" for row in rows if row["panel_group"] == "V71_DISORDER_FAILURE_REPLAY") == 31
    beta_rows = [row for row in rows if row["panel_group"] == "V71_CLOSED_BETA_TRACKING"]
    assert len(beta_rows) == 30
    assert all(row["tracking_only"] is True for row in beta_rows)
    assert all(row["score_label"] == "contradicted" for row in beta_rows)
    assert all(row["predicted_mechanism_class"] == "oligomerization_controlled_folding" for row in beta_rows)
