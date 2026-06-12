from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v73_uses_expected_e68_repair_panel_composition() -> None:
    manifest = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V73"
        / "v73_e68_beta_closure_repair_target_manifest.json"
    )
    expected = {
        "V72_CLOSED_BETA_FAILURE_REPLAY": 30,
        "SOLUBLE_BETA_BARREL_SANDWICH_POSITIVE": 35,
        "MEMBRANE_BETA_BARREL_POSITIVE": 30,
        "BETA_PROPELLER_BLADE_REPEAT_POSITIVE": 30,
        "BETA_HELIX_SOLENOID_REPEAT_STACK_POSITIVE": 25,
        "ALPHA_BETA_BARREL_SENTINEL": 20,
        "PRIORITY_CONTEXT_SENTINEL": 20,
        "AMBIGUOUS_BETA_RICH_CONTROL": 10,
    }
    assert manifest["batch_id"] == "V73_BETA_CLOSURE_TOPOLOGY_REPAIR_PANEL_200"
    assert manifest["engine_version_used"] == "E68"
    assert manifest["target_count_selected"] == 200
    assert manifest["composition_rule"] == expected
    assert Counter(row["panel_group"] for row in manifest["selected_targets"]) == expected


def test_v73_repairs_closed_beta_and_preserves_priorities() -> None:
    cert = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V73"
        / "v73_e68_beta_closure_repair_certificate.json"
    )
    assert cert["status"] == "V73_E68_BETA_CLOSURE_TOPOLOGY_REPAIR_PANEL_PASSED_REVIEW_REQUIRED"
    assert cert["controls_passed"] is True
    assert cert["targets_total"] == 200
    assert cert["accepted_count"] == 190
    assert cert["accepted_supported"] == 190
    assert cert["clean_abstain_supported"] == 10
    assert cert["failed_accepted_count"] == 0
    assert cert["targeted_failed_accepted_count"] == 0
    assert cert["accepted_accuracy"] == pytest.approx(1.0)
    assert cert["coverage"] == pytest.approx(0.95)
    assert cert["v72_closed_beta_failures_repaired"] == 30
    assert cert["sentinel_regression_count"] == 0
    assert cert["ambiguous_beta_abstained"] == 10
    assert cert["failed_accepted_by_failure_mode"] == {}
    assert cert["next_required_batch"] == "V74_RCSB_NONREDUNDANT_200_DISCOVERY_E68"


def test_v73_scoring_separates_accepted_support_from_clean_abstain() -> None:
    scoring = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V73"
        / "v73_e68_beta_closure_repair_scoring_report.json"
    )
    rows = scoring["rows"]
    replay = [row for row in rows if row["panel_group"] == "V72_CLOSED_BETA_FAILURE_REPLAY"]
    ambiguous = [row for row in rows if row["panel_group"] == "AMBIGUOUS_BETA_RICH_CONTROL"]
    sentinels = [row for row in rows if row["panel_group"] == "PRIORITY_CONTEXT_SENTINEL"]

    assert len(replay) == 30
    assert all(row["predicted_mechanism_class"] == "beta_closure_topology" for row in replay)
    assert all(row["score_label"] == "supported" for row in replay)
    assert all(row["acceptance_decision"] == "abstain_recommended" for row in ambiguous)
    assert all(row["score_label"] == "supported" for row in ambiguous)
    assert Counter(row["sentinel_family"] for row in sentinels) == {
        "TRUE_TM_PRIORITY": 5,
        "ASSEMBLY_PRIORITY": 5,
        "METAL_LIGAND_PRIORITY": 5,
        "DISORDER_BOUNDARY_PRIORITY": 5,
    }
    assert all(row["score_label"] == "supported" for row in sentinels)


def test_v73_failure_report_is_empty_after_e68_repair() -> None:
    failure_report = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V73"
        / "v73_e68_beta_closure_repair_failure_report.json"
    )
    assert failure_report["failure_count"] == 0
    assert failure_report["failed_accepted_by_failure_mode"] == {}
    assert failure_report["dominant_failure_mode"] is None
