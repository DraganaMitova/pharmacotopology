from __future__ import annotations

import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def v63_outputs() -> dict[str, object]:
    paths = {
        "target_manifest": ROOT / "data" / "protein_esperanto_engine" / "V63" / "v63_rcsb_500_target_manifest.json",
        "engine_declaration": ROOT / "data" / "protein_esperanto_engine" / "V63" / "v63_e61_engine_declaration.json",
        "scoring_report": ROOT / "data" / "protein_esperanto_engine" / "V63" / "v63_rcsb_500_scoring_report.json",
        "failure_report": ROOT / "data" / "protein_esperanto_engine" / "V63" / "v63_rcsb_500_failure_report.json",
        "certificate": ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V63_RCSB_500_DISCOVERY_BATCH" / "v63_rcsb_500_discovery_batch_certificate.json",
        "report": ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V63_RCSB_500_DISCOVERY_BATCH" / "V63_RCSB_500_DISCOVERY_BATCH_REPORT.md",
    }
    for path in paths.values():
        assert path.exists()
    return {
        "paths": paths,
        "certificate": json.loads(paths["certificate"].read_text(encoding="utf-8")),
        "manifest": json.loads(paths["target_manifest"].read_text(encoding="utf-8")),
        "scoring": json.loads(paths["scoring_report"].read_text(encoding="utf-8")),
        "failure_report": json.loads(paths["failure_report"].read_text(encoding="utf-8")),
    }


def test_v63_runs_500_target_e61_discovery_batch(v63_outputs: dict[str, object]) -> None:
    cert = v63_outputs["certificate"]
    manifest = v63_outputs["manifest"]
    assert isinstance(cert, dict)
    assert isinstance(manifest, dict)
    assert cert["batch_id"] == "V63_RCSB_500_DISCOVERY_BATCH"
    assert cert["batch_mode"] == "discovery"
    assert cert["engine_version_used"] == "E61"
    assert cert["targets_total"] == 500
    assert manifest["target_count_selected"] == 500
    assert manifest["target_selection_manual"] is False
    assert manifest["sequence_cluster_representative_selection"] is True
    assert cert["controls_passed"] is True
    assert cert["engine_modified_during_batch"] is False


def test_v63_reports_discovery_metrics_and_failure_modes(v63_outputs: dict[str, object]) -> None:
    cert = v63_outputs["certificate"]
    assert isinstance(cert, dict)
    assert cert["accepted_count"] == 500
    assert cert["supported_count"] == 238
    assert cert["failed_accepted_count"] == 262
    assert cert["abstain_count"] == 0
    assert cert["accepted_accuracy"] == pytest.approx(0.476)
    assert cert["raw_accuracy"] == pytest.approx(0.476)
    assert cert["coverage"] == pytest.approx(1.0)
    assert cert["failure_modes"] == {
        "disorder_misread": 15,
        "membrane_misread": 216,
        "oligomer_state_misread": 2,
        "wrong_regime": 29,
    }


def test_v63_preserves_discovery_claim_boundary(v63_outputs: dict[str, object]) -> None:
    cert = v63_outputs["certificate"]
    assert isinstance(cert, dict)
    assert cert["claim_allowed"] is False
    assert "discovery batch" in cert["claim_blocked_reason"]
    assert cert["folding_problem_solved"] is False
    assert cert["readme_check_skipped_by_user_instruction"] is True


def test_v63_failure_report_mines_missing_words(v63_outputs: dict[str, object]) -> None:
    failure_report = v63_outputs["failure_report"]
    assert isinstance(failure_report, dict)
    assert failure_report["failure_count"] == 262
    top = failure_report["missing_words_top_10"]
    assert top[0] == {"failure_type": "membrane_misread", "count": 216}
    assert top[1] == {"failure_type": "wrong_regime", "count": 29}
