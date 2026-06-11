from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v68_builds_requested_assembly_specialist_panel() -> None:
    manifest = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V68"
        / "v68_oligomer_assembly_panel_target_manifest.json"
    )
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V68_OLIGOMER_ASSEMBLY_PANEL_200"
        / "v68_e65_oligomer_assembly_panel_200_certificate.json"
    )
    expected = {
        "V67_ASSEMBLY_REQUIRED_REPLAY": 48,
        "BIOLOGICAL_OLIGOMER_ASSEMBLY_REQUIRED": 40,
        "COILED_COIL_HELIX_BUNDLE_ASSEMBLY": 30,
        "TRUE_TRANSMEMBRANE_OLIGOMER_CHANNEL": 30,
        "SOLUBLE_MONOMERIC_HYDROPHOBIC_CORE_SENTINEL": 25,
        "COFACTOR_LIGAND_SOLUBLE_SENTINEL": 15,
        "GENERIC_COMPLEX_NOT_ASSEMBLY_CONTROL": 12,
    }
    assert manifest["batch_id"] == "V68_OLIGOMER_ASSEMBLY_PANEL_200"
    assert manifest["engine_version_used"] == "E65"
    assert manifest["baseline_engine_version"] == "E64"
    assert manifest["target_count_selected"] == 200
    assert manifest["composition_rule"] == expected
    assert Counter(row["panel_group"] for row in manifest["selected_targets"]) == expected
    assert cert["batch_id"] == "V68_OLIGOMER_ASSEMBLY_PANEL_200"
    assert cert["engine_version_used"] == "E65"
    assert cert["targets_total"] == 200
    assert cert["controls_passed"] is True


def test_v68_collapses_v67_dominant_failure_without_sentinel_regressions() -> None:
    cert = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V68"
        / "v68_oligomer_assembly_panel_certificate.json"
    )
    assert cert["status"] == "V68_E65_OLIGOMER_ASSEMBLY_PANEL_PASSED_REVIEW_REQUIRED"
    assert cert["v67_dominant_failure_replay_count"] == 48
    assert cert["v67_dominant_failure_repaired_by_e65"] == 48
    assert cert["v67_dominant_failure_remaining"] == 0
    assert cert["sentinel_regression_count"] == 0
    assert cert["failed_accepted_count"] == 0
    assert cert["next_required_batch"] == "V69_RCSB_NONREDUNDANT_200_DISCOVERY_E65"


def test_v68_reports_required_panel_metrics() -> None:
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V68_OLIGOMER_ASSEMBLY_PANEL_200"
        / "v68_e65_oligomer_assembly_panel_200_certificate.json"
    )
    assert cert["assembly_required_correct"] == 118
    assert cert["true_TM_preserved"] == 30
    assert cert["soluble_hydrophobic_not_called_assembly"] == 25
    assert cert["cofactor_not_called_assembly"] == 15
    assert cert["generic_complex_not_called_assembly"] == 12
    assert cert["coiled_coil_register_detected"] == 30
    assert cert["domain_swap_candidates_detected"] == 3
    assert cert["abstain_count"] == 12
    assert cert["accepted_count"] == 188
    assert cert["accepted_accuracy"] == pytest.approx(1.0)
    assert cert["raw_accuracy"] == pytest.approx(1.0)
    assert cert["coverage"] == pytest.approx(0.94)


def test_v68_generic_complex_controls_abstain_cleanly() -> None:
    scoring = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V68"
        / "v68_oligomer_assembly_panel_scoring_report.json"
    )
    rows = [row for row in scoring["rows"] if row["panel_group"] == "GENERIC_COMPLEX_NOT_ASSEMBLY_CONTROL"]
    assert len(rows) == 12
    assert all(row["acceptance_decision"] == "abstain_recommended" for row in rows)
    assert all(row["predicted_mechanism_class"] == "insufficient_evidence_clean_abstain" for row in rows)
    assert all(row["score_label"] == "supported" for row in rows)


def test_v68_failure_report_is_empty_but_present() -> None:
    failure_report = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V68"
        / "v68_oligomer_assembly_panel_failure_report.json"
    )
    assert failure_report["failure_cases_reported"] is True
    assert failure_report["failure_count"] == 0
    assert failure_report["failure_modes"] == {}
    assert failure_report["failure_grammar_rows"] == []
