from __future__ import annotations

import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v65_builds_five_way_membrane_topology_panel() -> None:
    manifest = _read(ROOT / "data" / "protein_esperanto_engine" / "V65" / "v65_membrane_topology_panel_manifest.json")
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL"
        / "v65_membrane_topology_esperanto_panel_certificate.json"
    )
    assert manifest["panel_target_count"] == 70
    assert manifest["targets_per_group"] == 14
    assert manifest["panel_groups"] == {
        "A_TRUE_TRANSMEMBRANE_TOPOLOGY": 14,
        "B_SOLUBLE_HYDROPHOBIC_NO_TOPOLOGY": 14,
        "C_COFACTOR_BOUND_SOLUBLE_HYDROPHOBIC_POCKET": 14,
        "D_OLIGOMERIC_SOLUBLE_INTERFACE_HYDROPHOBICITY": 14,
        "E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE": 14,
    }
    assert cert["engine_version_used"] == "E62"
    assert cert["controls_passed"] is True


def test_v65_exposes_e62_topology_false_membrane_failures() -> None:
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL"
        / "v65_membrane_topology_esperanto_panel_certificate.json"
    )
    assert cert["status"] == "V65_MEMBRANE_TOPOLOGY_PANEL_FAILURES_REVIEW_REQUIRED"
    assert cert["targets_total"] == 70
    assert cert["supported_count"] == 49
    assert cert["failed_accepted_count"] == 21
    assert cert["abstain_count"] == 7
    assert cert["accepted_accuracy"] == pytest.approx(2 / 3)
    assert cert["raw_accuracy"] == pytest.approx(0.7)
    assert cert["false_membrane_call_count"] == 21
    assert cert["peripheral_false_transmembrane_count"] == 7
    assert cert["true_transmembrane_missed_count"] == 0
    assert cert["failure_modes"] == {
        "peripheral_misread_as_transmembrane": 7,
        "soluble_hydrophobic_false_membrane": 14,
    }
    assert cert["engine_revision_required"] is True
    assert cert["engine_revision_recommended"] == "E63_MEMBRANE_TOPOLOGY_GRAMMAR_REVISION"


def test_v65_distinctions_are_recorded_before_claims() -> None:
    manifest = _read(ROOT / "data" / "protein_esperanto_engine" / "V65" / "v65_membrane_topology_panel_manifest.json")
    distinctions = set(manifest["required_esperanto_distinctions"])
    assert "transmembrane helix" in distinctions
    assert "transmembrane beta barrel" in distinctions
    assert "amphipathic peripheral helix" in distinctions
    assert "soluble hydrophobic core" in distinctions
    assert "cofactor-buried hydrophobic pocket" in distinctions
    assert "oligomeric interface hydrophobicity" in distinctions
    assert "inside/outside topology" in distinctions


def test_v65_preserves_claim_boundary() -> None:
    cert = _read(
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL"
        / "v65_membrane_topology_esperanto_panel_certificate.json"
    )
    assert cert["claim_allowed"] is False
    assert cert["folding_problem_solved"] is False
    assert "grammar-mining gate" in cert["claim_blocked_reason"]
