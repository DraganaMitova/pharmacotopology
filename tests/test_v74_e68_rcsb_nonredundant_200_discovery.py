from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v74_uses_fresh_four_shard_discovery_composition() -> None:
    manifest = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V74"
        / "v74_rcsb_nonredundant_200_target_manifest.json"
    )
    expected = {
        "V74A_BROAD_RCSB_NONREDUNDANT": 50,
        "V74B_COILED_COIL_REPEAT_SOLENOID_ENRICHED": 50,
        "V74C_DISULFIDE_SECRETORY_EXTRACELLULAR_ENRICHED": 50,
        "V74D_MULTIDOMAIN_DOMAIN_SWAP_ALLOSTERY_UNUSUAL_ENRICHED": 50,
    }
    assert manifest["batch_id"] == "V74_RCSB_NONREDUNDANT_200_DISCOVERY_E68"
    assert manifest["engine_version_used"] == "E68"
    assert manifest["target_count_selected"] == 200
    assert manifest["composition_rule"] == expected
    assert Counter(row["shard"] for row in manifest["selected_targets"]) == expected
    assert all(not row["fresh_exclusion_batches"] for row in manifest["selected_targets"])


def test_v74_mines_multidomain_allostery_without_boundary_leakage() -> None:
    cert = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V74"
        / "v74_rcsb_nonredundant_200_certificate.json"
    )
    assert cert["status"] == "V74_E68_RCSB_NONREDUNDANT_200_DISCOVERY_FAILURES_MINED_REVIEW_REQUIRED"
    assert cert["targets_total"] == 200
    assert cert["accepted_count"] == 107
    assert cert["accepted_supported"] == 11
    assert cert["clean_abstain_supported"] == 93
    assert cert["failed_accepted_count"] == 96
    assert cert["accepted_accuracy"] == pytest.approx(0.102803738317757)
    assert cert["coverage"] == pytest.approx(0.535)
    assert cert["controls_passed"] is True
    assert cert["sentinel_regressions"] == 0
    assert cert["withheld_context_leakage_detected"] is False
    assert cert["dominant_failure_mode"] == "multidomain_allostery"
    assert cert["dominant_failure_count"] == 33
    assert cert["top_missing_esperanto_word"] == "multidomain_allostery"
    assert cert["recommended_next_engine_revision"] == "E69_MULTIDOMAIN_ALLOSTERIC_ARCHITECTURE_GRAMMAR"
    assert cert["next_required_batch"] == "E69_AND_V75_REPAIR_PANEL"
    assert cert["claim_allowed"] is False


def test_v74_failure_taxonomy_is_concentrated() -> None:
    failure_report = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V74"
        / "v74_rcsb_nonredundant_200_failure_report.json"
    )
    assert failure_report["failure_count"] == 96
    assert failure_report["failed_accepted_by_failure_mode"] == {
        "assembly_required_missed": 0,
        "closed_beta_topology": 2,
        "coiled_coil_register": 17,
        "disorder_misread": 0,
        "disulfide_secretory_redox_context": 28,
        "domain_swapping": 1,
        "knotted_topology": 0,
        "membrane_topology_missed_or_misread": 0,
        "metal_ligand_basin": 0,
        "multidomain_allostery": 33,
        "other": 0,
        "repeat_solenoid_topology": 15,
        "signal_peptide_vs_true_TM": 0,
    }
    top = failure_report["failure_grammar_rows"][0]["autopsy_sentence"]
    assert "The engine thought:" in top
    assert "Reality showed:" in top
    assert "Missing Esperanto word:" in top


def test_v74_dashboard_separates_accepted_support_and_clean_abstain() -> None:
    dashboard = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V74"
        / "v74_rcsb_nonredundant_200_dashboard.json"
    )["shards"]
    total = dashboard["TOTAL"]
    assert total["targets_total"] == 200
    assert total["accepted_count"] == 107
    assert total["accepted_supported"] == 11
    assert total["failed_accepted"] == 96
    assert total["clean_abstain_supported"] == 93
    assert total["top_failure_mode"] == "multidomain_allostery"

    assert dashboard["V74B_COILED_COIL_REPEAT_SOLENOID_ENRICHED"]["top_missing_esperanto_word"] == "coiled_coil_register"
    assert dashboard["V74C_DISULFIDE_SECRETORY_EXTRACELLULAR_ENRICHED"]["top_missing_esperanto_word"] == "disulfide_secretory_redox_context"
    assert dashboard["V74D_MULTIDOMAIN_DOMAIN_SWAP_ALLOSTERY_UNUSUAL_ENRICHED"]["top_missing_esperanto_word"] == "multidomain_allostery"
