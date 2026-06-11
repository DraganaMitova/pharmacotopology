from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v69_uses_fresh_four_shard_e65_discovery_targets() -> None:
    manifest = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V69"
        / "v69_rcsb_nonredundant_200_target_manifest.json"
    )
    expected = {
        "V69A_BROAD_RCSB_NONREDUNDANT": 50,
        "V69B_COFACTOR_LIGAND_METAL_ENRICHED": 50,
        "V69C_ASSEMBLY_COMPLEX_OLIGOMER_ENRICHED": 50,
        "V69D_HARD_TOPOLOGY_ENRICHED": 50,
    }
    assert manifest["batch_id"] == "V69_RCSB_NONREDUNDANT_200_DISCOVERY_E65"
    assert manifest["engine_version_used"] == "E65"
    assert manifest["target_count_selected"] == 200
    assert manifest["composition_rule"] == expected
    assert Counter(row["shard"] for row in manifest["selected_targets"]) == expected
    assert all(not row["fresh_exclusion_batches"] for row in manifest["selected_targets"])


def test_v69_mines_metal_cluster_geometry_as_next_missing_word() -> None:
    cert = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V69"
        / "v69_rcsb_nonredundant_200_certificate.json"
    )
    assert cert["status"] == "V69_E65_RCSB_NONREDUNDANT_200_DISCOVERY_FAILURES_MINED_REVIEW_REQUIRED"
    assert cert["controls_passed"] is True
    assert cert["sentinel_regressions"] == 0
    assert cert["targets_total"] == 200
    assert cert["accepted_count"] == 173
    assert cert["accepted_supported"] == 81
    assert cert["clean_abstain_supported"] == 27
    assert cert["failed_accepted_count"] == 92
    assert cert["accepted_accuracy"] == pytest.approx(0.4682080924855491)
    assert cert["coverage"] == pytest.approx(0.865)
    assert cert["top_missing_esperanto_word"] == "metal_cluster_geometry"
    assert cert["dominant_failure_count"] == 49
    assert cert["recommended_next_engine_revision"] == "E66_METAL_CLUSTER_AND_LIGAND_LOCKED_BASIN_GRAMMAR"
    assert cert["failed_accepted_by_failure_mode"]["metal_cluster_geometry"] == 49
    assert cert["failed_accepted_by_failure_mode"]["ligand_locked_basin"] == 14


def test_v69_dashboard_reports_each_shard_and_total() -> None:
    dashboard = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V69"
        / "v69_rcsb_nonredundant_200_dashboard.json"
    )
    assert set(dashboard["shards"]) == {
        "V69A_BROAD_RCSB_NONREDUNDANT",
        "V69B_COFACTOR_LIGAND_METAL_ENRICHED",
        "V69C_ASSEMBLY_COMPLEX_OLIGOMER_ENRICHED",
        "V69D_HARD_TOPOLOGY_ENRICHED",
        "TOTAL",
    }
    total = dashboard["shards"]["TOTAL"]
    assert total["top_missing_esperanto_word"] == "metal_cluster_geometry"
    assert total["failed_accepted_by_failure_mode"]["metal_cluster_geometry"] == 49
    assert dashboard["shards"]["V69B_COFACTOR_LIGAND_METAL_ENRICHED"]["top_missing_esperanto_word"] == "metal_cluster_geometry"
