from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v71_uses_fresh_four_shard_e66_discovery_targets() -> None:
    manifest = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V71"
        / "v71_rcsb_nonredundant_200_target_manifest.json"
    )
    expected = {
        "V71A_BROAD_RCSB_NONREDUNDANT": 50,
        "V71B_DISORDER_LOW_COMPLEXITY_FLEXIBLE_REGION_ENRICHED": 50,
        "V71C_BETA_BARREL_PROPELLER_REPEAT_SOLENOID_ENRICHED": 50,
        "V71D_COILED_COIL_HELIX_BUNDLE_MULTIDOMAIN_ENRICHED": 50,
    }
    assert manifest["batch_id"] == "V71_RCSB_NONREDUNDANT_200_DISCOVERY_E66"
    assert manifest["engine_version_used"] == "E66"
    assert manifest["target_count_selected"] == 200
    assert manifest["composition_rule"] == expected
    assert Counter(row["shard"] for row in manifest["selected_targets"]) == expected
    assert all(not row["fresh_exclusion_batches"] for row in manifest["selected_targets"])


def test_v71_mines_disorder_boundary_as_next_missing_word() -> None:
    cert = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V71"
        / "v71_rcsb_nonredundant_200_certificate.json"
    )
    assert cert["status"] == "V71_E66_RCSB_NONREDUNDANT_200_DISCOVERY_FAILURES_MINED_REVIEW_REQUIRED"
    assert cert["controls_passed"] is True
    assert cert["sentinel_regressions"] == 0
    assert cert["targets_total"] == 200
    assert cert["accepted_count"] == 181
    assert cert["accepted_supported"] == 111
    assert cert["clean_abstain_supported"] == 19
    assert cert["failed_accepted_count"] == 70
    assert cert["accepted_accuracy"] == pytest.approx(0.6132596685082873)
    assert cert["coverage"] == pytest.approx(0.905)
    assert cert["dominant_failure_mode"] == "disorder_misread"
    assert cert["dominant_failure_count"] == 31
    assert cert["failed_accepted_by_failure_mode"]["disorder_misread"] == 31
    assert cert["failed_accepted_by_failure_mode"]["closed_beta_topology"] == 30
    assert cert["recommended_next_engine_revision"] == "E67_DISORDER_BOUNDARY_AND_FOLD_UPON_BINDING_GRAMMAR"


def test_v71_dashboard_reports_each_shard_and_total() -> None:
    dashboard = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V71"
        / "v71_rcsb_nonredundant_200_dashboard.json"
    )
    assert set(dashboard["shards"]) == {
        "V71A_BROAD_RCSB_NONREDUNDANT",
        "V71B_DISORDER_LOW_COMPLEXITY_FLEXIBLE_REGION_ENRICHED",
        "V71C_BETA_BARREL_PROPELLER_REPEAT_SOLENOID_ENRICHED",
        "V71D_COILED_COIL_HELIX_BUNDLE_MULTIDOMAIN_ENRICHED",
        "TOTAL",
    }
    total = dashboard["shards"]["TOTAL"]
    assert total["top_failure_mode"] == "disorder_misread"
    assert total["top_missing_esperanto_word"] == "disorder_misread"
    assert total["failed_accepted_by_failure_mode"]["closed_beta_topology"] == 30
