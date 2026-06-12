from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V83P_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V83P"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v83p_certificate_passes_support_tier_selection_gate() -> None:
    cert = _read(V83P_ROOT / "v83p_real_physical_execution_selection_gate_certificate.json")

    assert cert["status"] == "V83P_REAL_PHYSICAL_EXECUTION_SELECTION_GATE_PASSED"
    assert cert["engine_version_used"] == "E77"
    assert cert["source_batch_id"] == "V83_COMPLEXITY_TOKEN_HEURISTIC_SCALING_AUDIT"
    assert cert["language_source_batch_id"] == "V82_COMPOSITIONAL_SENTENCE_PANEL_1500"
    assert cert["proxy_source_batch_id"] == "V82P_COMPOSITIONAL_PHYSICAL_HOLDOUT_GATE_768"
    assert cert["targets_total"] == 768
    assert cert["support_tier_counts"] == {
        "coarse_physical_proxy_support": 768,
        "independent_physical_holdout_support": 0,
        "language_support": 768,
    }
    assert cert["language_can_pass_without_physical_claim_count"] == 768
    assert cert["physical_proxy_can_support_language_count"] == 768
    assert cert["independent_experimental_or_coordinate_physical_proof_executed"] is False
    assert cert["physical_claim_blocked_until_independent_holdout"] is True
    assert cert["native_coordinates_used_before_seal"] is False
    assert cert["native_contacts_used_before_seal"] is False
    assert cert["coordinate_or_native_leakage_blocked"] is True
    assert cert["uses_static_observable_thresholds"] is False
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["failed_controls"] == []
    assert cert["failed_target_ids"] == []


def test_v83p_rows_keep_language_proxy_and_independent_physics_separate() -> None:
    rows = _read(V83P_ROOT / "v83p_real_physical_execution_selection_gate_rows.json")["rows"]

    assert len(rows) == 768
    assert all(row["language_support"] is True for row in rows)
    assert all(row["coarse_physical_proxy_support"] is True for row in rows)
    assert all(row["independent_physical_holdout_support"] is False for row in rows)
    assert all(row["independent_experimental_or_coordinate_physical_proof_executed"] is False for row in rows)
    assert all(row["language_can_pass_without_physical_claim"] is True for row in rows)
    assert all(row["physical_proxy_can_support_language"] is True for row in rows)
    assert all(row["physical_claim_blocked_until_independent_holdout"] is True for row in rows)
    assert all(row["native_coordinates_used_before_seal"] is False for row in rows)
    assert all(row["native_contacts_used_before_seal"] is False for row in rows)
    assert all(row["uses_static_observable_thresholds"] is False for row in rows)
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["protein_folding_solved"] is False for row in rows)
