from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V81P_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V81P"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v81p_certificate_passes_anti_tautology_holdout_gate() -> None:
    cert = _read(V81P_ROOT / "v81p_anti_tautology_physical_holdout_gate_512_certificate.json")

    assert cert["status"] == "V81P_ANTI_TAUTOLOGY_PHYSICAL_HOLDOUT_GATE_PASSED"
    assert cert["engine_version_used"] == "E75"
    assert cert["source_batch_id"] == "V81_PROTO_GRAMMAR_GENERALIZATION_PANEL"
    assert cert["targets_total"] == 512
    assert cert["runs_per_target"] == [
        "unbiased_execution",
        "selected_grammar_execution",
        "wrong_grammar_execution",
        "merge_grammar_execution",
        "masked_grammar_execution",
        "wrong_target_holdout_observable",
        "shuffled_observable_control",
    ]
    assert cert["postseal_observable_loaded_after_prediction_hash"] is True
    assert cert["observable_hash_independent_of_sealed_prediction_packet"] is True
    assert cert["selected_beats_wrong_grammar"] == 512
    assert cert["selected_beats_merge_grammar"] == 512
    assert cert["selected_beats_masked_grammar"] == 512
    assert cert["selected_beats_unbiased"] == 512
    assert cert["selected_predicts_independent_postseal_observable"] == 512
    assert cert["wrong_target_observable_control_fails"] == 512
    assert cert["shuffled_observable_control_fails"] == 512
    assert cert["anti_tautology_gate_passed"] is True
    assert cert["native_coordinates_used_before_seal"] is False
    assert cert["native_contacts_used_before_seal"] is False
    assert cert["coordinate_or_native_leakage_blocked"] is True
    assert cert["uses_static_observable_thresholds"] is False
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["failed_controls"] == []
    assert cert["failed_target_ids"] == []


def test_v81p_rows_compare_selected_wrong_merge_masked_and_wrong_observables() -> None:
    rows = _read(V81P_ROOT / "v81p_anti_tautology_physical_holdout_gate_512_rows.json")["rows"]

    assert len(rows) == 512
    assert all(row["sealed_prediction_hash"] != row["postseal_observable_hash"] for row in rows)
    assert all(row["postseal_observable_loaded_after_prediction_hash"] is True for row in rows)
    assert all(row["observable_hash_independent_of_sealed_prediction_packet"] is True for row in rows)
    assert all(row["selected_grammar_beats_wrong_grammar"] for row in rows)
    assert all(row["selected_grammar_beats_merge_grammar"] for row in rows)
    assert all(row["selected_grammar_beats_masked_grammar"] for row in rows)
    assert all(row["selected_predicts_independent_postseal_observable"] for row in rows)
    assert all(row["wrong_target_observable_control_fails"] for row in rows)
    assert all(row["shuffled_observable_control_fails"] for row in rows)
    assert all(row["anti_tautology_gate_passed"] for row in rows)
    assert all(row["native_coordinates_used_before_seal"] is False for row in rows)
    assert all(row["native_contacts_used_before_seal"] is False for row in rows)
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["protein_folding_solved"] is False for row in rows)
    assert all(
        row[execution]["uses_static_observable_thresholds"] is False
        for row in rows
        for execution in [
            "unbiased_execution",
            "selected_grammar_execution",
            "wrong_grammar_execution",
            "merge_grammar_execution",
            "masked_grammar_execution",
            "wrong_target_holdout_observable",
            "shuffled_observable_control",
        ]
    )
