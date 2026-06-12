from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V82P_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V82P"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v82p_certificate_passes_compositional_physical_holdout_gate() -> None:
    cert = _read(V82P_ROOT / "v82p_compositional_physical_holdout_gate_768_certificate.json")

    assert cert["status"] == "V82P_COMPOSITIONAL_PHYSICAL_HOLDOUT_GATE_PASSED"
    assert cert["engine_version_used"] == "E76"
    assert cert["source_batch_id"] == "V82_COMPOSITIONAL_SENTENCE_PANEL_1500"
    assert cert["targets_total"] == 768
    assert cert["runs_per_target"] == [
        "unbiased_execution",
        "selected_sentence_execution",
        "bag_of_words_execution",
        "wrong_order_execution",
        "wrong_head_execution",
        "masked_clause_execution",
        "wrong_target_observable_control",
    ]
    assert cert["postseal_observable_loaded_after_prediction_hash"] is True
    assert cert["observable_hash_independent_of_sealed_prediction_packet"] is True
    assert cert["selected_sentence_beats_bag_of_words"] == 768
    assert cert["selected_sentence_beats_wrong_order"] == 768
    assert cert["selected_sentence_beats_wrong_head"] == 768
    assert cert["selected_sentence_beats_masked_clause"] == 768
    assert cert["selected_sentence_beats_unbiased"] == 768
    assert cert["selected_sentence_predicts_independent_postseal_observable"] == 768
    assert cert["wrong_target_observable_control_fails"] == 768
    assert cert["compositional_physical_holdout_gate_passed"] is True
    assert cert["native_coordinates_used_before_seal"] is False
    assert cert["native_contacts_used_before_seal"] is False
    assert cert["coordinate_or_native_leakage_blocked"] is True
    assert cert["uses_static_observable_thresholds"] is False
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["failed_controls"] == []
    assert cert["failed_target_ids"] == []


def test_v82p_rows_compare_selected_sentence_against_composition_controls() -> None:
    rows = _read(V82P_ROOT / "v82p_compositional_physical_holdout_gate_768_rows.json")["rows"]

    assert len(rows) == 768
    assert all(row["sealed_prediction_hash"] != row["postseal_observable_hash"] for row in rows)
    assert all(row["postseal_observable_loaded_after_prediction_hash"] is True for row in rows)
    assert all(row["observable_hash_independent_of_sealed_prediction_packet"] is True for row in rows)
    assert all(row["selected_sentence_beats_bag_of_words"] for row in rows)
    assert all(row["selected_sentence_beats_wrong_order"] for row in rows)
    assert all(row["selected_sentence_beats_wrong_head"] for row in rows)
    assert all(row["selected_sentence_beats_masked_clause"] for row in rows)
    assert all(row["selected_sentence_predicts_independent_postseal_observable"] for row in rows)
    assert all(row["wrong_target_observable_control_fails"] for row in rows)
    assert all(row["compositional_physical_holdout_gate_passed"] for row in rows)
    assert all(row["native_coordinates_used_before_seal"] is False for row in rows)
    assert all(row["native_contacts_used_before_seal"] is False for row in rows)
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["protein_folding_solved"] is False for row in rows)
    assert all(
        row[execution]["uses_static_observable_thresholds"] is False
        for row in rows
        for execution in [
            "unbiased_execution",
            "selected_sentence_execution",
            "bag_of_words_execution",
            "wrong_order_execution",
            "wrong_head_execution",
            "masked_clause_execution",
            "wrong_target_observable_control",
        ]
    )
