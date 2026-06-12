from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V80P_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V80P"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v80p_certificate_passes_independent_postseal_holdout_gate() -> None:
    cert = _read(V80P_ROOT / "v80p_independent_physical_holdout_gate_256_certificate.json")

    assert cert["status"] == "V80P_INDEPENDENT_PHYSICAL_HOLDOUT_GATE_PASSED"
    assert cert["engine_version_used"] == "E74"
    assert cert["source_batch_id"] == "V80_CANDIDATE_WORD_LEXICON_DELTA_TRIAGE_77"
    assert cert["targets_total"] == 256
    assert cert["runs_per_target"] == [
        "unbiased_execution",
        "selected_grammar_execution",
        "wrong_grammar_execution",
        "merge_grammar_execution",
        "masked_grammar_execution",
    ]
    assert cert["selected_beats_wrong_grammar"] == 256
    assert cert["selected_beats_merge_grammar"] == 256
    assert cert["selected_beats_masked_grammar"] == 256
    assert cert["selected_beats_unbiased"] == 256
    assert cert["selected_predicts_independent_postseal_observable"] == 256
    assert cert["independent_postseal_observable_gate_passed"] is True
    assert cert["native_coordinates_used_before_seal"] is False
    assert cert["native_contacts_used_before_seal"] is False
    assert cert["coordinate_or_native_leakage_blocked"] is True
    assert cert["uses_static_observable_thresholds"] is False
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["failed_controls"] == []
    assert cert["failed_target_ids"] == []


def test_v80p_rows_run_selected_wrong_merge_masked_without_native_truth() -> None:
    rows = _read(V80P_ROOT / "v80p_independent_physical_holdout_gate_256_rows.json")["rows"]

    assert len(rows) == 256
    assert all(
        row["runs_per_target"]
        == [
            "unbiased_execution",
            "selected_grammar_execution",
            "wrong_grammar_execution",
            "merge_grammar_execution",
            "masked_grammar_execution",
        ]
        for row in rows
    )
    assert all(row["selected_grammar_beats_wrong_grammar"] for row in rows)
    assert all(row["selected_grammar_beats_merge_grammar"] for row in rows)
    assert all(row["selected_grammar_beats_masked_grammar"] for row in rows)
    assert all(row["selected_grammar_beats_unbiased"] for row in rows)
    assert all(row["selected_predicts_independent_postseal_observable"] for row in rows)
    assert all(row["independent_postseal_observable_gate_passed"] for row in rows)
    assert all(row["native_coordinates_used_before_seal"] is False for row in rows)
    assert all(row["native_contacts_used_before_seal"] is False for row in rows)
    assert all(row["coordinate_or_native_leakage_blocked"] is True for row in rows)
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["folding_problem_solved"] is False for row in rows)
    assert all(
        row[execution]["uses_static_observable_thresholds"] is False
        for row in rows
        for execution in [
            "unbiased_execution",
            "selected_grammar_execution",
            "wrong_grammar_execution",
            "merge_grammar_execution",
            "masked_grammar_execution",
        ]
    )
