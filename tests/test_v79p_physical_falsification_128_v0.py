from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V79P_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V79P"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v79p_certificate_passes_selected_vs_enemy_falsification() -> None:
    cert = _read(V79P_ROOT / "v79p_physical_falsification_128_certificate.json")

    assert cert["status"] == "V79P_PHYSICAL_FALSIFICATION_PASSED"
    assert cert["targets_total"] == 128
    assert cert["physical_falsification_execution_run"] is True
    assert cert["runs_per_target"] == [
        "unbiased_execution",
        "selected_grammar_biased_execution",
        "wrong_grammar_biased_execution",
        "masked_grammar_biased_execution",
    ]
    assert cert["selected_grammar_beats_wrong_grammar"] == 128
    assert cert["selected_grammar_beats_masked_grammar"] == 128
    assert cert["selected_support_is_falsification_not_unbiased_only"] is True
    assert cert["native_coordinates_used_before_seal"] is False
    assert cert["native_contacts_used_before_seal"] is False
    assert cert["coordinate_or_native_leakage_blocked"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["failed_controls"] == []


def test_v79p_rows_compare_selected_wrong_masked_without_native_truth() -> None:
    rows = _read(V79P_ROOT / "v79p_physical_falsification_128_rows.json")["rows"]

    assert len(rows) == 128
    assert all(row["physical_falsification_execution_run"] for row in rows)
    assert all(row["selected_grammar_beats_wrong_grammar"] for row in rows)
    assert all(row["selected_grammar_beats_masked_grammar"] for row in rows)
    assert all(row["selected_support_is_falsification_not_unbiased_only"] for row in rows)
    assert all(row["native_coordinates_used_before_seal"] is False for row in rows)
    assert all(row["native_contacts_used_before_seal"] is False for row in rows)
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["folding_problem_solved"] is False for row in rows)
