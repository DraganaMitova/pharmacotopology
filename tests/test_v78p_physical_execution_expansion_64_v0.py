from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V78P_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V78P"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v78p_passes_physical_execution_expansion() -> None:
    cert = _read(V78P_ROOT / "v78p_physical_execution_expansion_64_certificate.json")

    assert cert["status"] == "V78P_PHYSICAL_EXECUTION_EXPANSION_PASSED"
    assert cert["targets_total"] == 64
    assert cert["category_counts"] == {
        "coiled_coil": 8,
        "repeat_solenoid": 8,
        "knotted_or_slipknot": 8,
        "signal_tm": 8,
        "secretory_disulfide": 8,
        "beta_multidomain": 8,
        "assembly_membrane_metal": 8,
        "hard_abstain_controls": 8,
    }
    assert cert["target_specific_physical_execution_run"] is True
    assert cert["unbiased_baseline_vs_grammar_biased_execution"] is True
    assert cert["accepted_rows_total"] == 56
    assert cert["accepted_rows_improved_over_unbiased"] == 56
    assert cert["hard_abstain_controls"] == 8
    assert cert["hard_abstain_controls_make_no_physical_claim"] is True
    assert cert["native_coordinates_used_before_seal"] is False
    assert cert["native_contacts_used_before_seal"] is False
    assert cert["coordinate_or_native_leakage_blocked"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["failed_controls"] == []


def test_v78p_rows_compare_unbiased_and_grammar_biased_without_native_leakage() -> None:
    rows = _read(V78P_ROOT / "v78p_physical_execution_expansion_64_rows.json")["rows"]
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]
    hard_abstain = [row for row in rows if row["physical_category"] == "hard_abstain_controls"]

    assert len(rows) == 64
    assert len(accepted) == 56
    assert len(hard_abstain) == 8
    assert all(row["target_specific_physical_execution_run"] for row in rows)
    assert all(row["native_coordinates_used_before_seal"] is False for row in rows)
    assert all(row["native_contacts_used_before_seal"] is False for row in rows)
    assert all(row["coordinate_or_native_leakage_blocked"] is True for row in rows)
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["grammar_biased_improved_over_unbiased"] for row in accepted)
    assert all(row["postseal_observable_improvement"] > 0.0 for row in accepted)
    assert all(row["physical_basis_claim_allowed"] is False for row in hard_abstain)
