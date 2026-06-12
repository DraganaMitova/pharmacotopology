from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V76P_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V76P"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v76p_runs_target_specific_physical_execution_without_physical_claim() -> None:
    cert = _read(V76P_ROOT / "v76p_physical_execution_pilot_certificate.json")

    assert cert["status"] == "V76P_PHYSICAL_EXECUTION_PILOT_PASSED"
    assert cert["batch_id"] == "V76P_OPENMM_COARSE_TARGET_EXECUTION_PILOT"
    assert cert["engine_version_used"] == "E70"
    assert cert["targets_total"] == 8
    assert cert["accepted_learned_grammar_targets"] == 4
    assert cert["hard_abstain_controls"] == 4
    assert cert["target_specific_physical_execution_run"] is True
    assert cert["native_coordinates_used_before_seal"] is False
    assert cert["native_contacts_used_before_seal"] is False
    assert cert["accepted_target_improvements"] == 4
    assert cert["hard_abstain_controls_executed"] == 4
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["folding_problem_solved"] is False


def test_v76p_rows_compare_unbiased_and_grammar_biased_execution() -> None:
    cert = _read(V76P_ROOT / "v76p_physical_execution_pilot_certificate.json")
    rows = cert["rows"]
    accepted = [row for row in rows if row["pilot_role"] == "accepted_learned_grammar"]
    controls = [row for row in rows if row["pilot_role"] == "hard_abstain_control"]

    assert len(accepted) == 4
    assert len(controls) == 4
    assert all(row["target_specific_physical_execution_run"] is True for row in rows)
    assert all(row["native_coordinates_used_before_seal"] is False for row in rows)
    assert all(row["native_contacts_used_before_seal"] is False for row in rows)
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["folding_problem_solved"] is False for row in rows)
    assert all(row["unbiased_baseline"]["backend"] == row["execution_backend"] for row in rows)
    assert all(row["grammar_biased_execution"]["backend"] == row["execution_backend"] for row in rows)
    assert all(row["execution_backend"] in {"openmm_reference", "deterministic_equivalent"} for row in rows)
    assert all(row["grammar_biased_improved_over_unbiased"] for row in accepted)
    assert all(row["postseal_observable_improvement"] > 0.0 for row in accepted)
    assert all(row["acceptance_decision"] == "accepted" for row in accepted)
    assert all(row["acceptance_decision"] == "abstain_recommended" for row in controls)


def test_v76p_rows_file_matches_certificate_rows() -> None:
    cert = _read(V76P_ROOT / "v76p_physical_execution_pilot_certificate.json")
    rows_packet = _read(V76P_ROOT / "v76p_physical_execution_rows.json")

    assert rows_packet["kind"] == "V76P_PHYSICAL_EXECUTION_ROWS_v0"
    assert rows_packet["rows"] == cert["rows"]
