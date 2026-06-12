from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v75_uses_expected_self_decision_composition() -> None:
    manifest = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V75"
        / "v75_self_deciding_acceptance_cortex_target_manifest.json"
    )
    expected = {
        "V74_MULTIDOMAIN_ALLOSTERY_FAILURE_REPLAY": 33,
        "MULTIDOMAIN_HINGE_DOMAIN_SWAP_ALLOSTERY_POSITIVE": 40,
        "MODULAR_ARCHITECTURE_INTERDOMAIN_LOCK_POSITIVE": 25,
        "MONOMERIC_GLOBULAR_SENTINEL": 20,
        "ASSEMBLY_REQUIRED_SENTINEL": 20,
        "MEMBRANE_TM_SENTINEL": 20,
        "METAL_LIGAND_SENTINEL": 15,
        "DISORDER_BOUNDARY_SENTINEL": 15,
        "MISSING_WORD_CANDIDATE_CLEAN_ABSTAIN_CONTROL": 12,
    }
    assert manifest["batch_id"] == "V75_SELF_DECIDING_ACCEPTANCE_CORTEX"
    assert manifest["engine_version_used"] == "E69"
    assert manifest["target_count_selected"] == 200
    assert manifest["composition_rule"] == expected
    assert Counter(row["panel_group"] for row in manifest["selected_targets"]) == expected


def test_v75_passes_self_deciding_zero_failed_accepted_invariant() -> None:
    cert = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V75"
        / "v75_self_deciding_acceptance_cortex_certificate.json"
    )
    assert cert["status"] == "BATCH_PASSED_SELF_DECIDING_ZERO_FAILED_ACCEPTED"
    assert cert["zero_failed_accepted_required"] is True
    assert cert["self_decision_cortex_enabled"] is True
    assert cert["targets_total"] == 200
    assert cert["accepted_count"] == 187
    assert cert["accepted_supported"] == 187
    assert cert["clean_abstain_supported"] == 12
    assert cert["blocked_supported"] == 0
    assert cert["failed_accepted_count"] == 0
    assert cert["accepted_accuracy"] == pytest.approx(1.0)
    assert cert["coverage"] == pytest.approx(0.935)
    assert cert["v74_multidomain_failures_repaired"] == 33
    assert cert["sentinel_regressions"] == 0
    assert cert["missing_candidate_controls_cleanly_abstained"] == 12
    assert cert["accepted_operator_basis_stable"] == 187
    assert cert["accepted_cross_view_bound"] == 187
    assert cert["accepted_temporal_bound"] == 187
    assert cert["dominance_law_for_acceptance"] == "single_dominant_learned_mechanism_bound_across_views"
    assert cert["coefficient_probe_mode"] == "endogenous_observed_operator_permutations_no_static_scale_range"
    assert cert["candidate_grammars_implemented"] is False
    assert cert["candidate_grammar_acceptance_role"] == "clean_abstain_until_revision_implements_grammar"
    assert cert["unknown_word_abstentions_by_word"] == {
        "coiled_coil_register": 4,
        "disulfide_secretory_redox_context": 5,
        "repeat_solenoid_topology": 3,
    }
    assert cert["top_missing_esperanto_word"] == "disulfide_secretory_redox_context"
    assert cert["next_recommended_engine_revision"] == "E70_SECRETORY_DISULFIDE_REDOX_TOPOLOGY_GRAMMAR"
    assert cert["next_required_batch"] == "V76_SECRETORY_DISULFIDE_REPAIR_PANEL_200"
    assert cert["real_physical_calibration_inputs_used"] is True
    assert cert["real_physical_calibration_row_count"] == 8
    assert cert["real_physical_calibration_source_coordinate_database"] == "RCSB_PDB"
    assert cert["real_physical_calibration_target_native_excluded"] is True
    assert cert["real_physical_calibration_target_native_contacts_used_before_prediction"] is False
    assert cert["real_physical_calibration_coordinate_truth_used_as_prediction_input"] is False
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["physical_grounding_statuses"] == {
        "real_physical_calibration_inputs_loaded_openmm_available_but_target_physical_execution_not_run": 200
    }


def test_v75_missing_candidate_controls_abstain_by_self_decision() -> None:
    scoring = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V75"
        / "v75_self_deciding_acceptance_cortex_scoring_report.json"
    )
    rows = scoring["rows"]
    missing = [row for row in rows if row["panel_group"] == "MISSING_WORD_CANDIDATE_CLEAN_ABSTAIN_CONTROL"]
    assert len(missing) == 12
    assert all(row["acceptance_decision"] == "abstain_recommended" for row in missing)
    assert all(row["final_self_decision"] == "clean_abstain_missing_word" for row in missing)
    assert all(row["clean_abstain_supported"] for row in missing)
    assert all(row["score_label"] == "supported" for row in missing)
    assert all(row["missing_word_candidate"] in {
        "disulfide_secretory_redox_context",
        "coiled_coil_register",
        "repeat_solenoid_topology",
    } for row in missing)
    assert all(row["internal_consensus"] == "candidate_missing_word_competes_with_learned_grammar" for row in missing)


def test_v75_accepted_predictions_are_bound_and_operator_basis_stable() -> None:
    scoring = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V75"
        / "v75_self_deciding_acceptance_cortex_scoring_report.json"
    )
    accepted = [row for row in scoring["rows"] if row["acceptance_decision"] == "accepted"]
    assert len(accepted) == 187
    assert all(row["dominance_law"] == "single_dominant_learned_mechanism_bound_across_views" for row in accepted)
    assert all(row["operator_basis_stability"] == "stable_under_endogenous_operator_basis_probe" for row in accepted)
    assert all(row["coefficient_probe_mode"] == "endogenous_observed_operator_permutations_no_static_scale_range" for row in accepted)
    assert all(row["cross_view_missing_view_count"] == 0 for row in accepted)
    assert all(row["temporal_binding"] == "selected_observables_temporally_coherent" for row in accepted)
    assert all(row["real_physical_calibration_inputs_used"] is True for row in accepted)
    assert all(row["physical_basis_claim_allowed"] is False for row in accepted)
    assert all(
        row["physical_grounding_status"]
        == "real_physical_calibration_inputs_loaded_openmm_available_but_target_physical_execution_not_run"
        for row in accepted
    )


def test_v75_conservative_abstain_when_operator_basis_is_sensitive() -> None:
    scoring = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V75"
        / "v75_self_deciding_acceptance_cortex_scoring_report.json"
    )
    sensitive = [
        row
        for row in scoring["rows"]
        if row["operator_basis_stability"] == "coefficient_assignment_sensitive"
    ]
    assert [row["target_id"] for row in sensitive] == [
        "V75_049_MULTIDOMAIN_HINGE_DOMAIN_SWAP_ALLOSTERY_POSITIVE_1PLF_1"
    ]
    assert sensitive[0]["acceptance_decision"] == "abstain_recommended"
    assert sensitive[0]["final_self_decision"] == "clean_abstain_conflict"
    assert sensitive[0]["score_label"] == "abstained"
    assert sensitive[0]["accepted_supported"] is False


def test_v75_real_physical_calibration_inputs_are_built_from_locked_coordinate_data() -> None:
    calibration = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V75"
        / "physical_calibration"
        / "v75_real_physical_calibration_inputs.json"
    )
    assert calibration["kind"] == "PROTEIN_ESPERANTO_REAL_PHYSICAL_CALIBRATION_INPUTS_v0"
    assert calibration["source_coordinate_database"] == "RCSB_PDB"
    assert calibration["calibration_input_type"] == "coordinate_derived_native_contact_summary"
    assert calibration["row_count"] == 8
    assert calibration["target_native_excluded_from_calibration"] is True
    assert calibration["target_native_contacts_used_before_prediction"] is False
    assert calibration["coordinate_truth_used_as_prediction_input"] is False
    assert calibration["leave_one_target_out_calibration"] is True
    assert calibration["universal_physical_law_claim_allowed"] is False
    assert calibration["folding_problem_solved"] is False
    assert {row["source_accession"] for row in calibration["calibration_rows"]} == {
        "1MBN:A",
        "2LZM:A",
        "1TEN:A",
        "1CSP:A",
        "1TIM:A",
        "1PGA:A",
        "4AKE:A",
        "1CLL:A",
    }
    assert "native_contact_count" in calibration["observable_families"]
    assert calibration["calibration_hash"]


def test_v75_physical_gate_uses_real_calibration_but_blocks_physics_claims() -> None:
    cert = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V75"
        / "v75_self_deciding_acceptance_cortex_certificate.json"
    )
    scoring = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V75"
        / "v75_self_deciding_acceptance_cortex_scoring_report.json"
    )
    rows = scoring["rows"]
    assert cert["real_physical_calibration_inputs_used"] is True
    assert cert["physical_backend_available"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert all(row["real_physical_calibration_inputs_used"] is True for row in rows)
    assert all(row["real_physical_calibration_hash"] == cert["real_physical_calibration_hash"] for row in rows)
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["real_physical_calibration_row_count"] == 8 for row in rows)


def test_v75_dashboard_reports_required_fields() -> None:
    dashboard = _read(
        ROOT
        / "data"
        / "protein_esperanto_engine"
        / "V75"
        / "v75_self_deciding_acceptance_cortex_dashboard.json"
    )
    total = dashboard["shards"]["TOTAL"]
    for key in [
        "targets_total",
        "accepted_count",
        "accepted_supported",
        "clean_abstain_supported",
        "blocked_supported",
        "failed_accepted",
        "accepted_accuracy",
        "coverage",
        "unknown_word_abstentions",
        "top_missing_esperanto_word",
    ]:
        assert key in total
    assert total["failed_accepted"] == 0
    assert total["accepted_accuracy"] == pytest.approx(1.0)
