from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
V76_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V76"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v76_uses_expected_mixed_repair_panel_composition() -> None:
    manifest = _read(V76_ROOT / "v76_secretory_disulfide_target_manifest.json")
    expected = {
        "SECRETORY_DISULFIDE_POSITIVE": 60,
        "METAL_CYS_HIS_NEGATIVE": 30,
        "SIGNAL_TM_NEGATIVE": 20,
        "CYSTEINE_GLOBULAR_NEGATIVE": 20,
        "BETA_CYSTEINE_NEGATIVE": 20,
        "COILED_REPEAT_MISSING_WORD_ABSTAIN_CONTROL": 20,
        "V75_ACCEPTED_SENTINEL_REPLAY": 30,
    }

    assert manifest["batch_id"] == "V76_SECRETORY_DISULFIDE_REPAIR_PANEL_200"
    assert manifest["engine_version_used"] == "E70"
    assert manifest["baseline_engine_version"] == "E69"
    assert manifest["target_count_selected"] == 200
    assert manifest["composition_rule"] == expected
    assert manifest["matched_control_dominance_acceptance"] is True
    assert manifest["fixed_accepted_support_thresholds_used"] is False
    assert Counter(row["panel_group"] for row in manifest["selected_targets"]) == expected


def test_v76_passes_zero_failed_accepted_with_matched_control_dominance() -> None:
    cert = _read(V76_ROOT / "v76_secretory_disulfide_certificate.json")

    assert cert["status"] == "V76_E70_SECRETORY_DISULFIDE_REPAIR_PANEL_PASSED"
    assert cert["targets_total"] == 200
    assert cert["accepted_count"] == 160
    assert cert["accepted_supported"] == 160
    assert cert["clean_abstain_supported"] == 40
    assert cert["failed_accepted"] == 0
    assert cert["accepted_accuracy"] == pytest.approx(1.0)
    assert cert["coverage"] == pytest.approx(0.8)
    assert cert["matched_control_dominance_acceptance"] is True
    assert cert["fixed_accepted_support_thresholds_used"] is False
    assert cert["secretory_disulfide_positive_supported"] == 60
    assert cert["secretory_disulfide_positive_total"] == 60
    assert cert["missing_word_controls_cleanly_abstained"] == 20
    assert cert["sentinel_regressions"] == 0
    assert cert["controls_passed"] is True
    assert cert["real_physical_calibration_inputs_used"] is True
    assert cert["real_physical_calibration_source_coordinate_database"] == "RCSB_PDB"
    assert cert["real_physical_calibration_target_native_excluded"] is True
    assert cert["real_physical_calibration_target_native_contacts_used_before_prediction"] is False
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["folding_problem_solved"] is False
    assert cert["next_required_batch"] == "V76P_OPENMM_COARSE_TARGET_EXECUTION_PILOT"


def test_v76_disulfide_rows_use_paired_controls_not_static_thresholds() -> None:
    scoring = _read(V76_ROOT / "v76_secretory_disulfide_scoring_report.json")
    rows = scoring["rows"]
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]
    disulfide = [row for row in rows if row["panel_group"] == "SECRETORY_DISULFIDE_POSITIVE"]

    assert len(accepted) == 160
    assert all(row["accepted_supported"] for row in accepted)
    assert all(row["matched_control_dominance_passed"] for row in accepted)
    assert all(row["matched_control_dominance"]["uses_static_observable_thresholds"] is False for row in rows)
    assert len(disulfide) == 60

    required_controls = {
        "real_sequence_beats_shuffled_no_context_control",
        "real_metadata_beats_metadata_masked_source",
        "real_sequence_beats_cysteine_masked_control",
        "disulfide_redox_perturbation_decreases_selected_state",
        "wrong_grammar_challenge_fails",
    }
    for row in disulfide:
        matched = row["matched_control_dominance"]
        by_name = {control["control"]: control for control in matched["control_rows"]}
        assert set(by_name) == required_controls
        assert all(control["passed"] for control in by_name.values())
        assert by_name["real_sequence_beats_shuffled_no_context_control"]["real_value"] > by_name["real_sequence_beats_shuffled_no_context_control"]["control_value"]
        assert by_name["real_metadata_beats_metadata_masked_source"]["real_value"] > by_name["real_metadata_beats_metadata_masked_source"]["control_value"]
        assert by_name["real_sequence_beats_cysteine_masked_control"]["real_value"] > by_name["real_sequence_beats_cysteine_masked_control"]["control_value"]
        perturbation_rows = by_name["disulfide_redox_perturbation_decreases_selected_state"]["control_value"]
        assert perturbation_rows
        assert all(item["perturbed_value"] < item["baseline_value"] for item in perturbation_rows)


def test_v76_negative_and_missing_word_controls_do_not_get_stolen_by_disulfide() -> None:
    scoring = _read(V76_ROOT / "v76_secretory_disulfide_scoring_report.json")
    rows = scoring["rows"]
    counts = Counter(
        (
            row["panel_group"],
            row["acceptance_decision"],
            row["predicted_mechanism_class"],
            row["missing_word_candidate"],
        )
        for row in rows
    )

    assert counts[("METAL_CYS_HIS_NEGATIVE", "accepted", "metal_cluster_and_ligand_locked_basin", None)] == 30
    assert counts[("CYSTEINE_GLOBULAR_NEGATIVE", "accepted", "globular_closure", None)] == 20
    assert counts[("BETA_CYSTEINE_NEGATIVE", "accepted", "beta_closure_topology", None)] == 20
    assert counts[("V75_ACCEPTED_SENTINEL_REPLAY", "accepted", "multidomain_allosteric_architecture", None)] == 30
    assert counts[("SIGNAL_TM_NEGATIVE", "abstain_recommended", "membrane_multidomain_folding_proteostasis", "signal_peptide_vs_true_TM")] == 20
    assert counts[("COILED_REPEAT_MISSING_WORD_ABSTAIN_CONTROL", "abstain_recommended", "insufficient_evidence_clean_abstain", "coiled_coil_register")] == 12
    assert counts[("COILED_REPEAT_MISSING_WORD_ABSTAIN_CONTROL", "abstain_recommended", "insufficient_evidence_clean_abstain", "repeat_solenoid_topology")] == 8

    clean_abstains = [row for row in rows if row["acceptance_decision"] == "abstain_recommended"]
    assert len(clean_abstains) == 40
    assert all(row["clean_abstain_supported"] for row in clean_abstains)
    assert all(row["final_self_decision"] == "clean_abstain_missing_word" for row in clean_abstains)
