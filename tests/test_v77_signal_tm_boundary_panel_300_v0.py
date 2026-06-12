from __future__ import annotations

import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V77_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V77"
E71_CERT = ROOT / "data" / "protein_esperanto_engine" / "E71" / "e71_signal_peptide_true_tm_routing_grammar_certificate.json"


EXPECTED_COMPOSITION = {
    "SIGNAL_PEPTIDE_POSITIVE": 50,
    "TRUE_TM_NEGATIVE": 60,
    "SIGNAL_ANCHOR_DECOY_ABSTAIN": 10,
    "SECRETORY_DISULFIDE_SENTINEL": 40,
    "MEMBRANE_SENTINEL": 40,
    "COILED_REPEAT_KNOTTED_MISSING_WORD_ABSTAIN_CONTROL": 40,
    "SOLUBLE_GLOBULAR_SENTINEL": 60,
}


def _read(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def test_v77_manifest_is_real_300_target_signal_tm_panel() -> None:
    manifest = _read(V77_ROOT / "v77_signal_tm_target_manifest.json")
    assert manifest["batch_id"] == "V77_SIGNAL_TM_BOUNDARY_PANEL_300"
    assert manifest["engine_version_used"] == "E71"
    assert manifest["target_count_selected"] == 300
    assert manifest["composition_rule"] == EXPECTED_COMPOSITION
    assert Counter(row["panel_group"] for row in manifest["selected_targets"]) == EXPECTED_COMPOSITION
    assert manifest["matched_control_dominance_acceptance"] is True
    assert manifest["fixed_accepted_support_thresholds_used"] is False
    assert all(row["protein_id"] and row["sequence"] for row in manifest["selected_targets"])


def test_v77_certificate_passes_with_zero_failed_accepted_and_no_static_thresholds() -> None:
    cert = _read(V77_ROOT / "v77_signal_tm_certificate.json")
    assert cert["status"] == "V77_E71_SIGNAL_TM_BOUNDARY_PANEL_PASSED"
    assert cert["targets_total"] == 300
    assert cert["accepted_count"] == 250
    assert cert["accepted_supported"] == 250
    assert cert["clean_abstain_supported"] == 50
    assert cert["failed_accepted_count"] == 0
    assert cert["accepted_accuracy"] == 1.0
    assert cert["controls_passed"] is True
    assert cert["sentinel_regressions"] == 0
    assert cert["withheld_context_leakage_detected"] is False
    assert cert["uses_static_observable_thresholds"] is False
    assert cert["matched_control_dominance_acceptance"] is True
    assert cert["real_physical_calibration_inputs_used"] is True
    assert cert["real_physical_calibration_row_count"] > 0
    assert cert["real_physical_calibration_target_native_contacts_used_before_prediction"] is False
    assert cert["real_physical_calibration_coordinate_truth_used_as_prediction_input"] is False
    assert cert["physical_basis_claim_allowed"] is False


def test_v77_scoring_separates_accepted_supported_from_clean_abstain_supported() -> None:
    rows = _read(V77_ROOT / "v77_signal_tm_scoring_report.json")["rows"]
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]
    abstained = [row for row in rows if row["acceptance_decision"] == "abstain_recommended"]
    assert len(accepted) == 250
    assert len(abstained) == 50
    assert all(row["accepted_supported"] for row in accepted)
    assert all(row["matched_control_dominance_passed"] for row in accepted)
    assert all(row["clean_abstain_supported"] for row in abstained)
    assert all(row["matched_control_dominance"]["uses_static_observable_thresholds"] is False for row in rows)


def test_v77_signal_positive_controls_are_required_for_acceptance() -> None:
    rows = _read(V77_ROOT / "v77_signal_tm_scoring_report.json")["rows"]
    signal_rows = [row for row in rows if row["panel_group"] == "SIGNAL_PEPTIDE_POSITIVE"]
    assert len(signal_rows) == EXPECTED_COMPOSITION["SIGNAL_PEPTIDE_POSITIVE"]
    required_controls = {
        "real_sequence_beats_n_terminal_shuffled_control",
        "real_sequence_beats_n_terminal_hydrophobic_masked_control",
        "real_metadata_beats_metadata_masked_source",
        "real_signal_beats_true_tm_decoy",
        "real_signal_beats_signal_anchor_decoy",
        "signal_peptide_perturbation_changes_route_more_than_true_tm_control",
        "true_tm_control_does_not_collapse_into_signal_peptide",
        "signal_anchor_decoy_abstains",
        "wrong_grammar_challenge_fails",
    }
    for row in signal_rows:
        controls = {control["control"] for control in row["matched_control_dominance"]["control_rows"]}
        assert required_controls.issubset(controls)
        assert row["known_signal_peptide_word"] == "signal_peptide_vs_true_TM"
        assert row["predicted_mechanism_class"] == "signal_peptide_vs_true_tm_routing"
        assert row["accepted_supported"] is True


def test_v77_negative_and_sentinel_groups_do_not_regress() -> None:
    rows = _read(V77_ROOT / "v77_signal_tm_scoring_report.json")["rows"]
    by_group = {group: [row for row in rows if row["panel_group"] == group] for group in EXPECTED_COMPOSITION}
    assert all(row["predicted_mechanism_class"] == "membrane_multidomain_folding_proteostasis" and row["accepted_supported"] for row in by_group["TRUE_TM_NEGATIVE"])
    assert all(row["predicted_mechanism_class"] == "membrane_multidomain_folding_proteostasis" and row["accepted_supported"] for row in by_group["MEMBRANE_SENTINEL"])
    assert all(row["predicted_mechanism_class"] == "secretory_disulfide_redox_topology" and row["accepted_supported"] for row in by_group["SECRETORY_DISULFIDE_SENTINEL"])
    assert all(row["clean_abstain_supported"] for row in by_group["SIGNAL_ANCHOR_DECOY_ABSTAIN"])
    assert all(row["clean_abstain_supported"] for row in by_group["COILED_REPEAT_KNOTTED_MISSING_WORD_ABSTAIN_CONTROL"])


def test_e71_certificate_points_to_v77_proof() -> None:
    cert = _read(E71_CERT)
    assert cert["proof_batch"] == "V77_SIGNAL_TM_BOUNDARY_PANEL_300"
    assert cert["proof_status"] == "V77_E71_SIGNAL_TM_BOUNDARY_PANEL_PASSED"
    assert cert["v77_failed_accepted_count"] == 0
    assert cert["v77_withheld_context_leakage_detected"] is False
