from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v37_mechanism_question_probes_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v37_probes", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_baseline_fixture_assigns_all_three_mechanism_classes() -> None:
    mod = _load_module()
    cert = mod.build_v37(mod._ready_fixture_dossiers())
    assert cert["control_status"] == "V37_MECHANISM_QUESTION_PROBES_PASSED_CLAIM_DISABLED"
    assert cert["assigned_target_count"] == 3
    assert cert["partial_target_count"] == 0
    assert cert["blocked_target_count"] == 0
    assert cert["mechanism_classes_assigned"] == {
        "KcsA": "membrane_pore_filter_oligomeric_ion_selectivity",
        "XCL1_lymphotactin": "metamorphic_two_state_fold_switch",
        "alpha_synuclein_SNCA": "intrinsic_disorder_contextual_ensemble",
    }
    assert cert["target_name_masking_passed"] is True
    assert cert["swapped_dossier_detection_passed"] is True
    assert cert["coordinate_derived_source_count"] == 0
    assert cert["internal_runtime_source_count"] == 0
    assert cert["claim_allowed"] is False
    assert cert["new_MD_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["positive_folding_evidence_found"] is False
    assert cert["folding_problem_solved"] is False


def test_required_controls_all_pass() -> None:
    mod = _load_module()
    controls = mod.run_v37_controls(mod._ready_fixture_dossiers())
    assert len(controls) == 12
    assert [row["passed"] for row in controls] == [True] * 12
    assert {row["control_id"] for row in controls} == {
        "baseline_v36_dossiers_assign_all_three_mechanism_classes",
        "target_name_masked_dossiers_still_assign_from_evidence_content",
        "swapped_dossiers_are_detected_not_renamed_into_target_grammar",
        "kcsa_without_ion_filter_evidence_is_partial",
        "xcl1_without_state_b_is_partial",
        "xcl1_forced_mixed_state_pooling_is_blocked",
        "snca_without_disorder_evidence_is_partial",
        "snca_forced_compact_single_fold_is_blocked",
        "generic_only_annotations_do_not_assign_high_confidence_grammar",
        "coordinate_derived_source_supplied_to_v37_is_blocked",
        "internal_runtime_source_supplied_to_v37_is_blocked",
        "placeholder_source_row_supplied_to_v37_is_blocked",
    }


def test_masked_target_names_still_assign_from_content() -> None:
    mod = _load_module()
    packages = mod._mask_dossier_names(mod._ready_fixture_dossiers())
    assert mod.target_name_masking_passes(packages) is True
    cert = mod.build_v37(packages)
    assert cert["control_status"] == "V37_MECHANISM_QUESTION_PROBES_PASSED_CLAIM_DISABLED"


def test_swapped_dossiers_are_detected() -> None:
    mod = _load_module()
    packages = mod._ready_fixture_dossiers()
    swapped = {
        "KcsA": packages["alpha_synuclein_SNCA"],
        "XCL1_lymphotactin": packages["XCL1_lymphotactin"],
        "alpha_synuclein_SNCA": packages["KcsA"],
    }
    swapped["KcsA"]["target"] = "KcsA"
    swapped["alpha_synuclein_SNCA"]["target"] = "alpha_synuclein_SNCA"
    cert = mod._evaluate_control_packages(swapped)
    assert cert["control_status"] == "V37_BLOCKED_FORCED_WRONG_GRAMMAR"
    assert any("mechanism_class_mismatch" in check for check in cert["failed_checks"])


def test_kcsa_ion_filter_removed_is_partial() -> None:
    mod = _load_module()
    packages = mod._ready_fixture_dossiers()
    packages["KcsA"]["bucket_status"]["present_buckets"] = ["sequence_or_family_identity", "membrane_topology_context"]
    cert = mod._evaluate_control_packages(packages)
    assert cert["control_status"] == "V37_PARTIAL_MECHANISM_ASSIGNMENT_CLEAN_ABSTAIN"
    assert any("KcsA:missing_required_operator:ion_selectivity_operator" in check for check in cert["failed_checks"])


def test_xcl1_state_b_removed_partial_and_mixed_pooling_blocked() -> None:
    mod = _load_module()
    packages = mod._ready_fixture_dossiers()
    packages["XCL1_lymphotactin"]["bucket_status"]["present_buckets"] = [
        "state_A_function_context",
        "metamorphic_two_state_context",
        "no_mixed_state_pooling_rule",
    ]
    cert = mod._evaluate_control_packages(packages)
    assert cert["control_status"] == "V37_PARTIAL_MECHANISM_ASSIGNMENT_CLEAN_ABSTAIN"
    assert any("state_B_beta_sandwich_dimer_operator" in check for check in cert["failed_checks"])

    packages = mod._ready_fixture_dossiers()
    packages["XCL1_lymphotactin"]["grammar_rules"]["mixed_state_pooling_allowed"] = True
    cert = mod._evaluate_control_packages(packages)
    assert cert["control_status"] == "V37_BLOCKED_FORCED_WRONG_GRAMMAR"
    assert any("single_consensus_fold_operator" in check for check in cert["failed_checks"])


def test_snca_disorder_removed_partial_and_single_fold_blocked() -> None:
    mod = _load_module()
    packages = mod._ready_fixture_dossiers()
    packages["alpha_synuclein_SNCA"]["bucket_status"]["present_buckets"] = [
        "ensemble_context",
        "disorder_to_order_context",
        "no_single_native_fold_rule",
    ]
    cert = mod._evaluate_control_packages(packages)
    assert cert["control_status"] == "V37_PARTIAL_MECHANISM_ASSIGNMENT_CLEAN_ABSTAIN"
    assert any("intrinsic_disorder_operator" in check for check in cert["failed_checks"])

    packages = mod._ready_fixture_dossiers()
    packages["alpha_synuclein_SNCA"]["grammar_rules"]["single_native_fold_model_allowed"] = True
    cert = mod._evaluate_control_packages(packages)
    assert cert["control_status"] == "V37_BLOCKED_FORCED_WRONG_GRAMMAR"
    assert any("compact_native_single_fold_operator" in check for check in cert["failed_checks"])


def test_generic_only_annotations_do_not_assign_high_confidence() -> None:
    mod = _load_module()
    generic = mod._ready_fixture_dossiers()
    generic["KcsA"] = mod._minimal_dossier("KcsA", ["sequence_or_family_identity"], {"native_metrics_before_selection_allowed": False}, "generic protein annotation")
    cert = mod._evaluate_control_packages(generic)
    assert cert["control_status"] == "V37_PARTIAL_MECHANISM_ASSIGNMENT_CLEAN_ABSTAIN"
    assert any("no_high_confidence_mechanism_assignment" in check for check in cert["failed_checks"])


def test_coordinate_internal_placeholder_controls_block() -> None:
    mod = _load_module()
    packages = mod._ready_fixture_dossiers()
    packages["KcsA"]["source_rows"] = [{"source_name": "1BL8 contacts", "coordinate_derived": True}]
    cert = mod._evaluate_control_packages(packages)
    assert cert["control_status"] == "V37_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"
    assert any("coordinate_or_prediction_leakage" in check for check in cert["failed_checks"])

    packages = mod._ready_fixture_dossiers()
    packages["KcsA"]["source_rows"] = [{"source_name": "runtime", "source_url_or_citation": "first_contact_clean_pharmacotopology_layer_run/V36/report.json"}]
    cert = mod._evaluate_control_packages(packages)
    assert cert["control_status"] == "V37_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"
    assert any("internal_runtime_source" in check for check in cert["failed_checks"])

    packages = mod._ready_fixture_dossiers()
    packages["KcsA"]["source_rows"] = [{"source_name": "TODO placeholder", "source_url_or_citation": "citation needed"}]
    cert = mod._evaluate_control_packages(packages)
    assert cert["control_status"] == "V37_BLOCKED_FORCED_WRONG_GRAMMAR"
    assert any("placeholder_or_missing_source_name" in check for check in cert["failed_checks"])


def test_writer_outputs_maps_and_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_v37(mod._ready_fixture_dossiers())
    paths = mod.write_outputs(tmp_path / "out", tmp_path / "maps", cert)
    for path in paths.values():
        assert path.exists()
    written = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert written["next_decision"]["claim_allowed"] is False
    assert written["artifacts"]["answer_table"] == str(paths["answer_table"])
