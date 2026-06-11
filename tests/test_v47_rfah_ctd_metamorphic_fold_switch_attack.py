from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_v47_rfah_ctd_metamorphic_sources_v0.py"
RUNNER = ROOT / "scripts" / "run_v47_rfah_ctd_metamorphic_fold_switch_attack_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_builder_records_rfah_sources_without_opening_holdouts() -> None:
    builder = _load(BUILDER, "v47_builder_sources")
    summary = builder.build_sources()
    assert summary["sequence_length"] == 162
    assert summary["ctd_length"] == 62
    assert summary["ctd_sequence"].startswith("PKDIVD")
    assert summary["domain_count"] == 2
    assert summary["prediction_input_source_count"] == 4
    assert summary["holdouts_created"] is False
    manifest = json.loads((ROOT / "data" / "live_unsolved_targets" / "V47" / "RfaH_CTD" / "sources" / "source_manifest.json").read_text(encoding="utf-8"))
    assert manifest["prediction_inputs_separated_from_validation_holdouts"] is True
    assert manifest["answer_key_available_to_prediction"] is False
    assert manifest["coordinate_sources_available_before_prediction"] is False


def test_v47_passes_live_fold_switch_packet_but_not_absolute_solved_flag() -> None:
    builder = _load(BUILDER, "v47_builder_pass")
    runner = _load(RUNNER, "v47_runner_pass")
    builder.build_sources()
    paths = runner.run_v47()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["control_status"] == "V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_PASSED_REVIEW_REQUIRED"
    assert cert["live_fold_switch_solution_packet"] is True
    assert cert["protein_folding_solved_candidate_strengthened"] is True
    assert cert["folding_problem_solved"] is False
    assert cert["prediction_sealed_before_holdout"] is True
    assert cert["v46_committed"] is True
    assert cert["v46_commit_hash"]
    assert cert["v47_committed"] is False


def test_v47_operator_buckets_are_metamorphic_specific() -> None:
    builder = _load(BUILDER, "v47_builder_buckets")
    runner = _load(RUNNER, "v47_runner_buckets")
    builder.build_sources()
    cert = json.loads(runner.run_v47()["certificate"].read_text(encoding="utf-8"))
    required = {
        "alpha_helical_hairpin_state_operator",
        "beta_barrel_or_beta_roll_state_operator",
        "NTD_bound_autoinhibition_operator",
        "CTD_release_switch_operator",
        "transcription_activation_context_operator",
        "translation_coupling_context_operator",
        "partner_context_refolding_operator",
        "no_single_consensus_fold_operator",
    }
    forbidden = {
        "compact_single_native_fold_operator",
        "generic_two_state_annotation_only_operator",
        "intrinsic_disorder_phase_separation_operator",
        "membrane_channel_operator",
        "solved_atomic_structure_operator",
        "coordinate_contact_operator",
        "AlphaFold_confidence_proxy_operator",
    }
    assert required.issubset(set(cert["operator_buckets"]))
    assert not forbidden.intersection(cert["operator_buckets"])
    assert cert["forbidden_single_fold_rejection_passed"] is True
    assert cert["forbidden_wrong_class_rejection_passed"] is True


def test_v47_support_rates_and_counts_pass() -> None:
    builder = _load(BUILDER, "v47_builder_support")
    runner = _load(RUNNER, "v47_runner_support")
    builder.build_sources()
    cert = json.loads(runner.run_v47()["certificate"].read_text(encoding="utf-8"))
    assert cert["perturbation_prediction_count"] >= 10
    assert cert["validated_prediction_count"] >= 7
    assert cert["contradicted_prediction_count"] <= 2
    assert cert["operator_support_rate"] >= 0.85
    assert cert["perturbation_support_rate"] >= 0.80
    assert cert["state_separation_support_rate"] >= 0.70
    assert cert["partner_context_support_rate"] >= 0.60


def test_v47_leakage_and_shortcut_controls_pass() -> None:
    builder = _load(BUILDER, "v47_builder_controls")
    runner = _load(RUNNER, "v47_runner_controls")
    builder.build_sources()
    cert = json.loads(runner.run_v47()["certificate"].read_text(encoding="utf-8"))
    assert cert["coordinate_derived_source_count_before_prediction"] == 0
    assert cert["internal_runtime_source_count_for_prediction"] == 0
    assert cert["holdout_leakage_detected"] is False
    assert cert["native_metrics_used_before_prediction"] is False
    assert cert["coordinate_truth_used_before_prediction"] is False
    assert cert["passed_control_count"] == cert["control_count"]
    control_ids = {control["control_id"] for control in cert["controls"]}
    required_controls = {
        "prediction_sealed_before_holdout_validation",
        "holdout_files_unavailable_to_prediction_function",
        "pdb_coordinates_before_sealing_blocked",
        "alphafold_esmfold_coordinates_before_sealing_blocked",
        "internal_runtime_source_as_evidence_blocked",
        "target_name_only_assignment_blocked",
        "generic_transcription_factor_annotation_alone_cannot_make_full_packet",
        "one_state_only_rfah_evidence_becomes_partial_or_invalid",
        "single_consensus_fold_forcing_blocked",
        "removing_alpha_state_evidence_weakens_support",
        "removing_beta_state_evidence_weakens_support",
        "removing_ntd_partner_context_weakens_switch_grammar",
        "swapped_xcl1_kaib_mad2_evidence_does_not_validate_rfah_specific_predictions",
        "failed_predictions_remain_failed_not_repaired",
        "fewer_than_four_holdout_classes_status_partial",
        "all_generic_metamorphic_claims_status_fails",
    }
    assert required_controls.issubset(control_ids)


def test_v47_status_degrades_for_partial_generic_or_one_state_cases() -> None:
    builder = _load(BUILDER, "v47_builder_degrade")
    runner = _load(RUNNER, "v47_runner_degrade")
    builder.build_sources()
    manifest = runner.load_source_manifest()
    packet = runner.predict_solution_packet(manifest)
    holdouts = runner._postseal_holdouts(packet)
    validation = runner._validate_predictions(packet, holdouts)
    partial_validation = {**validation, "holdout_source_classes": validation["holdout_source_classes"][:3]}
    assert runner._status_for_conditions(packet, partial_validation) == "V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_PARTIAL_CLEAN_ABSTAIN"
    generic_packet = runner.predict_solution_packet({
        **manifest,
        "prediction_sources": [{"source_id": "GENERIC_TRANSCRIPTION_FACTOR_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    assert generic_packet["mechanism_class"] == "insufficient_evidence_clean_abstain"
    assert runner._status_for_conditions(generic_packet, validation) == "V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_FAILED_PREDICTIONS"
    alpha_only_packet = runner.predict_solution_packet({
        **manifest,
        "prediction_sources": [{"source_id": "RFAH_ALPHA_STATE_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    beta_only_packet = runner.predict_solution_packet({
        **manifest,
        "prediction_sources": [{"source_id": "RFAH_BETA_STATE_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    assert alpha_only_packet["mechanism_class"] == "insufficient_evidence_clean_abstain"
    assert beta_only_packet["mechanism_class"] == "insufficient_evidence_clean_abstain"


def test_v47_certificate_has_required_fields() -> None:
    builder = _load(BUILDER, "v47_builder_fields")
    runner = _load(RUNNER, "v47_runner_fields")
    builder.build_sources()
    cert = json.loads(runner.run_v47()["certificate"].read_text(encoding="utf-8"))
    required = [
        "kind",
        "run_mode",
        "target",
        "sequence_region_scope",
        "mechanism_class",
        "sealed_prediction_hash",
        "prediction_sealed_before_holdout",
        "prediction_input_source_count",
        "holdout_source_count",
        "operator_bucket_count",
        "perturbation_prediction_count",
        "validated_prediction_count",
        "contradicted_prediction_count",
        "operator_support_rate",
        "perturbation_support_rate",
        "state_separation_support_rate",
        "partner_context_support_rate",
        "forbidden_single_fold_rejection_passed",
        "forbidden_wrong_class_rejection_passed",
        "live_fold_switch_solution_packet",
        "protein_folding_solved_candidate_strengthened",
        "folding_problem_solved",
        "coordinate_derived_source_count_before_prediction",
        "internal_runtime_source_count_for_prediction",
        "holdout_leakage_detected",
        "native_metrics_used_before_prediction",
        "coordinate_truth_used_before_prediction",
        "claim_allowed",
        "allowed_claim_text",
        "forbidden_claims",
        "control_count",
        "passed_control_count",
        "failed_checks",
        "next_action",
    ]
    for key in required:
        assert key in cert
    assert cert["failed_checks"] == []
