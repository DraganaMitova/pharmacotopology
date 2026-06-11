from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_v48_sars2_orf6_viral_accessory_sources_v0.py"
RUNNER = ROOT / "scripts" / "run_v48_sars2_orf6_viral_accessory_host_hijacking_attack_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_builder_records_orf6_sources_without_opening_holdouts() -> None:
    builder = _load(BUILDER, "v48_builder_sources")
    summary = builder.build_sources()
    assert summary["sequence_length"] == 61
    assert summary["c_terminal_region"] == "38-61"
    assert summary["c_terminal_sequence"].endswith("QPMEID")
    assert summary["minimal_c_terminal_motif"] == "SQLDEEQPMEID"
    assert summary["region_window_count"] == 5
    assert summary["prediction_input_source_count"] == 4
    assert summary["holdouts_created"] is False
    manifest = json.loads((ROOT / "data" / "live_unsolved_targets" / "V48" / "SARS2_ORF6" / "sources" / "source_manifest.json").read_text(encoding="utf-8"))
    assert manifest["prediction_inputs_separated_from_validation_holdouts"] is True
    assert manifest["answer_key_available_to_prediction"] is False
    assert manifest["coordinate_sources_available_before_prediction"] is False


def test_v48_passes_live_viral_packet_but_not_absolute_solved_flag() -> None:
    builder = _load(BUILDER, "v48_builder_pass")
    runner = _load(RUNNER, "v48_runner_pass")
    builder.build_sources()
    paths = runner.run_v48()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["control_status"] == "V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_PASSED_REVIEW_REQUIRED"
    assert cert["live_viral_host_hijacking_solution_packet"] is True
    assert cert["protein_folding_solved_candidate_strengthened"] is True
    assert cert["folding_problem_solved"] is False
    assert cert["prediction_sealed_before_holdout"] is True
    assert cert["v47_committed"] is True
    assert cert["v47_commit_hash"]
    assert cert["v48_committed"] is False


def test_v48_operator_buckets_are_orf6_specific() -> None:
    builder = _load(BUILDER, "v48_builder_buckets")
    runner = _load(RUNNER, "v48_runner_buckets")
    builder.build_sources()
    cert = json.loads(runner.run_v48()["certificate"].read_text(encoding="utf-8"))
    required = {
        "C_terminal_host_interaction_operator",
        "RAE1_NUP98_binding_context_operator",
        "nuclear_transport_disruption_operator",
        "interferon_antagonism_context_operator",
        "short_linear_motif_or_MoRF_operator",
        "disorder_to_interface_operator",
        "localization_context_operator",
        "no_globular_single_fold_operator",
    }
    forbidden = {
        "compact_single_native_fold_operator",
        "generic_viral_accessory_annotation_only_operator",
        "membrane_channel_operator",
        "IDP_phase_separation_operator",
        "metamorphic_alpha_beta_fold_switch_operator",
        "solved_atomic_structure_operator",
        "coordinate_contact_operator",
        "AlphaFold_confidence_proxy_operator",
    }
    assert required.issubset(set(cert["operator_buckets"]))
    assert not forbidden.intersection(cert["operator_buckets"])
    assert cert["forbidden_globular_fold_rejection_passed"] is True
    assert cert["forbidden_generic_viral_annotation_rejection_passed"] is True


def test_v48_support_rates_and_counts_pass() -> None:
    builder = _load(BUILDER, "v48_builder_support")
    runner = _load(RUNNER, "v48_runner_support")
    builder.build_sources()
    cert = json.loads(runner.run_v48()["certificate"].read_text(encoding="utf-8"))
    assert cert["perturbation_prediction_count"] >= 10
    assert cert["validated_prediction_count"] >= 7
    assert cert["contradicted_prediction_count"] <= 2
    assert cert["operator_support_rate"] >= 0.85
    assert cert["perturbation_support_rate"] >= 0.80
    assert cert["host_interaction_support_rate"] >= 0.60
    assert cert["functional_consequence_support_rate"] >= 0.60


def test_v48_leakage_and_shortcut_controls_pass() -> None:
    builder = _load(BUILDER, "v48_builder_controls")
    runner = _load(RUNNER, "v48_runner_controls")
    builder.build_sources()
    cert = json.loads(runner.run_v48()["certificate"].read_text(encoding="utf-8"))
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
        "generic_viral_protein_annotation_alone_cannot_make_full_packet",
        "compact_single_fold_forcing_blocked",
        "removing_c_terminal_evidence_weakens_support",
        "removing_rae1_nup98_evidence_weakens_support",
        "removing_nuclear_transport_evidence_weakens_functional_grammar",
        "swapped_orf8_orf3a_nsp_evidence_does_not_validate_orf6_specific_predictions",
        "failed_predictions_remain_failed_not_repaired",
        "fewer_than_four_holdout_classes_status_partial",
        "all_generic_viral_host_claims_status_fails",
        "ifn_antagonism_not_promoted_to_atomic_fold_proof",
        "orf6_not_overpromoted_into_solved_globular_structure",
    }
    assert required_controls.issubset(control_ids)


def test_v48_status_degrades_for_partial_generic_or_cterm_only_cases() -> None:
    builder = _load(BUILDER, "v48_builder_degrade")
    runner = _load(RUNNER, "v48_runner_degrade")
    builder.build_sources()
    manifest = runner.load_source_manifest()
    packet = runner.predict_solution_packet(manifest)
    holdouts = runner._postseal_holdouts(packet)
    validation = runner._validate_predictions(packet, holdouts)
    partial_validation = {**validation, "holdout_source_classes": validation["holdout_source_classes"][:3]}
    assert runner._status_for_conditions(packet, partial_validation) == "V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_PARTIAL_CLEAN_ABSTAIN"
    generic_packet = runner.predict_solution_packet({
        **manifest,
        "prediction_sources": [{"source_id": "GENERIC_VIRAL_ACCESSORY_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    assert generic_packet["mechanism_class"] == "insufficient_evidence_clean_abstain"
    assert runner._status_for_conditions(generic_packet, validation) == "V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_FAILED_PREDICTIONS"
    cterm_only_packet = runner.predict_solution_packet({
        **manifest,
        "prediction_sources": [{"source_id": "ORF6_C_TERMINAL_ONLY_NO_HOST", "allowed_use": "prediction_input_before_sealing"}],
    })
    no_cterm_packet = runner.predict_solution_packet({
        **manifest,
        "prediction_sources": [{"source_id": "ORF6_RAE1_NUP98_ONLY_NO_CTERM", "allowed_use": "prediction_input_before_sealing"}],
    })
    assert cterm_only_packet["mechanism_class"] == "insufficient_evidence_clean_abstain"
    assert no_cterm_packet["mechanism_class"] == "insufficient_evidence_clean_abstain"


def test_v48_certificate_has_required_fields() -> None:
    builder = _load(BUILDER, "v48_builder_fields")
    runner = _load(RUNNER, "v48_runner_fields")
    builder.build_sources()
    cert = json.loads(runner.run_v48()["certificate"].read_text(encoding="utf-8"))
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
        "host_interaction_support_rate",
        "functional_consequence_support_rate",
        "forbidden_globular_fold_rejection_passed",
        "forbidden_generic_viral_annotation_rejection_passed",
        "live_viral_host_hijacking_solution_packet",
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
