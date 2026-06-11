from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_v46_cftr_f508del_membrane_multidomain_sources_v0.py"
RUNNER = ROOT / "scripts" / "run_v46_cftr_f508del_membrane_multidomain_attack_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_builder_records_cftr_sources_without_opening_holdouts() -> None:
    builder = _load(BUILDER, "v46_builder_sources")
    summary = builder.build_sources()
    assert summary["sequence_length"] == 1480
    assert summary["f508_residue"] == "F"
    assert summary["domain_count"] == 5
    assert summary["transmembrane_segment_count"] == 12
    assert summary["prediction_input_source_count"] == 3
    assert summary["holdouts_created"] is False
    manifest = json.loads((ROOT / "data" / "live_unsolved_targets" / "V46" / "CFTR_F508del" / "sources" / "source_manifest.json").read_text(encoding="utf-8"))
    assert manifest["prediction_inputs_separated_from_validation_holdouts"] is True
    assert manifest["answer_key_available_to_prediction"] is False


def test_v46_passes_live_membrane_packet_but_not_absolute_solved_flag() -> None:
    builder = _load(BUILDER, "v46_builder_pass")
    runner = _load(RUNNER, "v46_runner_pass")
    builder.build_sources()
    paths = runner.run_v46()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["control_status"] == "V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_SOLUTION_PACKET_PASSED_REVIEW_REQUIRED"
    assert cert["live_membrane_solution_packet"] is True
    assert cert["protein_folding_solved_candidate_strengthened"] is True
    assert cert["folding_problem_solved"] is False
    assert cert["prediction_sealed_before_holdout"] is True


def test_v46_operator_buckets_are_membrane_multidomain_specific() -> None:
    builder = _load(BUILDER, "v46_builder_buckets")
    runner = _load(RUNNER, "v46_runner_buckets")
    builder.build_sources()
    cert = json.loads(runner.run_v46()["certificate"].read_text(encoding="utf-8"))
    required = {
        "membrane_domain_operator",
        "NBD1_stability_operator",
        "F508del_local_destabilization_operator",
        "interdomain_interface_coupling_operator",
        "trafficking_quality_control_operator",
        "corrector_or_rescue_context_operator",
        "multidomain_assembly_operator",
        "channel_function_context_operator",
    }
    forbidden = {
        "generic_channel_annotation_only_operator",
        "single_local_mutation_only_operator",
        "compact_single_domain_fold_operator",
        "solved_atomic_structure_operator",
        "coordinate_contact_operator",
        "AlphaFold_confidence_proxy_operator",
    }
    assert required.issubset(set(cert["operator_buckets"]))
    assert not forbidden.intersection(cert["operator_buckets"])
    assert cert["forbidden_single_local_explanation_rejection_passed"] is True
    assert cert["forbidden_generic_channel_rejection_passed"] is True


def test_v46_support_rates_and_counts_pass() -> None:
    builder = _load(BUILDER, "v46_builder_support")
    runner = _load(RUNNER, "v46_runner_support")
    builder.build_sources()
    cert = json.loads(runner.run_v46()["certificate"].read_text(encoding="utf-8"))
    assert cert["perturbation_prediction_count"] >= 10
    assert cert["validated_prediction_count"] >= 7
    assert cert["contradicted_prediction_count"] <= 2
    assert cert["operator_support_rate"] >= 0.85
    assert cert["perturbation_support_rate"] >= 0.80
    assert cert["interdomain_coupling_support_rate"] >= 0.60
    assert cert["rescue_logic_support_rate"] >= 0.60


def test_v46_leakage_and_shortcut_controls_pass() -> None:
    builder = _load(BUILDER, "v46_builder_controls")
    runner = _load(RUNNER, "v46_runner_controls")
    builder.build_sources()
    cert = json.loads(runner.run_v46()["certificate"].read_text(encoding="utf-8"))
    assert cert["coordinate_derived_source_count_before_prediction"] == 0
    assert cert["internal_runtime_source_count_for_prediction"] == 0
    assert cert["holdout_leakage_detected"] is False
    assert cert["native_metrics_used_before_prediction"] is False
    assert cert["coordinate_truth_used_before_prediction"] is False
    assert cert["passed_control_count"] == cert["control_count"]
    control_ids = {control["control_id"] for control in cert["controls"]}
    assert "generic_channel_annotation_alone_cannot_make_full_packet" in control_ids
    assert "single_local_mutation_only_explanation_blocked" in control_ids
    assert "swapped_kcsa_cftr_evidence_does_not_validate_specific_predictions" in control_ids


def test_v46_status_degrades_for_partial_generic_or_single_local_cases() -> None:
    builder = _load(BUILDER, "v46_builder_degrade")
    runner = _load(RUNNER, "v46_runner_degrade")
    builder.build_sources()
    manifest = runner.load_source_manifest()
    packet = runner.predict_solution_packet(manifest)
    holdouts = runner._postseal_holdouts(packet)
    validation = runner._validate_predictions(packet, holdouts)
    partial_validation = {**validation, "holdout_source_classes": validation["holdout_source_classes"][:3]}
    assert runner._status_for_conditions(packet, partial_validation) == "V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_PARTIAL_CLEAN_ABSTAIN"
    generic_packet = runner.predict_solution_packet({
        **manifest,
        "prediction_sources": [{"source_id": "GENERIC_CHANNEL_ANNOTATION_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    assert generic_packet["mechanism_class"] == "insufficient_evidence_clean_abstain"
    assert runner._status_for_conditions(generic_packet, validation) == "V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_FAILED_PREDICTIONS"
    single_packet = runner.predict_solution_packet({
        **manifest,
        "prediction_sources": [{"source_id": "UNIPROT_P13569_F508DEL_VARIANT_NONCOORDINATE", "allowed_use": "prediction_input_before_sealing"}],
    })
    assert single_packet["mechanism_class"] == "insufficient_evidence_clean_abstain"


def test_v46_certificate_has_required_fields() -> None:
    builder = _load(BUILDER, "v46_builder_fields")
    runner = _load(RUNNER, "v46_runner_fields")
    builder.build_sources()
    cert = json.loads(runner.run_v46()["certificate"].read_text(encoding="utf-8"))
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
        "interdomain_coupling_support_rate",
        "rescue_logic_support_rate",
        "forbidden_single_local_explanation_rejection_passed",
        "forbidden_generic_channel_rejection_passed",
        "live_membrane_solution_packet",
        "protein_folding_solved_candidate_strengthened",
        "folding_problem_solved",
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
