from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_v44_fus_lc_live_unsolved_sources_v0.py"
RUNNER = ROOT / "scripts" / "run_v44_fus_lc_live_unsolved_mechanism_attack_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_builder_records_fus_lc_sources_without_opening_holdouts() -> None:
    builder = _load(BUILDER, "v44_builder_sources")
    summary = builder.build_sources()
    assert summary["sequence_length"] == 214
    assert summary["qgsy_fraction"] > 0.80
    assert summary["prediction_input_source_count"] == 3
    assert summary["holdouts_created"] is False
    manifest = json.loads((ROOT / "data" / "live_unsolved_targets" / "V44" / "FUS_LC" / "sources" / "source_manifest.json").read_text(encoding="utf-8"))
    assert manifest["prediction_inputs_separated_from_validation_holdouts"] is True
    assert manifest["answer_key_available_to_prediction"] is False


def test_v44_passes_live_solution_packet_but_not_absolute_solved_flag() -> None:
    builder = _load(BUILDER, "v44_builder_pass")
    runner = _load(RUNNER, "v44_runner_pass")
    builder.build_sources()
    paths = runner.run_v44()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["control_status"] == "V44_FUS_LC_LIVE_UNSOLVED_SOLUTION_PACKET_PASSED_REVIEW_REQUIRED"
    assert cert["live_unsolved_target_solution_packet"] is True
    assert cert["protein_folding_solved_candidate_strengthened"] is True
    assert cert["folding_problem_solved"] is False
    assert cert["prediction_sealed_before_holdout"] is True
    assert cert["claim_allowed"] is True


def test_v44_operator_buckets_and_forbidden_buckets() -> None:
    builder = _load(BUILDER, "v44_builder_buckets")
    runner = _load(RUNNER, "v44_runner_buckets")
    builder.build_sources()
    cert = json.loads(runner.run_v44()["certificate"].read_text(encoding="utf-8"))
    required = {
        "low_complexity_disorder_operator",
        "aromatic_sticker_pattern_operator",
        "polar_spacer_context_operator",
        "phosphorylation_shift_operator",
        "LLPS_self_association_operator",
        "fibril_or_gel_state_shift_operator",
        "context_bound_RNA_or_protein_interaction_operator",
    }
    forbidden = {
        "compact_single_native_fold_operator",
        "solved_atomic_structure_operator",
        "de_novo_global_coordinate_solution_operator",
        "AlphaFold_confidence_proxy_operator",
        "coordinate_contact_operator",
    }
    assert required.issubset(set(cert["operator_buckets"]))
    assert not forbidden.intersection(cert["operator_buckets"])
    assert cert["forbidden_single_fold_rejection_passed"] is True


def test_v44_perturbation_support_and_context_support() -> None:
    builder = _load(BUILDER, "v44_builder_support")
    runner = _load(RUNNER, "v44_runner_support")
    builder.build_sources()
    cert = json.loads(runner.run_v44()["certificate"].read_text(encoding="utf-8"))
    assert cert["perturbation_prediction_count"] >= 8
    assert cert["validated_prediction_count"] >= 5
    assert cert["contradicted_prediction_count"] <= 2
    assert cert["operator_support_rate"] >= 0.85
    assert cert["perturbation_support_rate"] >= 0.80
    assert cert["context_shift_support_rate"] >= 0.75


def test_v44_leakage_and_shortcut_controls_pass() -> None:
    builder = _load(BUILDER, "v44_builder_controls")
    runner = _load(RUNNER, "v44_runner_controls")
    builder.build_sources()
    cert = json.loads(runner.run_v44()["certificate"].read_text(encoding="utf-8"))
    assert cert["coordinate_derived_source_count_before_prediction"] == 0
    assert cert["internal_runtime_source_count_for_prediction"] == 0
    assert cert["holdout_leakage_detected"] is False
    assert cert["native_metrics_used_before_prediction"] is False
    assert cert["coordinate_truth_used_before_prediction"] is False
    assert cert["passed_control_count"] == cert["control_count"]
    control_ids = {control["control_id"] for control in cert["controls"]}
    assert "target_name_only_assignment_blocked" in control_ids
    assert "generic_idp_annotation_alone_cannot_make_full_packet" in control_ids
    assert "swapped_snca_tdp43_htt_evidence_does_not_validate_fus_specific_predictions" in control_ids


def test_v44_status_degrades_for_partial_or_generic_cases() -> None:
    builder = _load(BUILDER, "v44_builder_degrade")
    runner = _load(RUNNER, "v44_runner_degrade")
    builder.build_sources()
    manifest = runner.load_source_manifest()
    packet = runner.predict_solution_packet(manifest)
    holdouts = runner._postseal_holdouts(packet)
    validation = runner._validate_predictions(packet, holdouts)
    partial_validation = {**validation, "holdout_source_classes": validation["holdout_source_classes"][:3]}
    assert runner._status_for_conditions(packet, partial_validation) == "V44_FUS_LC_LIVE_UNSOLVED_PARTIAL_CLEAN_ABSTAIN"
    generic_packet = runner.predict_solution_packet({
        **manifest,
        "prediction_sources": [{"source_id": "GENERIC_IDP_ANNOTATION_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    assert generic_packet["mechanism_class"] == "insufficient_evidence_clean_abstain"
    assert runner._status_for_conditions(generic_packet, validation) == "V44_FUS_LC_LIVE_UNSOLVED_FAILED_PREDICTIONS"


def test_v44_certificate_has_required_fields() -> None:
    builder = _load(BUILDER, "v44_builder_fields")
    runner = _load(RUNNER, "v44_runner_fields")
    builder.build_sources()
    cert = json.loads(runner.run_v44()["certificate"].read_text(encoding="utf-8"))
    required = [
        "kind",
        "run_mode",
        "target",
        "sequence_region_scope",
        "mechanism_class",
        "sealed_prediction_hash",
        "prediction_input_source_count",
        "holdout_source_count",
        "operator_bucket_count",
        "perturbation_prediction_count",
        "validated_prediction_count",
        "operator_support_rate",
        "perturbation_support_rate",
        "context_shift_support_rate",
        "live_unsolved_target_solution_packet",
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
