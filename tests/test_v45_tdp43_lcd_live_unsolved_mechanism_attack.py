from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_v45_tdp43_lcd_live_unsolved_sources_v0.py"
RUNNER = ROOT / "scripts" / "run_v45_tdp43_lcd_live_unsolved_mechanism_attack_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_builder_records_tdp43_sources_without_opening_holdouts() -> None:
    builder = _load(BUILDER, "v45_builder_sources")
    summary = builder.build_sources()
    assert summary["sequence_length"] == 141
    assert summary["gnqsyfm_fraction"] > 0.75
    assert summary["aromatic_count"] == 12
    assert summary["prediction_input_source_count"] == 3
    assert summary["holdouts_created"] is False
    manifest = json.loads((ROOT / "data" / "live_unsolved_targets" / "V45" / "TDP43_LCD" / "sources" / "source_manifest.json").read_text(encoding="utf-8"))
    assert manifest["prediction_inputs_separated_from_validation_holdouts"] is True
    assert manifest["answer_key_available_to_prediction"] is False


def test_v45_passes_live_solution_packet_but_not_absolute_solved_flag() -> None:
    builder = _load(BUILDER, "v45_builder_pass")
    runner = _load(RUNNER, "v45_runner_pass")
    builder.build_sources()
    paths = runner.run_v45()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["control_status"] == "V45_TDP43_LCD_LIVE_UNSOLVED_SOLUTION_PACKET_PASSED_REVIEW_REQUIRED"
    assert cert["live_unsolved_target_solution_packet"] is True
    assert cert["protein_folding_solved_candidate_strengthened"] is True
    assert cert["folding_problem_solved"] is False
    assert cert["prediction_sealed_before_holdout"] is True
    assert cert["claim_allowed"] is True


def test_v45_operator_buckets_are_tdp_specific_and_not_fus_transfer() -> None:
    builder = _load(BUILDER, "v45_builder_buckets")
    runner = _load(RUNNER, "v45_runner_buckets")
    builder.build_sources()
    cert = json.loads(runner.run_v45()["certificate"].read_text(encoding="utf-8"))
    required = {
        "low_complexity_disorder_operator",
        "glycine_rich_prion_like_lcd_operator",
        "sparse_aromatic_sticker_operator",
        "methionine_hydrophobic_cluster_operator",
        "alpha_helical_llps_segment_operator",
        "phosphorylation_pathology_shift_operator",
        "LLPS_self_association_operator",
        "amyloid_or_solid_state_maturation_operator",
        "nucleic_acid_context_switch_operator",
        "disease_mutation_modulation_operator",
    }
    forbidden = {
        "compact_single_native_fold_operator",
        "solved_atomic_structure_operator",
        "de_novo_global_coordinate_solution_operator",
        "AlphaFold_confidence_proxy_operator",
        "coordinate_contact_operator",
        "fus_lc_tyrosine_ladder_transfer_operator",
    }
    assert required.issubset(set(cert["operator_buckets"]))
    assert not forbidden.intersection(cert["operator_buckets"])
    assert cert["forbidden_single_fold_rejection_passed"] is True
    assert cert["fus_transfer_rejection_passed"] is True


def test_v45_perturbation_support_and_context_support() -> None:
    builder = _load(BUILDER, "v45_builder_support")
    runner = _load(RUNNER, "v45_runner_support")
    builder.build_sources()
    cert = json.loads(runner.run_v45()["certificate"].read_text(encoding="utf-8"))
    assert cert["perturbation_prediction_count"] >= 9
    assert cert["validated_prediction_count"] >= 6
    assert cert["contradicted_prediction_count"] <= 2
    assert cert["operator_support_rate"] >= 0.90
    assert cert["perturbation_support_rate"] >= 0.85
    assert cert["context_shift_support_rate"] >= 0.80


def test_v45_leakage_and_shortcut_controls_pass() -> None:
    builder = _load(BUILDER, "v45_builder_controls")
    runner = _load(RUNNER, "v45_runner_controls")
    builder.build_sources()
    cert = json.loads(runner.run_v45()["certificate"].read_text(encoding="utf-8"))
    assert cert["coordinate_derived_source_count_before_prediction"] == 0
    assert cert["internal_runtime_source_count_for_prediction"] == 0
    assert cert["holdout_leakage_detected"] is False
    assert cert["native_metrics_used_before_prediction"] is False
    assert cert["coordinate_truth_used_before_prediction"] is False
    assert cert["passed_control_count"] == cert["control_count"]
    control_ids = {control["control_id"] for control in cert["controls"]}
    assert "fus_lc_answer_transfer_blocked" in control_ids
    assert "generic_idp_annotation_alone_cannot_make_full_packet" in control_ids
    assert "nucleic_acid_context_removed_weakens_context_support" in control_ids


def test_v45_status_degrades_for_partial_generic_or_fus_transfer_cases() -> None:
    builder = _load(BUILDER, "v45_builder_degrade")
    runner = _load(RUNNER, "v45_runner_degrade")
    builder.build_sources()
    manifest = runner.load_source_manifest()
    packet = runner.predict_solution_packet(manifest)
    holdouts = runner._postseal_holdouts(packet)
    validation = runner._validate_predictions(packet, holdouts)
    partial_validation = {**validation, "holdout_source_classes": validation["holdout_source_classes"][:3]}
    assert runner._status_for_conditions(packet, partial_validation) == "V45_TDP43_LCD_LIVE_UNSOLVED_PARTIAL_CLEAN_ABSTAIN"
    generic_packet = runner.predict_solution_packet({
        **manifest,
        "prediction_sources": [{"source_id": "GENERIC_IDP_ANNOTATION_ONLY", "allowed_use": "prediction_input_before_sealing"}],
    })
    assert generic_packet["mechanism_class"] == "insufficient_evidence_clean_abstain"
    assert runner._status_for_conditions(generic_packet, validation) == "V45_TDP43_LCD_LIVE_UNSOLVED_FAILED_PREDICTIONS"
    fus_packet = runner.predict_solution_packet({
        **manifest,
        "prediction_sources": [{"source_id": "FUS_LC_SEQUENCE_PATTERN_FEATURES_DERIVED_FROM_ALLOWED_INPUTS", "allowed_use": "prediction_input_before_sealing"}],
    })
    assert fus_packet["mechanism_class"] == "insufficient_evidence_clean_abstain"


def test_v45_certificate_has_required_fields() -> None:
    builder = _load(BUILDER, "v45_builder_fields")
    runner = _load(RUNNER, "v45_runner_fields")
    builder.build_sources()
    cert = json.loads(runner.run_v45()["certificate"].read_text(encoding="utf-8"))
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
