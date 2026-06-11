from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_v43_solved_flag_trial_panel_v0.py"
RUNNER = ROOT / "scripts" / "run_v43_protein_folding_solved_flag_trial_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_builder_imports_v42_panel_with_contact_requirements() -> None:
    builder = _load(BUILDER, "v43_builder")
    summary = builder.build_panel()
    assert summary["panel_target_count"] == 24
    assert summary["contact_scorable_target_count"] >= 12
    assert summary["protein_folding_solved_candidate"] is False
    panel = json.loads((ROOT / "data" / "solved_flag_trial" / "V43" / "panel" / "panel_manifest.json").read_text(encoding="utf-8"))
    assert panel["holdouts_created"] is False
    assert panel["answer_key_available_to_prediction"] is False


def test_v43_can_set_solved_candidate_but_not_absolute_solved_flag() -> None:
    builder = _load(BUILDER, "v43_builder_pass")
    runner = _load(RUNNER, "v43_runner_pass")
    builder.build_panel()
    paths = runner.run_v43()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["control_status"] == "V43_PROTEIN_FOLDING_SOLVED_CANDIDATE_PASSED_REVIEW_REQUIRED"
    assert cert["protein_folding_solved_candidate"] is True
    assert cert["positive_folding_evidence_found"] is True
    assert cert["claim_allowed"] is True
    assert cert["folding_problem_solved"] is False
    assert cert["prediction_sealed_before_holdout"] is True


def test_v43_thresholds_and_baselines_pass() -> None:
    builder = _load(BUILDER, "v43_builder_thresholds")
    runner = _load(RUNNER, "v43_runner_thresholds")
    builder.build_panel()
    cert = json.loads(runner.run_v43()["certificate"].read_text(encoding="utf-8"))
    assert all(row["passed"] for row in cert["thresholds"])
    assert sum(1 for row in cert["thresholds"] if row["passed"]) == 22
    assert cert["beats_random_baseline"] is True
    assert cert["beats_keyword_baseline"] is True
    assert cert["beats_majority_baseline"] is True
    assert cert["beats_simple_sequence_feature_baseline"] is True
    assert cert["top_L_long_range_contact_precision"] >= 0.40
    assert cert["contact_enrichment_vs_random"] >= 2.0


def test_v43_leakage_shortcuts_are_blocked() -> None:
    builder = _load(BUILDER, "v43_builder_controls")
    runner = _load(RUNNER, "v43_runner_controls")
    builder.build_panel()
    panel = runner.load_panel()
    target = panel["targets"][0]
    bad_coord = {**target, "prediction_sources": [{"source_id": "bad", "source_type": "PDB coordinates", "source_url_or_citation": "PDB coordinate", "coordinate_derived": True}]}
    assert runner.predict_target(bad_coord)["coordinate_derived_source_count_before_prediction"] > 0
    bad_internal = {**target, "prediction_sources": [{"source_id": "bad", "source_type": "runtime", "source_url_or_citation": "first_contact_clean_pharmacotopology_layer_run/report", "internal_runtime_source": True}]}
    assert runner.predict_target(bad_internal)["internal_runtime_source_count_for_prediction"] > 0
    name_only = {**target, "prediction_sources": []}
    assert runner.predict_target(name_only)["mechanism_class"] == "insufficient_evidence_clean_abstain"


def test_v43_disorder_and_multistate_are_not_forced_single_fold() -> None:
    builder = _load(BUILDER, "v43_builder_classes")
    runner = _load(RUNNER, "v43_runner_classes")
    builder.build_panel()
    panel = runner.load_panel()
    idp = next(t for t in panel["targets"] if t["panel_group"] == "disordered_or_partially_disordered")
    multi = next(t for t in panel["targets"] if t["panel_group"] == "multistate_allosteric_metamorphic_binding")
    membrane_generic = {**panel["targets"][0], "prediction_sources": [{"source_id": "membrane_only", "feature_tags": ["membrane"], "region_hints": []}]}
    assert runner.predict_target(idp)["mechanism_class"] != "compact_single_fold_core_closure"
    assert runner.predict_target(multi)["mechanism_class"] != "compact_single_fold_core_closure"
    assert runner.predict_target(membrane_generic)["mechanism_class"] != "membrane_pore_filter_or_transport"


def test_v43_contradicting_holdout_fails_not_rescued() -> None:
    builder = _load(BUILDER, "v43_builder_contradict")
    runner = _load(RUNNER, "v43_runner_contradict")
    builder.build_panel()
    panel = runner.load_panel()
    packet = runner.predict_target(panel["targets"][0])
    holdout = runner._postseal_holdout(panel["targets"][0], packet)
    holdout["contradicts_prediction"] = True
    result = runner._score_target(packet, holdout)
    assert result["validation_level"] == "failed_contradicted_holdout"


def test_v43_certificate_has_required_fields_and_controls() -> None:
    builder = _load(BUILDER, "v43_builder_fields")
    runner = _load(RUNNER, "v43_runner_fields")
    builder.build_panel()
    cert = json.loads(runner.run_v43()["certificate"].read_text(encoding="utf-8"))
    required = [
        "region_overlap_f1",
        "motif_or_signature_recovery_rate",
        "top_L_long_range_contact_precision",
        "top_L_over_2_long_range_contact_precision",
        "contact_enrichment_vs_random",
        "false_contact_promotion_rate",
        "compact_vs_disorder_accuracy",
        "membrane_vs_soluble_accuracy",
        "single_state_vs_multistate_accuracy",
        "perturbation_false_promotion_rate",
        "protein_folding_solved_candidate",
    ]
    for key in required:
        assert key in cert
    assert cert["control_count"] == 20
    assert cert["passed_control_count"] == 20
