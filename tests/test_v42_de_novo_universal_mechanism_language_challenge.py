from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_v42_de_novo_universal_mechanism_challenge_panel_v0.py"
RUNNER = ROOT / "scripts" / "run_v42_de_novo_universal_mechanism_language_challenge_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_builder_creates_24_target_panel() -> None:
    builder = _load(BUILDER, "v42_builder")
    summary = builder.build_panel()
    assert summary["panel_target_count"] == 24
    panel = json.loads((ROOT / "data" / "de_novo_mechanism_language" / "V42" / "panel" / "panel_manifest.json").read_text(encoding="utf-8"))
    assert panel["holdouts_created"] is False
    assert panel["answer_key_available_to_prediction"] is False
    assert panel["panel_groups"]["known_anchor"] == 4
    assert panel["panel_groups"]["membrane_channel_transporter"] == 4


def test_v42_passes_and_beats_required_baselines() -> None:
    builder = _load(BUILDER, "v42_builder_pass")
    runner = _load(RUNNER, "v42_runner_pass")
    builder.build_panel()
    paths = runner.run_v42()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["control_status"] == "V42_DE_NOVO_MECHANISM_LANGUAGE_CHALLENGE_PASSED_CLAIM_DISABLED"
    assert cert["panel_target_count"] == 24
    assert cert["sealed_prediction_count"] == 24
    assert cert["prediction_sealed_before_holdout"] is True
    assert cert["beats_random_baseline"] is True
    assert cert["beats_keyword_baseline"] is True
    assert cert["beats_majority_baseline"] is True
    assert cert["claim_allowed"] is False
    assert cert["folding_problem_solved"] is False


def test_sealed_predictions_have_hash_and_holdouts_reference_hash() -> None:
    builder = _load(BUILDER, "v42_builder_seal")
    runner = _load(RUNNER, "v42_runner_seal")
    builder.build_panel()
    runner.run_v42()
    pred_root = ROOT / "data" / "de_novo_mechanism_language" / "V42" / "predictions_sealed"
    holdout_root = ROOT / "data" / "de_novo_mechanism_language" / "V42" / "holdouts"
    packet = json.loads((pred_root / "TARGET_001" / "sealed_prediction_packet.json").read_text(encoding="utf-8"))
    holdout = json.loads((holdout_root / "TARGET_001" / "post_seal_holdout_manifest.json").read_text(encoding="utf-8"))
    assert packet["prediction_hash"]
    assert packet["sealed_before_holdout"] is True
    assert holdout["holdout_opened_after_prediction_hash"] == packet["prediction_hash"]


def test_leakage_and_shortcut_controls_block_prediction_inputs() -> None:
    builder = _load(BUILDER, "v42_builder_controls")
    runner = _load(RUNNER, "v42_runner_controls")
    builder.build_panel()
    panel = runner.load_panel()
    target = panel["targets"][0]
    bad_coord = {**target, "prediction_sources": [{"source_id": "bad", "source_type": "PDB coordinates", "source_url_or_citation": "PDB coordinate", "coordinate_derived": True}]}
    assert runner.predict_target(bad_coord)["coordinate_derived_source_count_for_prediction"] > 0
    bad_internal = {**target, "prediction_sources": [{"source_id": "bad", "source_type": "runtime", "source_url_or_citation": "first_contact_clean_pharmacotopology_layer_run/report", "internal_runtime_source": True}]}
    assert runner.predict_target(bad_internal)["internal_runtime_source_count_for_prediction"] > 0
    name_only = {**target, "prediction_sources": []}
    assert runner.predict_target(name_only)["mechanism_class"] == "insufficient_evidence_clean_abstain"


def test_generic_hard_class_shortcuts_do_not_promote() -> None:
    builder = _load(BUILDER, "v42_builder_shortcuts")
    runner = _load(RUNNER, "v42_runner_shortcuts")
    builder.build_panel()
    panel = runner.load_panel()
    target = panel["targets"][0]
    membrane_only = {**target, "prediction_sources": [{"source_id": "membrane_only", "feature_tags": ["membrane"], "region_hints": []}]}
    assert runner.predict_target(membrane_only)["mechanism_class"] != "membrane_pore_filter_or_transport"
    multi_no_states = {**target, "prediction_sources": [{"source_id": "multi", "feature_tags": ["metamorphic"], "region_hints": []}]}
    assert runner.predict_target(multi_no_states)["mechanism_class"] != "metamorphic_or_multistate_switch"
    idp = next(t for t in panel["targets"] if t["panel_group"] == "disordered_or_partially_disordered")
    assert runner.predict_target(idp)["mechanism_class"] != "compact_single_fold_core_closure"


def test_contradicting_holdout_fails_not_rescued() -> None:
    builder = _load(BUILDER, "v42_builder_contradict")
    runner = _load(RUNNER, "v42_runner_contradict")
    builder.build_panel()
    panel = runner.load_panel()
    sealed = runner.predict_target(panel["targets"][0])
    holdout = runner._holdout_for_target(panel["targets"][0], sealed)
    holdout["contradicts_prediction"] = True
    result = runner._score_target(sealed, holdout)
    assert result["validation_level"] == "failed_contradicted_holdout"


def test_certificate_has_required_fields_and_controls() -> None:
    builder = _load(BUILDER, "v42_builder_fields")
    runner = _load(RUNNER, "v42_runner_fields")
    builder.build_panel()
    paths = runner.run_v42()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    required = [
        "mechanism_class_accuracy",
        "hard_class_precision",
        "hard_class_recall",
        "operator_region_support_rate",
        "perturbation_prediction_support_rate",
        "baseline_scores",
        "coordinate_derived_source_count_for_prediction",
        "internal_runtime_source_count_for_prediction",
        "holdout_leakage_detected",
    ]
    for key in required:
        assert key in cert
    assert cert["control_count"] == 15
    assert cert["passed_control_count"] == 15
