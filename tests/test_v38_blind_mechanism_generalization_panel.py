from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_v38_blind_mechanism_generalization_panel_v0.py"
RUNNER = ROOT / "scripts" / "run_v38_blind_mechanism_generalization_panel_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_builder_creates_masked_panel_and_separate_answer_key(tmp_path: Path) -> None:
    builder = _load(BUILDER, "v38_builder")
    summary = builder.write_panel(tmp_path / "panel")
    assert summary["panel_target_count"] == 9
    assert summary["known_positive_count"] == 3
    assert summary["decoy_count"] == 6
    assert summary["answer_key_used_for_assignment"] is False
    assert len(list((tmp_path / "panel" / "masked_inputs").glob("TARGET_*.json"))) == 9
    assert (tmp_path / "panel" / "answer_key.json").exists()


def test_blind_panel_passes_generalization(tmp_path: Path) -> None:
    builder = _load(BUILDER, "v38_builder_pass")
    runner = _load(RUNNER, "v38_runner_pass")
    panel_root = tmp_path / "panel"
    builder.write_panel(panel_root)
    cert = runner.build_v38(panel_root)
    assert cert["control_status"] == "V38_BLIND_MECHANISM_GENERALIZATION_PASSED_CLAIM_DISABLED"
    assert cert["correct_known_positive_count"] == 3
    assert cert["known_positive_count"] == 3
    assert cert["false_positive_hard_grammar_count"] == 0
    assert cert["false_negative_known_positive_count"] == 0
    assert cert["masked_panel_used"] is True
    assert cert["answer_key_used_for_assignment"] is False
    assert cert["coordinate_derived_source_count"] == 0
    assert cert["internal_runtime_source_count"] == 0
    assert cert["claim_allowed"] is False
    assert cert["folding_problem_solved"] is False


def test_required_controls_all_pass(tmp_path: Path) -> None:
    builder = _load(BUILDER, "v38_builder_controls")
    runner = _load(RUNNER, "v38_runner_controls")
    panel_root = tmp_path / "panel"
    builder.write_panel(panel_root)
    masked = runner.load_masked_panel(panel_root)
    unmasked = runner.load_unmasked_panel(panel_root)
    key = json.loads((panel_root / "answer_key.json").read_text(encoding="utf-8"))
    controls = runner.run_controls(masked, key, unmasked)
    assert len(controls) == 12
    assert [row["passed"] for row in controls] == [True] * 12


def test_decoys_do_not_promote_to_hard_grammar(tmp_path: Path) -> None:
    builder = _load(BUILDER, "v38_builder_decoys")
    runner = _load(RUNNER, "v38_runner_decoys")
    panel_root = tmp_path / "panel"
    builder.write_panel(panel_root)
    masked = runner.load_masked_panel(panel_root)
    assignments = {row["masked_id"]: runner.assign_masked_dossier(row)["assigned_mechanism_class"] for row in masked}
    assert assignments["TARGET_004"] == "other_membrane_or_transport_context"
    assert assignments["TARGET_005"] == "other_membrane_or_transport_context"
    assert assignments["TARGET_006"] == "soluble_single_or_contextual_fold_not_metamorphic"
    assert assignments["TARGET_008"] == "soluble_single_or_contextual_fold_not_metamorphic"


def test_generic_coordinate_internal_placeholder_and_key_leakage_block() -> None:
    runner = _load(RUNNER, "v38_runner_blocks")
    generic = {"masked_id": "TARGET_X", "evidence_buckets": ["generic"], "mechanism_evidence": ["generic protein"], "grammar_rules": {}, "source_provenance": []}
    assert runner.assign_masked_dossier(generic)["assigned_mechanism_class"] == "insufficient_evidence_clean_abstain"
    coord = {"masked_id": "TARGET_X", "evidence_buckets": [], "mechanism_evidence": [], "grammar_rules": {}, "source_provenance": [{"coordinate_derived": True, "source_name": "1BL8"}]}
    assert runner.assign_masked_dossier(coord)["assignment_status"] == "blocked_leakage"
    internal = {"masked_id": "TARGET_X", "evidence_buckets": [], "mechanism_evidence": [], "grammar_rules": {}, "source_provenance": [{"source_url_or_citation": "first_contact_clean_pharmacotopology_layer_run/report.json"}]}
    assert runner.assign_masked_dossier(internal)["assignment_status"] == "blocked_leakage"
    placeholder = {"masked_id": "TARGET_X", "evidence_buckets": [], "mechanism_evidence": [], "grammar_rules": {}, "source_provenance": [{"source_name": "TODO placeholder"}]}
    assert runner.assign_masked_dossier(placeholder)["assignment_status"] == "blocked_key"
    key = {"masked_id": "TARGET_X", "expected_mechanism_class": "x", "evidence_buckets": [], "mechanism_evidence": [], "grammar_rules": {}, "source_provenance": []}
    assert runner.assign_masked_dossier(key)["assignment_status"] == "blocked_key"
    name = {"masked_id": "TARGET_X", "evidence_buckets": [], "mechanism_evidence": ["KcsA"], "grammar_rules": {}, "source_provenance": []}
    assert runner.assign_masked_dossier(name)["assignment_status"] == "blocked_key"


def test_writer_outputs_artifacts(tmp_path: Path) -> None:
    builder = _load(BUILDER, "v38_builder_writer")
    runner = _load(RUNNER, "v38_runner_writer")
    panel_root = tmp_path / "panel"
    builder.write_panel(panel_root)
    cert = runner.build_v38(panel_root)
    paths = runner.write_outputs(tmp_path / "out", panel_root, cert)
    for path in paths.values():
        assert path.exists()
    written = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert written["next_decision"]["claim_allowed"] is False
