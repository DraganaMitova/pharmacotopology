from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "run_v58_real_sequence_time_blind_folding_replication_gate_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_v58_passes_real_sequence_gate_from_cached_intake() -> None:
    runner = _load(RUNNER, "v58_runner_pass")
    paths = runner.run_v58()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["status"] == "V58_REAL_SEQUENCE_TIME_BLIND_FOLDING_REPLICATION_PASSED_REVIEW_REQUIRED"
    assert cert["target_count"] == 20
    assert cert["target_selection_manual"] is False
    assert cert["level1_regime_selection_supported_count"] > cert["baseline_scores"]["majority_class_accuracy"] * cert["target_count"]
    assert cert["level3_topology_or_observable_supported_count"] == cert["level1_regime_selection_supported_count"]
    assert cert["passed_control_count"] == cert["control_count"]
    assert cert["folding_problem_solved"] is False
    assert cert["atomistic_md_executed"] is False
    assert cert["coordinate_truth_used_before_seal"] is False
    assert cert["alphafold_used_before_seal"] is False
    assert cert["failure_cases_reported"] is True
    assert paths["report"].exists()


def test_v58_raw_intake_is_cached_real_rcsb_recent_release_data() -> None:
    runner = _load(RUNNER, "v58_runner_cache")
    paths = runner.run_v58()
    raw = json.loads(paths["raw_candidate_cache"].read_text(encoding="utf-8"))
    assert raw["kind"] == "V58_RCSB_RECENT_RELEASE_RAW_CANDIDATE_ENTITIES_v0"
    assert raw["candidate_entity_count"] >= 20
    assert "RCSB" in raw["source"]
    assert raw["recent_release_start"] == runner.RECENT_RELEASE_START
    first = raw["candidates"][0]
    assert first["entry_id"]
    assert first["entity_id"]
    assert set(first["sequence"]) <= set("ACDEFGHIKLMNPQRSTVWY")


def test_v58_target_selection_is_automatic_and_has_twenty_targets() -> None:
    runner = _load(RUNNER, "v58_runner_manifest")
    paths = runner.run_v58()
    manifest = json.loads(paths["target_manifest"].read_text(encoding="utf-8"))
    assert manifest["kind"] == "V58_REAL_SEQUENCE_AUTOMATIC_TARGET_SELECTION_v0"
    assert manifest["target_selection_manual"] is False
    assert manifest["target_count_selected"] == 20
    assert "deterministic cyclic fill" in manifest["selection_rule"]
    assert len({row["target_id"] for row in manifest["selected_targets"]}) == 20


def test_v58_sequence_only_and_annotation_scoring_are_both_emitted() -> None:
    runner = _load(RUNNER, "v58_runner_scoring")
    paths = runner.run_v58()
    sequence_only = json.loads(paths["sequence_only_scoring"].read_text(encoding="utf-8"))
    annotation = json.loads(paths["annotation_scoring"].read_text(encoding="utf-8"))
    assert len(sequence_only["rows"]) == 20
    assert len(annotation["rows"]) == 20
    assert sum(1 for row in annotation["rows"] if row["level1_regime_selection"]) >= sum(
        1 for row in sequence_only["rows"] if row["level1_regime_selection"]
    )
    assert all(row["level4_process_replication"] == "not_claimed_v59_required" for row in annotation["rows"])


def test_v58_sealed_packets_do_not_use_coordinate_or_alphafold_inputs() -> None:
    runner = _load(RUNNER, "v58_runner_seals")
    runner.run_v58()
    packet_root = ROOT / "data" / "protein_esperanto_engine" / "V58" / "sealed_predictions" / "sequence_plus_annotation"
    holdout_root = ROOT / "data" / "protein_esperanto_engine" / "V58" / "holdouts_postseal" / "sequence_plus_annotation"
    packets = list(packet_root.glob("*/sealed_simulation_packet.json"))
    assert len(packets) == 20
    for packet_path in packets:
        packet = json.loads(packet_path.read_text(encoding="utf-8"))
        holdout = json.loads((holdout_root / packet["target_id"] / "postseal_holdout_manifest.json").read_text(encoding="utf-8"))
        assert packet["sealed_before_holdout"] is True
        assert packet["evidence_manifest"]["coordinate_derived_source_count_before_prediction"] == 0
        assert packet["evidence_manifest"]["internal_runtime_source_count_for_prediction"] == 0
        assert holdout["holdout_opened_after_prediction_hash"] == packet["prediction_hash"]
        assert packet["folding_problem_solved"] is False


def test_v58_failure_report_is_explicit_not_repaired_away() -> None:
    runner = _load(RUNNER, "v58_runner_failures")
    paths = runner.run_v58()
    report = json.loads(paths["failure_report"].read_text(encoding="utf-8"))
    assert report["failure_cases_reported"] is True
    assert len(report["failure_cases"]) >= 1
    for row in report["failure_cases"]:
        assert row["target_id"]
        assert row["predicted"]
        assert row["expected"]


def test_v58_controls_include_leakage_wrong_grammar_and_baseline_checks() -> None:
    runner = _load(RUNNER, "v58_runner_controls")
    cert = json.loads(runner.run_v58()["certificate"].read_text(encoding="utf-8"))
    controls = {row["control_id"]: row for row in cert["controls"]}
    for control_id in [
        "v58_target_selection_automatic",
        "v58_target_count_20",
        "v58_engine_source_frozen",
        "v58_no_new_language",
        "v58_all_sealed_before_holdout",
        "v58_all_wrong_grammars_fail",
        "v58_coordinate_leakage_blocks",
        "v58_internal_runtime_blocks",
        "v58_random_sequence_abstains",
        "v58_annotation_accuracy_beats_majority_baseline",
        "v58_failure_cases_reported",
        "v58_folding_problem_solved_never_true",
    ]:
        assert controls[control_id]["passed"] is True
    assert cert["baseline_scores"]["annotation_engine_accuracy"] > cert["baseline_scores"]["majority_class_accuracy"]
