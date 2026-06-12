from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
V59_RUNNER = ROOT / "scripts" / "run_v59_real_folding_process_replication_panel_v0.py"
V60_RUNNER = ROOT / "scripts" / "run_v60_bounded_protein_esperanto_solved_claim_gate_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def v59_outputs() -> dict[str, object]:
    runner = _load(V59_RUNNER, "v59_runner_for_tests")
    paths = runner.run_v59()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    scoring = json.loads(paths["scoring_report"].read_text(encoding="utf-8"))
    return {"paths": paths, "certificate": cert, "scoring": scoring}


def test_v59_process_panel_passes_after_explicit_engine_hardening(v59_outputs: dict[str, object]) -> None:
    cert = v59_outputs["certificate"]
    assert isinstance(cert, dict)
    assert cert["status"] == "V59_PARTIAL_PROCESS_REPLICATION_WITH_ABSTENTION_PASSED_REVIEW_REQUIRED"
    assert cert["target_count"] == 10
    assert cert["accepted_count"] == 9
    assert cert["accepted_with_caution_count"] == 3
    assert cert["abstain_recommended_count"] == 1
    assert cert["accepted_process_accuracy"] == 1.0
    assert cert["raw_process_accuracy"] == 0.9
    assert cert["passed_control_count"] == cert["control_count"] == 14
    assert cert["engine_biology_modified"] is True
    assert cert["folding_problem_solved"] is False
    assert cert["universal_folding_solved"] is False
    assert cert["atomistic_md_executed"] is False
    assert cert["coordinate_truth_used_before_seal"] is False
    assert cert["alphafold_used_before_seal"] is False


def test_v59_engine_revision_is_bounded_not_operator_expansion(v59_outputs: dict[str, object]) -> None:
    cert = v59_outputs["certificate"]
    controls = {row["control_id"]: row for row in cert["controls"]}
    revision = controls["engine_revision_explicit_and_bounded"]["observed"]
    assert revision["engine_revision_required_by_v59"] is True
    assert revision["engine_operator_set_modified"] is False
    assert revision["engine_mechanism_class_set_modified"] is False
    assert revision["frozen_after_v59_revision"] is True
    assert "hardened mechanism grammar selection" in " ".join(revision["engine_revision_scope"])


def test_v59_scoring_covers_process_classes_and_all_rows_are_supported(v59_outputs: dict[str, object]) -> None:
    rows = v59_outputs["scoring"]["rows"]
    assert len(rows) == 10
    assert {row["expected_process_class"] for row in rows} >= {"two_state", "intermediate_bearing", "multi_basin"}
    supported = [row for row in rows if row["score_label"] == "supported"]
    abstained = [row for row in rows if row["score_label"] == "abstained"]
    assert len(supported) == 9
    assert [row["target_id"] for row in abstained] == ["V59_VILLIN_HP35"]
    assert all(row["p1_process_class"] and row["p2_early_late_ordering"] and row["p3_mutation_sensitivity"] for row in supported)
    assert all(row["sealed_prediction_hash"] == row["holdout_opened_after_prediction_hash"] for row in rows)


def test_v59_engine_hardening_routes_previous_misses_to_globular_process_packets(v59_outputs: dict[str, object]) -> None:
    packet_root = ROOT / "data" / "protein_esperanto_engine" / "V59" / "sealed_predictions"
    for target_id in ["V59_CI2", "V59_UBIQUITIN", "V59_TRP_CAGE"]:
        engine_packet = json.loads((packet_root / target_id / "engine_packet.json").read_text(encoding="utf-8"))
        process_packet = json.loads((packet_root / target_id / "sealed_process_prediction_packet.json").read_text(encoding="utf-8"))
        assert engine_packet["selected_mechanism_grammar"]["mechanism_class"] == "globular_closure"
        assert process_packet["process_prediction"]["process_decision"] == "accepted"


def test_v60_computes_bounded_claim_not_universal_solved(v59_outputs: dict[str, object]) -> None:
    runner = _load(V60_RUNNER, "v60_runner_for_tests")
    paths = runner.run_v60()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["status"] == "V60_SOLVED_CLAIM_BLOCKED_ENGINE_REVISION_REQUIRED"
    assert cert["protein_esperanto_engine_claim_allowed"] is False
    assert cert["bounded_accepted_target_process_replication_supported"] is False
    assert cert["bounded_real_sequence_behavior_supported"] is False
    assert cert["universal_folding_problem_solved"] is False
    assert cert["atomistic_folding_solved"] is False
    assert cert["all_sequence_prediction_solved"] is False
    assert cert["folding_problem_solved"] is False
    assert cert["alphafold_replaced"] is False
    assert cert["engine_revision_after_v58_required"] is True
    assert cert["engine_operator_set_modified_in_v59"] is False
    assert cert["engine_mechanism_class_set_modified_in_v59"] is False
    assert "not yet claim-ready" in cert["allowed_claim"]
    assert any("Universal protein folding is solved" in claim for claim in cert["forbidden_claims"])
