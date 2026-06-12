from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "run_v57_blind_protein_esperanto_generalization_gate_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_v57_passes_with_frozen_engine_and_writes_certificate() -> None:
    runner = _load(RUNNER, "v57_runner_pass")
    paths = runner.run_v57()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["status"] == "V57_GENERALIZATION_BLOCKED_ENGINE_NEEDS_REVISION"
    assert cert["target_count"] == 5
    assert cert["fresh_regime_count"] >= 4
    assert cert["supported_validation_count"] == 1
    assert cert["abstained_targets"] == ["V57_GB1_DOMAIN"]
    assert cert["passed_control_count"] < cert["control_count"]
    assert cert["engine_modified_for_v57"] is False
    assert cert["folding_problem_solved"] is False
    assert cert["atomistic_md_performed"] is False
    assert cert["readme_touched"] is False
    assert paths["report"].exists()


def test_v57_frozen_engine_declaration_has_no_new_language() -> None:
    runner = _load(RUNNER, "v57_runner_frozen")
    paths = runner.run_v57()
    declaration = json.loads(paths["engine_declaration"].read_text(encoding="utf-8"))
    assert declaration["engine_modified_for_v57"] is False
    assert "new operators" in declaration["disallowed_v57_changes"]
    assert "new mechanism classes" in declaration["disallowed_v57_changes"]
    assert "README updates" in declaration["disallowed_v57_changes"]
    assert declaration["frozen_operator_names"] == runner.UNIVERSAL_OPERATORS
    assert declaration["frozen_mechanism_classes"] == runner.MECHANISM_CLASSES


def test_v57_blind_target_manifest_uses_fresh_regime_slots() -> None:
    runner = _load(RUNNER, "v57_runner_targets")
    paths = runner.run_v57()
    manifest = json.loads(paths["target_manifest"].read_text(encoding="utf-8"))
    slots = {row["regime_slot"] for row in manifest["targets"]}
    assert manifest["answer_key_available_to_prediction"] is False
    assert manifest["holdouts_available_before_seal"] is False
    assert manifest["readme_update_allowed"] is False
    assert slots == {
        "globular_closure",
        "idp_phase_or_disorder",
        "membrane_proteostasis",
        "metamorphic_switch",
        "short_interface_hijacking",
    }
    for forbidden in ["V44", "V45", "V46", "V47", "V48"]:
        assert all(forbidden not in row["target_id"] for row in manifest["targets"])


def test_v57_sealed_packets_select_expected_fresh_grammars() -> None:
    runner = _load(RUNNER, "v57_runner_grammars")
    runner.run_v57()
    packet_root = ROOT / "data" / "protein_esperanto_engine" / "V57" / "sealed_simulation_packets"
    expected = {
            "V57_GB1_DOMAIN": "insufficient_evidence_clean_abstain",
        "V57_HNRNPA1_LCD": "intrinsic_disorder_phase_separation",
        "V57_RHODOPSIN_P23H": "membrane_multidomain_folding_proteostasis",
        "V57_XCL1_SWITCH": "metamorphic_fold_switching",
        "V57_HIV_TAT": "short_region_host_interface_hijacking",
    }
    for target_id, mechanism_class in expected.items():
        packet = json.loads((packet_root / target_id / "sealed_simulation_packet.json").read_text(encoding="utf-8"))
        assert packet["sealed_before_holdout"] is True
        assert packet["selected_mechanism_grammar"]["mechanism_class"] == mechanism_class
        assert packet["evidence_manifest"]["coordinate_derived_source_count_before_prediction"] == 0
        assert packet["folding_problem_solved"] is False


def test_v57_wrong_grammar_controls_fail_or_abstain_for_every_target() -> None:
    runner = _load(RUNNER, "v57_runner_wrong")
    runner.run_v57()
    wrong_root = ROOT / "data" / "protein_esperanto_engine" / "V57" / "wrong_grammar_controls"
    for target_id in runner.FRESH_TARGETS:
        packet = json.loads((wrong_root / target_id / "wrong_grammar_packet.json").read_text(encoding="utf-8"))
        assert packet["selected_mechanism_grammar"]["forced_grammar_rejected"] is True
        assert packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain"


def test_v57_perturbation_response_separates_correct_from_controls() -> None:
    runner = _load(RUNNER, "v57_runner_perturb")
    paths = runner.run_v57()
    table = json.loads(paths["perturbation_table"].read_text(encoding="utf-8"))
    rows = [row for target in table["targets"] for row in target["rows"]]
    assert rows
    assert any(not row["direction_passed"] for row in rows)
    assert any(row["perturbation_id"] == "GB1_CORE_HYDROPHOBIC_DAMAGE" and row["observed_direction"] == "unchanged" for row in rows)
    assert any(row["perturbation_id"] == "RHO_CHAPERONE_RESCUE_CONTEXT" and row["observed_direction"] == "increase" for row in rows)
    assert any(row["perturbation_id"] == "TAT_BASIC_REGION_DAMAGE" and row["observed_direction"] == "decrease" for row in rows)
    assert all(row["observed_direction"] == "unchanged" for row in rows if row["perturbation_id"].endswith("CONTROL"))


def test_v57_holdouts_open_after_hash_and_validate() -> None:
    runner = _load(RUNNER, "v57_runner_holdouts")
    paths = runner.run_v57()
    validation_report = json.loads(paths["validation_report"].read_text(encoding="utf-8"))
    assert validation_report["supported_validation_count"] == 1
    for validation in validation_report["validations"]:
        assert validation["holdout_opened_after_prediction_hash"]
    supported = [row for row in validation_report["validations"] if row["score_label"] == "supported"]
    assert [row["target_id"] for row in supported] == ["V57_HNRNPA1_LCD"]


def test_v57_controls_cover_leakage_random_sequence_and_holdout_order() -> None:
    runner = _load(RUNNER, "v57_runner_controls")
    cert = json.loads(runner.run_v57()["certificate"].read_text(encoding="utf-8"))
    controls = {row["control_id"]: row for row in cert["controls"]}
    for control_id in [
        "v57_engine_source_frozen",
        "v57_no_new_operators_or_mechanism_classes",
        "v57_all_wrong_grammars_fail",
        "v57_random_sequence_abstains",
        "v57_coordinate_leakage_blocks",
        "v57_holdout_before_seal_blocks",
        "v57_all_holdouts_reference_hash",
        "v57_folding_problem_solved_never_true",
    ]:
        assert controls[control_id]["passed"] is True
    assert controls["V57_GB1_DOMAIN_selected_expected_grammar"]["passed"] is False
