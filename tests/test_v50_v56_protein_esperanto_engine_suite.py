from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "run_v50_v56_protein_esperanto_engine_suite_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_v50_v56_suite_passes_and_writes_certificate() -> None:
    runner = _load(RUNNER, "v50_v56_runner_pass")
    paths = runner.run_v50_v56()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["status"] == "V50_TO_V56_PROTEIN_ESPERANTO_ENGINE_PASSED_REVIEW_REQUIRED"
    assert cert["v50_grammar_extraction_passed"] is True
    assert cert["v51_engine_contract_frozen"] is True
    assert cert["v52_coarse_simulator_mvp_passed"] is True
    assert cert["v53_hard_class_battery_passed"] is True
    assert cert["v54_perturbation_engine_passed"] is True
    assert cert["v55_postseal_validation_passed"] is True
    assert cert["v56_openmm_bridge_spec_passed"] is True
    assert cert["passed_control_count"] == cert["control_count"]
    assert cert["folding_problem_solved"] is False
    assert cert["atomistic_md_executed"] is False
    assert paths["report"].exists()


def test_v50_grammar_rewrites_every_required_hard_class() -> None:
    runner = _load(RUNNER, "v50_v56_runner_grammar")
    grammar = runner.extract_v50_grammar()
    sentences = {row["mechanism_class"]: row["grammar_sentence"] for row in grammar["transition_rules"]}
    for mechanism_class in [
        "intrinsic_disorder_phase_separation",
        "membrane_multidomain_folding_proteostasis",
        "metamorphic_fold_switching",
        "short_region_host_interface_hijacking",
    ]:
        assert mechanism_class in sentences
        assert " -> " in sentences[mechanism_class]
        assert len(sentences[mechanism_class].split(" -> ")) == 3
    for section in [
        "universal_marks",
        "universal_operators",
        "mechanism_classes",
        "state_variables",
        "environmental_pressures",
        "allowed_evidence_classes",
        "forbidden_evidence_classes",
        "transition_rules",
        "perturbation_prediction_rules",
        "falsification_rules",
        "null_controls",
        "simulation_readiness_decision",
    ]:
        assert section in grammar["required_sections"]


def test_evidence_boundary_blocks_leakage_and_requires_spatial_proxy_tags() -> None:
    runner = _load(RUNNER, "v50_v56_runner_gate")
    gate = runner.evidence_boundary_gate([
        {
            "source_id": "SEQ_ONLY",
            "source_role": "prediction_input",
            "source_class": "pure_non_coordinate",
            "coordinate_derived": False,
            "internal_runtime_source": False,
        },
        {
            "source_id": "BAD_COORD",
            "source_role": "prediction_input",
            "source_class": "coordinate_derived",
            "coordinate_derived": True,
        },
        {
            "source_id": "BAD_RUNTIME",
            "source_role": "prediction_input",
            "source_class": "internal_runtime",
            "internal_runtime": True,
        },
        {
            "source_id": "UNTAGGED_FRET",
            "source_role": "prediction_input",
            "source_type": "FRET distance summary",
            "coordinate_derived": False,
        },
        {
            "source_id": "TAGGED_DCA",
            "source_role": "prediction_input",
            "source_class": "spatial_proxy_non_coordinate",
            "spatial_proxy": True,
        },
    ])
    assert gate["allowed_initialization_source_ids"] == ["SEQ_ONLY", "TAGGED_DCA"]
    assert gate["coordinate_derived_source_count_before_prediction"] == 1
    assert gate["internal_runtime_source_count_for_prediction"] == 1
    assert gate["spatial_proxy_untagged_or_misused_count"] == 1


def test_sequence_field_contract_contains_v51_state_vector_marks() -> None:
    runner = _load(RUNNER, "v50_v56_runner_sequence_field")
    profile = runner.load_target_profile("V48_SARS2_ORF6")
    packet = runner.build_sealed_simulation_packet(
        target_id=profile["target_id"],
        target_name=profile["target_name"],
        sequence=profile["sequence"],
        sources=profile["sources"],
        focus_regions=profile["focus_regions"],
        perturbations=profile["perturbations"],
    )
    residue = packet["initial_sequence_field_map"]["residues"][0]
    for key in [
        "position_index",
        "residue_identity",
        "segment_id",
        "charge_mark",
        "hydrophobic_mark",
        "disorder_mark",
        "secondary_propensity_mark",
        "interface_mark",
        "membrane_mark",
        "operator_activations",
        "current_state",
        "state_confidence",
    ]:
        assert key in residue
    assert packet["initial_sequence_field_map"]["coordinate_truth_used"] is False


def test_hard_class_battery_outputs_different_regime_trajectories() -> None:
    runner = _load(RUNNER, "v50_v56_runner_battery")
    paths = runner.run_v50_v56()
    battery = json.loads(paths["battery"].read_text(encoding="utf-8"))
    rows = {row["target_id"]: row for row in battery["rows"]}
    assert rows["V44_FUS_LC"]["mechanism_class"] == "intrinsic_disorder_phase_separation"
    assert rows["V45_TDP43_LCD"]["mechanism_class"] == "intrinsic_disorder_phase_separation"
    assert rows["V46_CFTR_F508DEL"]["mechanism_class"] == "membrane_multidomain_folding_proteostasis"
    assert rows["V47_RFAH_CTD"]["mechanism_class"] == "metamorphic_fold_switching"
    assert rows["V48_SARS2_ORF6"]["mechanism_class"] == "short_region_host_interface_hijacking"
    assert rows["V44_FUS_LC"]["compact_single_fold_probability"] <= 0.12
    assert rows["V47_RFAH_CTD"]["final_state_basin_occupancy"]["alpha_context_basin"] >= 0.30
    assert rows["V47_RFAH_CTD"]["final_state_basin_occupancy"]["beta_released_basin"] >= 0.30
    assert rows["V48_SARS2_ORF6"]["final_state_basin_occupancy"]["host_interface_engaged"] >= 0.60


def test_controls_cover_random_shuffled_wrong_grammar_and_holdout_order() -> None:
    runner = _load(RUNNER, "v50_v56_runner_controls")
    cert = json.loads(runner.run_v50_v56()["certificate"].read_text(encoding="utf-8"))
    controls = {row["control_id"]: row for row in cert["controls"]}
    for control_id in [
        "random_sequence_control",
        "shuffled_sequence_control",
        "swapped_evidence_control",
        "wrong_target_control",
        "generic_annotation_only_control",
        "coordinate_leakage_control",
        "internal_runtime_leakage_control",
        "forced_wrong_grammar_control",
        "failed_prediction_not_repaired_after_holdout",
        "holdout_opened_before_seal_control",
        "wild_type_mutant_direction_control",
    ]:
        assert controls[control_id]["passed"] is True


def test_v54_perturbations_move_in_expected_direction() -> None:
    runner = _load(RUNNER, "v50_v56_runner_perturb")
    paths = runner.run_v50_v56()
    data = json.loads(paths["perturbations"].read_text(encoding="utf-8"))
    assert data["passed"] is True
    rows = [row for target in data["targets"] for row in target["rows"]]
    assert rows
    assert all(row["direction_passed"] for row in rows)
    assert any(row["perturbation_id"] == "CFTR_F508DEL_DAMAGE" and row["observed_direction"] == "decrease" for row in rows)
    assert any(row["perturbation_id"] == "ORF6_C_TERMINAL_DISRUPTION" and row["observed_direction"] == "decrease" for row in rows)
