from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_v40_mechanism_perturbation_sources_v0.py"
RUNNER = ROOT / "scripts" / "run_v40_mechanism_perturbation_pressure_tests_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_builder_writes_v40_sources_and_packets() -> None:
    builder = _load(BUILDER, "v40_builder")
    summary = builder.build_all()
    assert summary["target_count"] == 3
    assert summary["perturbation_packet_count"] == 3
    assert summary["perturbation_count_by_target"]["KcsA"] >= 3
    assert summary["perturbation_count_by_target"]["XCL1_lymphotactin"] >= 3
    assert summary["perturbation_count_by_target"]["alpha_synuclein_SNCA"] >= 3
    assert summary["answer_key_used_for_prediction"] is False
    assert summary["folding_problem_solved"] is False


def test_v40_baseline_passes_after_builder() -> None:
    builder = _load(BUILDER, "v40_builder_pass")
    runner = _load(RUNNER, "v40_runner_pass")
    builder.build_all()
    cert = runner.build_v40()
    assert cert["control_status"] == "V40_MECHANISM_PERTURBATION_PRESSURE_PASSED_CLAIM_DISABLED"
    assert cert["supported_target_count"] == 3
    assert cert["coordinate_derived_source_count"] == 0
    assert cert["internal_runtime_source_count"] == 0
    assert cert["answer_key_used_for_prediction"] is False
    assert cert["claim_allowed"] is False
    assert cert["folding_problem_solved"] is False


def test_required_controls_all_pass_and_record_ablation_levels() -> None:
    builder = _load(BUILDER, "v40_builder_controls")
    runner = _load(RUNNER, "v40_runner_controls")
    builder.build_all()
    controls = runner.run_controls(runner.load_inputs())
    assert len(controls) == 14
    assert [row["passed"] for row in controls] == [True] * 14
    by_id = {row["control_id"]: row for row in controls}
    assert by_id["kcsa_filter_ion_selectivity_removed_weakens_or_invalidates"]["observed"]["validation_level"] in {
        "perturbation_partially_supported_clean_abstain",
        "perturbation_contradicted_or_failed",
    }
    assert by_id["xcl1_state_a_removed_becomes_partial_or_invalid"]["observed"]["validation_level"] == "perturbation_partially_supported_clean_abstain"
    assert by_id["xcl1_state_b_removed_becomes_partial_or_invalid"]["observed"]["validation_level"] == "perturbation_partially_supported_clean_abstain"
    assert by_id["snca_disorder_evidence_removed_becomes_partial_or_invalid"]["observed"]["validation_level"] == "perturbation_partially_supported_clean_abstain"
    assert by_id["placeholder_citations_source_rows_are_blocked"]["observed"]["validation_level"] == "blocked_for_leakage"


def test_required_perturbation_bucket_ablations_are_not_supported() -> None:
    builder = _load(BUILDER, "v40_builder_ablate")
    runner = _load(RUNNER, "v40_runner_ablate")
    builder.build_all()
    inputs = runner.load_inputs()

    packet, manifest = inputs["KcsA"]
    packet = json.loads(json.dumps(packet))
    packet["perturbation_predictions"] = [
        row for row in packet["perturbation_predictions"]
        if row["perturbation_id"] not in {"KCSA_V40_P1_FILTER_SIGNATURE", "KCSA_V40_P2_ION_SELECTIVITY"}
    ]
    kcsa = runner.validate_target(packet, manifest)
    assert kcsa["validation_level"] in {
        "perturbation_partially_supported_clean_abstain",
        "perturbation_contradicted_or_failed",
    }
    assert "filter_signature_perturbation" in kcsa["missing_required_perturbation_buckets"]
    assert "ion_selectivity_perturbation" in kcsa["missing_required_perturbation_buckets"]

    packet, manifest = inputs["XCL1_lymphotactin"]
    packet = json.loads(json.dumps(packet))
    packet["perturbation_predictions"] = [
        row for row in packet["perturbation_predictions"] if row["perturbation_id"] != "XCL1_V40_P1_STATE_A_LOSS"
    ]
    xcl1 = runner.validate_target(packet, manifest)
    assert xcl1["validation_level"] == "perturbation_partially_supported_clean_abstain"
    assert "state_A_loss_or_weakening" in xcl1["missing_required_perturbation_buckets"]

    packet, manifest = inputs["alpha_synuclein_SNCA"]
    packet = json.loads(json.dumps(packet))
    packet["perturbation_predictions"] = [
        row for row in packet["perturbation_predictions"] if row["perturbation_id"] != "SNCA_V40_P1_DISORDER_LOSS"
    ]
    snca = runner.validate_target(packet, manifest)
    assert snca["validation_level"] == "perturbation_partially_supported_clean_abstain"
    assert "disorder_evidence_loss" in snca["missing_required_perturbation_buckets"]


def test_leakage_answer_key_and_placeholder_are_detected() -> None:
    builder = _load(BUILDER, "v40_builder_leak")
    runner = _load(RUNNER, "v40_runner_leak")
    builder.build_all()
    inputs = runner.load_inputs()
    packet, manifest = inputs["KcsA"]

    leaked_manifest = json.loads(json.dumps(manifest))
    leaked_manifest["validation_holdout_sources"].append({
        "source_id": "bad_coordinate",
        "source_name": "1BL8 coordinate contact",
        "source_url_or_citation": "pdb coordinate evidence",
        "source_date_or_version": "accessed_2026-06-11",
        "coordinate_derived": True,
        "independent_from_prediction_sources": True,
    })
    assert runner.validate_target(packet, leaked_manifest)["coordinate_derived_source_count"] > 0

    placeholder_manifest = json.loads(json.dumps(manifest))
    placeholder_manifest["validation_holdout_sources"].append({
        "source_id": "placeholder",
        "source_name": "TODO placeholder",
        "source_url_or_citation": "citation needed",
        "source_date_or_version": "accessed_2026-06-11",
        "independent_from_prediction_sources": True,
    })
    assert runner.validate_target(packet, placeholder_manifest)["placeholder_source_count"] > 0

    leaked_packet = json.loads(json.dumps(packet))
    leaked_packet["answer_key_used_for_prediction"] = True
    assert runner.validate_target(leaked_packet, manifest)["answer_key_source_count"] > 0


def test_swapped_perturbation_holdouts_fail_validation() -> None:
    builder = _load(BUILDER, "v40_builder_swap")
    runner = _load(RUNNER, "v40_runner_swap")
    builder.build_all()
    inputs = runner.load_inputs()
    result = runner.validate_target(inputs["KcsA"][0], inputs["alpha_synuclein_SNCA"][1])
    assert result["validation_level"] != "perturbation_supported_by_holdout"
    assert any("manifest_target_mismatch" in check for check in result["failed_checks"])


def test_writer_outputs_certificate_and_target_results(tmp_path: Path) -> None:
    builder = _load(BUILDER, "v40_builder_writer")
    runner = _load(RUNNER, "v40_runner_writer")
    builder.build_all()
    paths = runner.write_outputs(tmp_path / "out")
    for path in paths.values():
        assert path.exists()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["next_decision"]["claim_allowed"] is False
    assert cert["control_count"] == 14
    for target in ["KcsA", "XCL1_lymphotactin", "alpha_synuclein_SNCA"]:
        result_path = ROOT / "data" / "mechanism_perturbations" / "V40" / "validation" / target / "perturbation_validation_result.json"
        answer_path = ROOT / "data" / "mechanism_perturbations" / "V40" / "validation" / target / "scientist_question_answer.md"
        assert result_path.exists()
        assert answer_path.exists()
