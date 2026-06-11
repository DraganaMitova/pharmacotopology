from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_v39_mechanism_prediction_holdouts_v0.py"
RUNNER = ROOT / "scripts" / "run_v39_mechanism_to_falsifiable_prediction_validation_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_builder_writes_prediction_and_holdout_files() -> None:
    builder = _load(BUILDER, "v39_builder")
    summary = builder.build_all()
    assert summary["target_count"] == 3
    assert summary["prediction_packet_count"] == 3
    assert summary["holdout_source_count"] >= 12
    assert summary["answer_key_used_for_prediction"] is False


def test_v39_baseline_passes_after_builder() -> None:
    builder = _load(BUILDER, "v39_builder_pass")
    runner = _load(RUNNER, "v39_runner_pass")
    builder.build_all()
    cert = runner.build_v39()
    assert cert["control_status"] == "V39_MECHANISM_PREDICTIONS_VALIDATED_CLAIM_DISABLED"
    assert cert["validated_target_count"] == 3
    assert cert["coordinate_derived_source_count"] == 0
    assert cert["internal_runtime_source_count"] == 0
    assert cert["answer_key_used_for_prediction"] is False
    assert cert["claim_allowed"] is False
    assert cert["folding_problem_solved"] is False


def test_required_controls_all_pass() -> None:
    builder = _load(BUILDER, "v39_builder_controls")
    runner = _load(RUNNER, "v39_runner_controls")
    builder.build_all()
    controls = runner.run_controls(runner.load_inputs())
    assert len(controls) == 12
    assert [row["passed"] for row in controls] == [True] * 12
    by_id = {row["control_id"]: row for row in controls}
    assert by_id["kcsa_filter_ion_removed_partial_or_invalid"]["observed"]["validation_level"] == "partially_supported_clean_abstain"
    assert by_id["xcl1_one_state_only_partial_or_invalid"]["observed"]["validation_level"] == "partially_supported_clean_abstain"
    assert by_id["snca_disorder_removed_partial_or_invalid"]["observed"]["validation_level"] == "partially_supported_clean_abstain"
    assert by_id["placeholder_citations_blocked"]["observed"]["validation_level"] == "blocked_for_leakage"


def test_required_mechanism_bucket_ablations_are_not_high_confidence() -> None:
    builder = _load(BUILDER, "v39_builder_ablate")
    runner = _load(RUNNER, "v39_runner_ablate")
    builder.build_all()
    inputs = runner.load_inputs()

    packet, holdout = inputs["KcsA"]
    packet = json.loads(json.dumps(packet))
    packet["falsifiable_predictions"] = [
        row for row in packet["falsifiable_predictions"] if row["prediction_id"] != "KCSA_P1_FILTER"
    ]
    kcsa = runner.validate_target(packet, holdout)
    assert kcsa["validation_level"] == "partially_supported_clean_abstain"
    assert "filter_or_signature_holdout" in kcsa["missing_required_prediction_buckets"]

    packet, holdout = inputs["XCL1_lymphotactin"]
    packet = json.loads(json.dumps(packet))
    packet["falsifiable_predictions"] = [
        row for row in packet["falsifiable_predictions"] if row["prediction_id"] != "XCL1_P2_TWO_STATES"
    ]
    xcl1 = runner.validate_target(packet, holdout)
    assert xcl1["validation_level"] == "partially_supported_clean_abstain"
    assert "state_A_state_B_holdout" in xcl1["missing_required_prediction_buckets"]

    packet, holdout = inputs["alpha_synuclein_SNCA"]
    packet = json.loads(json.dumps(packet))
    packet["falsifiable_predictions"] = [
        row for row in packet["falsifiable_predictions"] if row["prediction_id"] != "SNCA_P1_DISORDER"
    ]
    snca = runner.validate_target(packet, holdout)
    assert snca["validation_level"] == "partially_supported_clean_abstain"
    assert "disorder_holdout" in snca["missing_required_prediction_buckets"]


def test_leakage_and_placeholder_are_detected() -> None:
    builder = _load(BUILDER, "v39_builder_leak")
    runner = _load(RUNNER, "v39_runner_leak")
    builder.build_all()
    inputs = runner.load_inputs()
    packet, holdout = inputs["KcsA"]
    holdout = json.loads(json.dumps(holdout))
    holdout["holdout_validation_evidence"].append({"source_id": "bad", "source_name": "1BL8 coordinate contact", "coordinate_derived": True})
    assert runner.validate_target(packet, holdout)["coordinate_derived_source_count"] > 0
    holdout = json.loads(json.dumps(inputs["KcsA"][1]))
    holdout["holdout_validation_evidence"].append({"source_id": "bad", "source_name": "TODO placeholder", "source_url_or_citation": "citation needed"})
    assert runner.validate_target(packet, holdout)["placeholder_source_count"] > 0
    packet = json.loads(json.dumps(inputs["KcsA"][0]))
    packet["answer_key_used_for_prediction"] = True
    assert runner.validate_target(packet, inputs["KcsA"][1])["answer_key_source_count"] > 0


def test_swapped_holdout_fails_validation() -> None:
    builder = _load(BUILDER, "v39_builder_swap")
    runner = _load(RUNNER, "v39_runner_swap")
    builder.build_all()
    inputs = runner.load_inputs()
    result = runner.validate_target(inputs["KcsA"][0], inputs["alpha_synuclein_SNCA"][1])
    assert result["validation_level"] != "supported_by_independent_holdout"


def test_writer_outputs_certificate_and_validation_results(tmp_path: Path) -> None:
    builder = _load(BUILDER, "v39_builder_writer")
    runner = _load(RUNNER, "v39_runner_writer")
    builder.build_all()
    paths = runner.write_outputs(tmp_path / "out")
    for path in paths.values():
        assert path.exists()
    cert = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert cert["next_decision"]["claim_allowed"] is False
