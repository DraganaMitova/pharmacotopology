from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_v36_real_evidence_dossiers_v0.py"
GATE = ROOT / "scripts" / "run_v36_real_evidence_dossier_gate_v0.py"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_builder_writes_three_target_dossiers(tmp_path: Path) -> None:
    builder = _load(BUILDER, "v36_builder")
    summary = builder.build_all_dossiers(tmp_path / "external_evidence_dossiers")
    assert summary["target_count"] == 3
    assert summary["source_counts_by_target"] == {
        "KcsA": 4,
        "XCL1_lymphotactin": 4,
        "alpha_synuclein_SNCA": 4,
    }
    for target in summary["targets"]:
        for name in ["source_manifest", "evidence_dossier", "evidence_table", "acquisition_log", "rejected_sources"]:
            assert Path(summary["artifacts"][target][name]).exists()


def test_gate_ready_with_builder_dossiers(tmp_path: Path) -> None:
    builder = _load(BUILDER, "v36_builder_ready")
    gate = _load(GATE, "v36_gate_ready")
    data_root = tmp_path / "external_evidence_dossiers"
    builder.build_all_dossiers(data_root)
    packages = gate.load_packages(data_root)
    cert = gate.build_v36(packages)
    assert cert["control_status"] == "V36_REAL_EVIDENCE_DOSSIERS_READY_CLAIM_DISABLED"
    assert cert["targets_ready"] == ["KcsA", "XCL1_lymphotactin", "alpha_synuclein_SNCA"]
    assert cert["targets_partial"] == []
    assert cert["targets_blocked"] == []
    assert cert["external_source_count"] == 12
    assert cert["coordinate_derived_source_count"] == 0
    assert cert["internal_runtime_source_count"] == 0
    assert cert["placeholder_source_count"] == 0
    assert cert["claim_allowed"] is False
    assert cert["new_MD_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["native_metrics_used_for_selection"] is False
    assert cert["coordinate_truth_used_before_selection"] is False
    assert cert["positive_folding_evidence_found"] is False
    assert cert["folding_problem_solved"] is False


def test_required_controls_all_pass() -> None:
    gate = _load(GATE, "v36_gate_controls")
    controls = gate.run_v36_controls()
    assert len(controls) == 9
    assert [row["passed"] for row in controls] == [True] * 9
    assert {row["control_id"] for row in controls} == {
        "kcsa_v33_v34_coordinate_csvs_blocked",
        "kcsa_generic_channel_annotation_partial",
        "xcl1_one_state_only_partial",
        "xcl1_mixed_state_pooling_blocked",
        "snca_single_fold_grammar_blocked",
        "snca_without_disorder_evidence_partial",
        "placeholder_citation_blocked",
        "internal_runtime_source_blocked",
        "native_coordinate_metrics_before_selection_blocked",
    }


def test_kcsa_coordinate_csv_control_blocks() -> None:
    gate = _load(GATE, "v36_gate_kcsa_block")
    packages = gate._ready_fixture_packages()
    row = packages["KcsA"]["source_manifest"]["rows"][0]
    row["source_type"] = "RCSB coordinate contacts"
    row["source_url_or_citation"] = "data/external_constraints/KcsA/pore_filter/kcsa_1bl8_pore_filter_external_contacts.csv"
    row["coordinate_derived"] = True
    cert = gate._evaluate_control_packages(packages)
    assert cert["control_status"] == "V36_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"
    assert "KcsA:row_0:coordinate_or_prediction_leakage" in cert["failed_checks"]


def test_kcsa_generic_channel_annotation_is_partial() -> None:
    gate = _load(GATE, "v36_gate_kcsa_partial")
    packages = gate._ready_fixture_packages()
    packages["KcsA"]["source_manifest"]["rows"] = [gate._minimal_row("KcsA", "sequence_or_family_identity")]
    cert = gate._evaluate_control_packages(packages)
    assert cert["control_status"] == "V36_PARTIAL_EVIDENCE_DOSSIERS_CLEAN_ABSTAIN"
    assert "KcsA:missing_bucket:ion_selectivity_context" in cert["failed_checks"]


def test_xcl1_one_state_partial_and_mixed_pooling_blocked() -> None:
    gate = _load(GATE, "v36_gate_xcl1")
    packages = gate._ready_fixture_packages()
    packages["XCL1_lymphotactin"]["source_manifest"]["rows"] = [
        row for row in packages["XCL1_lymphotactin"]["source_manifest"]["rows"]
        if row["evidence_bucket"] != "state_B_function_context"
    ]
    cert = gate._evaluate_control_packages(packages)
    assert cert["control_status"] == "V36_PARTIAL_EVIDENCE_DOSSIERS_CLEAN_ABSTAIN"
    assert "XCL1_lymphotactin:missing_bucket:state_B_function_context" in cert["failed_checks"]

    packages = gate._ready_fixture_packages()
    packages["XCL1_lymphotactin"]["evidence_dossier"]["grammar_rules"]["mixed_state_pooling_allowed"] = True
    cert = gate._evaluate_control_packages(packages)
    assert cert["control_status"] == "V36_BLOCKED_PLACEHOLDER_OR_UNTRUSTED_SOURCES"
    assert "XCL1_lymphotactin:mixed_state_pooling_forbidden" in cert["failed_checks"]


def test_snca_controls_for_idp_grammar() -> None:
    gate = _load(GATE, "v36_gate_snca")
    packages = gate._ready_fixture_packages()
    packages["alpha_synuclein_SNCA"]["evidence_dossier"]["grammar_rules"]["single_native_fold_model_allowed"] = True
    cert = gate._evaluate_control_packages(packages)
    assert cert["control_status"] == "V36_BLOCKED_PLACEHOLDER_OR_UNTRUSTED_SOURCES"
    assert "alpha_synuclein_SNCA:single_native_fold_grammar_forbidden_for_SNCA" in cert["failed_checks"]

    packages = gate._ready_fixture_packages()
    packages["alpha_synuclein_SNCA"]["source_manifest"]["rows"] = [
        row for row in packages["alpha_synuclein_SNCA"]["source_manifest"]["rows"]
        if row["evidence_bucket"] != "intrinsic_disorder_context"
    ]
    cert = gate._evaluate_control_packages(packages)
    assert cert["control_status"] == "V36_PARTIAL_EVIDENCE_DOSSIERS_CLEAN_ABSTAIN"
    assert "alpha_synuclein_SNCA:missing_bucket:intrinsic_disorder_context" in cert["failed_checks"]


def test_placeholder_internal_and_native_metric_blocks() -> None:
    gate = _load(GATE, "v36_gate_blocks")

    packages = gate._ready_fixture_packages()
    packages["KcsA"]["source_manifest"]["rows"][0]["source_name"] = "TODO placeholder"
    packages["KcsA"]["source_manifest"]["rows"][0]["source_url_or_citation"] = "citation needed"
    cert = gate._evaluate_control_packages(packages)
    assert cert["control_status"] == "V36_BLOCKED_PLACEHOLDER_OR_UNTRUSTED_SOURCES"
    assert "KcsA:row_0:placeholder_or_missing_source_name" in cert["failed_checks"]

    packages = gate._ready_fixture_packages()
    packages["KcsA"]["source_manifest"]["rows"][0]["file_path"] = "first_contact_clean_pharmacotopology_layer_run/V35/report.json"
    packages["KcsA"]["source_manifest"]["rows"][0]["internal_runtime_source"] = True
    cert = gate._evaluate_control_packages(packages)
    assert cert["control_status"] == "V36_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"
    assert "KcsA:row_0:internal_runtime_source_supplied" in cert["failed_checks"]

    packages = gate._ready_fixture_packages()
    packages["KcsA"]["source_manifest"]["rows"][0]["native_metrics_used_for_selection"] = True
    cert = gate._evaluate_control_packages(packages)
    assert cert["control_status"] == "V36_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"
    assert "KcsA:row_0:native_metrics_used_for_selection" in cert["failed_checks"]


def test_writer_outputs_artifacts(tmp_path: Path) -> None:
    gate = _load(GATE, "v36_gate_writer")
    cert = gate.build_v36(gate._ready_fixture_packages())
    paths = gate.write_outputs(tmp_path / "out", cert)
    for path in paths.values():
        assert path.exists()
    written = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert written["next_decision"]["claim_allowed"] is False
    assert written["artifacts"]["certificate"] == str(paths["certificate"])
