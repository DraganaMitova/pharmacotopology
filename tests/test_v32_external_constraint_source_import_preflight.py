from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v32_external_constraint_source_import_preflight_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v32_import_preflight", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _v31_clean_abstain() -> dict:
    return {
        "preflight_status": "V31_CLEAN_ABSTAIN_NO_REAL_EXTERNAL_CONSTRAINTS_FOR_SELECTED_TARGETS",
        "claim_allowed": False,
        "new_md_executed": False,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "native_metrics_used_for_selection": False,
        "selected_V31_targets": ["XCL1_lymphotactin", "KcsA"],
    }


def test_v32_clean_abstains_when_no_import_manifest_rows() -> None:
    mod = _load_module()
    cert = mod.build_v32(_v31_clean_abstain(), {"manifest_present": False, "rows": []})
    assert cert["preflight_status"] == "V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED"
    assert cert["selected_V33_target"] is None
    assert cert["selected_V33_panel"] is None
    assert cert["claim_allowed"] is False
    assert cert["new_MD_allowed"] is False
    assert cert["new_MD_recommended"] is False
    assert cert["provenance_clean"] is True
    assert all(row["target_ready_for_V33"] is False for row in cert["target_rows"])


def test_v32_blocks_internal_runtime_report_as_external_source() -> None:
    mod = _load_module()
    manifest = {
        "manifest_present": True,
        "rows": [
            {
                "target": "KcsA",
                "file_path": "first_contact_clean_pharmacotopology_layer_run/V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION/v24_kcsa_external_coupling_availability.json",
                "evidence_type": "external_coupling",
                "state_or_context": "pore_filter_interface",
                "source_name": "internal_report",
                "source_url_or_citation": "internal runtime output",
            }
        ],
    }
    cert = mod.build_v32(_v31_clean_abstain(), manifest)
    assert cert["preflight_status"] == "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_BLOCKED"
    assert cert["provenance_clean"] is False
    assert "internal_generated_artifact_supplied_as_external_source" in cert["preflight_failed_checks"]
    assert cert["selected_V33_target"] is None
    assert cert["claim_allowed"] is False


def test_v32_imported_kcsa_constraints_can_select_kcsa_for_v33(tmp_path: Path) -> None:
    mod = _load_module()
    mod.REPO_ROOT = tmp_path
    f1 = tmp_path / "data" / "external_constraints" / "KcsA" / "pore_filter" / "kcsa_test_pore.csv"
    f2 = tmp_path / "data" / "external_constraints" / "KcsA" / "assembly_interface" / "kcsa_test_interface.csv"
    f1.parent.mkdir(parents=True, exist_ok=True)
    f2.parent.mkdir(parents=True, exist_ok=True)
    f1.write_text("i,j,score\n1,2,0.9\n", encoding="utf-8")
    f2.write_text("chain_a,chain_b,score\nA,B,0.8\n", encoding="utf-8")
    manifest = {
        "manifest_present": True,
        "rows": [
            {
                "target": "KcsA",
                "file_path": "data/external_constraints/KcsA/pore_filter/kcsa_test_pore.csv",
                "evidence_type": "pore_filter_coupling",
                "state_or_context": "TVGYG_selectivity_filter_pore_coupling",
                "source_name": "unit_test_external_source",
                "source_url_or_citation": "unit-test-citation",
                "source_date_or_version": "test",
            },
            {
                "target": "KcsA",
                "file_path": "data/external_constraints/KcsA/assembly_interface/kcsa_test_interface.csv",
                "evidence_type": "assembly_interface_constraint",
                "state_or_context": "tetramer_chain_interface_context",
                "source_name": "unit_test_external_source",
                "source_url_or_citation": "unit-test-citation",
                "source_date_or_version": "test",
            },
        ],
    }
    cert = mod.build_v32(_v31_clean_abstain(), manifest)
    assert cert["preflight_status"] == "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED"
    assert cert["selected_V33_target"] == "KcsA"
    krow = next(row for row in cert["target_rows"] if row["target"] == "KcsA")
    assert krow["target_ready_for_V33"] is True
    assert len(krow["valid_real_external_constraint_rows"]) == 2
    assert cert["new_MD_allowed"] is False
    assert cert["claim_allowed"] is False


def test_v32_xcl1_requires_both_state_buckets(tmp_path: Path) -> None:
    mod = _load_module()
    mod.REPO_ROOT = tmp_path
    f1 = tmp_path / "data" / "external_constraints" / "XCL1_lymphotactin" / "state_A" / "xcl1_test_a.csv"
    f2 = tmp_path / "data" / "external_constraints" / "XCL1_lymphotactin" / "state_B" / "xcl1_test_b.csv"
    f1.parent.mkdir(parents=True, exist_ok=True)
    f2.parent.mkdir(parents=True, exist_ok=True)
    f1.write_text("i,j,score\n1,2,0.9\n", encoding="utf-8")
    f2.write_text("i,j,score\n3,4,0.8\n", encoding="utf-8")

    one_state = {
        "manifest_present": True,
        "rows": [
            {
                "target": "XCL1_lymphotactin",
                "file_path": "data/external_constraints/XCL1_lymphotactin/state_A/xcl1_test_a.csv",
                "evidence_type": "state_specific_constraint",
                "state_or_context": "state_A_chemokine_like",
                "source_name": "unit_test_external_source",
                "source_url_or_citation": "unit-test-citation",
            }
        ],
    }
    cert_one = mod.build_v32(_v31_clean_abstain(), one_state)
    xrow_one = next(row for row in cert_one["target_rows"] if row["target"] == "XCL1_lymphotactin")
    assert xrow_one["target_ready_for_V33"] is False
    assert "state_B_real_external_constraint" in xrow_one["missing_for_V33"]

    two_state = {
        "manifest_present": True,
        "rows": [
            one_state["rows"][0],
            {
                "target": "XCL1_lymphotactin",
                "file_path": "data/external_constraints/XCL1_lymphotactin/state_B/xcl1_test_b.csv",
                "evidence_type": "state_specific_constraint",
                "state_or_context": "state_B_dimer_beta_sandwich_like",
                "source_name": "unit_test_external_source",
                "source_url_or_citation": "unit-test-citation",
            },
        ],
    }
    cert_two = mod.build_v32(_v31_clean_abstain(), two_state)
    assert cert_two["preflight_status"] == "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED"
    assert cert_two["selected_V33_target"] == "XCL1_lymphotactin"


def test_v32_writes_template_and_outputs(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_v32(_v31_clean_abstain(), {"manifest_present": False, "rows": []})
    manifest = tmp_path / "data" / "external_constraints" / "v32_external_constraint_source_import_manifest.json"
    paths = mod.write_outputs(tmp_path / "out", cert, manifest)
    for key in ["certificate", "import_template", "classified_import_rows", "target_rows", "target_rows_csv", "import_requirements", "decision", "report"]:
        assert paths[key].exists(), key
    template = json.loads(paths["import_template"].read_text())
    assert template["claim_allowed"] is False
    assert len(template["rows"]) >= 4
    written = json.loads(paths["certificate"].read_text())
    assert written["artifacts"]["decision"] == str(paths["decision"])
