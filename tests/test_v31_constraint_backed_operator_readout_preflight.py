from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v31_constraint_backed_operator_readout_preflight_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v31_preflight", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _target_row(target: str, candidates: list[str], annotations: list[str] | None = None) -> dict:
    return {
        "target": target,
        "pressure_type": f"{target}_pressure",
        "selected_signal": f"{target}_signal",
        "acquisition_status": "local_constraint_or_coupling_files_present_ready_for_next_readout_preflight",
        "candidate_coupling_or_constraint_files": candidates,
        "candidate_annotation_files": annotations or [],
        "raw_local_external_annotation_files": [],
        "ready_for_MD": False,
        "new_MD_allowed": False,
        "claim_allowed": False,
    }


def _v30_no_real_selected_scope() -> dict:
    return {
        "sprint_status": "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT_LOCKED",
        "claim_allowed": False,
        "new_md_executed": False,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "native_metrics_used_for_selection": False,
        "selected_V31_targets": ["XCL1_lymphotactin", "KcsA"],
        "target_rows": [
            _target_row(
                "XCL1_lymphotactin",
                ["first_contact_clean_pharmacotopology_layer_run/V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION/v27_xcl1_state_specific_coupling_availability.json"],
                ["first_contact_clean_pharmacotopology_layer_run/V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT/v20_xcl1_state_A_readout.json"],
            ),
            _target_row(
                "KcsA",
                ["first_contact_clean_pharmacotopology_layer_run/V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION/v24_kcsa_external_coupling_availability.json"],
                ["first_contact_clean_pharmacotopology_layer_run/V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT/v19_kcsa_pore_filter_readout.json"],
            ),
            _target_row("4AKE", ["data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json"]),
        ],
    }


def test_v31_clean_abstains_when_selected_targets_only_have_internal_reports() -> None:
    mod = _load_module()
    cert = mod.build_v31(_v30_no_real_selected_scope())
    assert cert["preflight_status"] == "V31_CLEAN_ABSTAIN_NO_REAL_EXTERNAL_CONSTRAINTS_FOR_SELECTED_TARGETS"
    assert cert["selected_V32_target"] is None
    assert cert["selected_V32_panel"] is None
    assert cert["claim_allowed"] is False
    assert cert["new_MD_allowed"] is False
    assert cert["new_MD_recommended"] is False
    assert cert["provenance_clean"] is True
    xcl1 = next(row for row in cert["target_rows"] if row["target"] == "XCL1_lymphotactin")
    assert xcl1["constraint_backed_operator_readout_allowed"] is False
    assert xcl1["generated_internal_report_files"]
    assert all(item["allowed_use"] == "allowed_for_audit_only" for item in xcl1["classified_candidate_files"])


def test_v31_passes_and_auto_selects_scope_target_with_real_external_constraint() -> None:
    mod = _load_module()
    v30 = _v30_no_real_selected_scope()
    k = next(row for row in v30["target_rows"] if row["target"] == "KcsA")
    k["candidate_coupling_or_constraint_files"] = ["data/external_constraints/KcsA/kcsa_external_couplings.csv"]
    cert = mod.build_v31(v30)
    assert cert["preflight_status"] == "V31_CONSTRAINT_BACKED_OPERATOR_READOUT_PREFLIGHT_PASSED_CLAIM_DISABLED"
    assert cert["selected_V32_target"] == "KcsA"
    assert cert["selected_V32_panel"] == "V32_FIRST_CONSTRAINT_BACKED_OPERATOR_READOUT"
    assert cert["real_external_constraint_targets_in_selected_scope"] == ["KcsA"]
    krow = next(row for row in cert["target_rows"] if row["target"] == "KcsA")
    assert krow["real_external_constraint_or_coupling_files"] == ["data/external_constraints/KcsA/kcsa_external_couplings.csv"]
    assert krow["allowed_use"] == "allowed_for_constraint_backed_operator_readout"


def test_v31_annotation_only_is_role_context_not_constraint_claim() -> None:
    mod = _load_module()
    c = mod.classify_candidate_file("data/external_annotations/XCL1/xcl1_state_annotation.json")
    assert c["evidence_source_class"] == "annotation_only_external_context"
    assert c["allowed_use"] == "allowed_for_role_context_only"
    assert c["counts_as_real_external_constraint"] is False


def test_v31_alignment_source_is_derivation_preflight_only() -> None:
    mod = _load_module()
    c = mod.classify_candidate_file("external_msa/4ake_pfam00406/PF00406_full.sto")
    assert c["evidence_source_class"] == "real_external_alignment_source"
    assert c["allowed_use"] == "allowed_for_constraint_derivation_preflight_only"
    assert c["counts_as_real_external_constraint"] is False


def test_v31_writes_all_required_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_v31(_v30_no_real_selected_scope())
    paths = mod.write_outputs(tmp_path / "out", cert)
    for key in ["certificate", "classification", "allowed_use_matrix", "decision", "rows", "rows_csv", "report"]:
        assert paths[key].exists(), key
    written = json.loads(paths["certificate"].read_text())
    assert written["artifacts"]["decision"] == str(paths["decision"])
    assert json.loads(paths["decision"].read_text())["claim_allowed"] is False
