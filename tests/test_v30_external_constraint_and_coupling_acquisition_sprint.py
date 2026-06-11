from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v30_external_constraint_and_coupling_acquisition_sprint_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v30_acquisition", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _row(target: str, missing: list[str], ops: list[str] | None = None) -> dict:
    return {
        "target": target,
        "pressure_type": f"{target}_pressure",
        "role_detected": True,
        "selected_signal": f"{target}_signal",
        "missing_next_evidence": missing,
        "mechanism_operator_used": ops or ["role_detect", "evidence_assign", "pollution_guard"],
        "claim_allowed": False,
    }


def _v29() -> dict:
    targets = [
        _row("XCL1_lymphotactin", ["state_specific_external_couplings_or_constraints_if_available", "condition_labels_if_available"], ["role_detect", "context_detect", "evidence_assign", "state_separation", "pollution_guard"]),
        _row("KcsA", ["external_couplings_if_available", "pore_filter_coupling_support"], ["role_detect", "context_detect", "evidence_assign", "topology_context", "interface_context", "pollution_guard"]),
        _row("p53_TAD_MDM2", ["external_couplings_if_available"], ["role_detect", "context_detect", "evidence_assign", "interface_context", "clean_abstain", "pollution_guard"]),
        _row("4AKE", ["interdomain_hinge_or_closure_evidence"]),
        _row("1UBQ", ["DCA_background_enrichment"]),
        _row("1CLL", ["interdomain_hinge_evidence_if_claim_needed"]),
    ]
    md_rows = [
        {"target": r["target"], "md_blockers": r["missing_next_evidence"], "new_MD_allowed": False, "new_MD_recommended": False}
        for r in targets
    ]
    return {
        "summary_status": "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_LOCKED",
        "md_readiness_decision_status": "V29_NO_NEW_MD_READY_EXTERNAL_CONSTRAINTS_REQUIRED",
        "claim_allowed": False,
        "new_md_executed": False,
        "positive_pressure_evidence_targets": ["p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"],
        "positive_folding_evidence_targets": [],
        "mechanism_operator_panel_summary": {"targets": targets},
        "md_readiness_rows": md_rows,
    }


def test_v30_locks_acquisition_without_md_or_folding_claim(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_v30(_v29(), tmp_path)
    assert cert["sprint_status"] == "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT_LOCKED"
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["new_MD_allowed"] is False
    assert cert["new_MD_recommended"] is False
    assert cert["positive_folding_evidence_targets"] == []
    assert cert["md_ready_targets"] == []
    assert cert["selected_next_panel"] == "V31_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_AND_PREFLIGHT_SPRINT"
    assert cert["selected_V31_targets"] == ["XCL1_lymphotactin", "KcsA"]


def test_v30_scans_target_specific_constraint_files(tmp_path: Path) -> None:
    mod = _load_module()
    d = tmp_path / "data" / "external_constraints" / "XCL1_lymphotactin"
    d.mkdir(parents=True)
    f = d / "xcl1_state_specific_constraints.csv"
    f.write_text("target,XCL1,state_specific_constraints\n", encoding="utf-8")
    cert = mod.build_v30(_v29(), tmp_path)
    xcl1 = next(row for row in cert["target_rows"] if row["target"] == "XCL1_lymphotactin")
    assert xcl1["local_constraints_or_couplings_present"] is True
    assert any("xcl1_state_specific_constraints.csv" in path for path in xcl1["candidate_coupling_or_constraint_files"])
    assert "XCL1_lymphotactin" in cert["constraint_or_coupling_ready_targets"]


def test_v30_writes_all_required_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_v30(_v29(), tmp_path)
    paths = mod.write_outputs(tmp_path / "out", cert)
    for key in ["certificate", "manifest", "preflight", "coupling_scan", "decision", "rows", "rows_csv", "report"]:
        assert paths[key].exists(), key
    written = json.loads(paths["certificate"].read_text())
    assert written["artifacts"]["decision"] == str(paths["decision"])
    assert json.loads(paths["decision"].read_text())["new_MD_allowed"] is False
