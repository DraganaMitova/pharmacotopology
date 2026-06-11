from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v29_mechanism_operator_panel_summary_and_md_readiness_decision_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v29_md_readiness", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _operator_row(target: str, ops: list[str], missing: list[str] | None = None) -> dict:
    return {
        "target": target,
        "pressure_type": f"{target}_pressure",
        "role_detected": True,
        "evidence_type": f"{target}_evidence",
        "selected_signal": f"{target}_signal",
        "abstain_reason_if_any": None,
        "missing_next_evidence": missing or [],
        "mechanism_operator_used": ops,
        "claim_allowed": False,
    }


def _v25() -> dict:
    rows = [
        _operator_row("4AKE", ["role_detect", "context_detect", "evidence_assign", "reachability_check", "replica_support", "chemical_confidence", "topology_context", "pollution_guard"], ["interdomain_hinge_or_closure_evidence"]),
        _operator_row("1UBQ", ["role_detect", "evidence_assign", "replica_support", "chemical_confidence", "pollution_guard"], ["DCA_background_enrichment"]),
        _operator_row("1CLL", ["role_detect", "context_detect", "evidence_assign", "replica_support", "chemical_confidence", "topology_context", "pollution_guard"], ["interdomain_hinge_evidence_if_claim_needed"]),
        _operator_row("p53_TAD_MDM2", ["role_detect", "context_detect", "evidence_assign", "interface_context", "clean_abstain", "pollution_guard"], ["external_couplings_if_available"]),
        _operator_row("KcsA", ["role_detect", "context_detect", "evidence_assign", "topology_context", "interface_context", "clean_abstain", "pollution_guard"], ["external_couplings_if_available", "pore_filter_coupling_support"]),
        _operator_row("XCL1_lymphotactin", ["role_detect", "context_detect", "evidence_assign", "state_separation", "interface_context", "clean_abstain", "pollution_guard"], ["state_specific_external_couplings_or_constraints_if_available"]),
    ]
    return {
        "sprint_status": "V25_FAST_MECHANISM_EVIDENCE_SPRINT_LOCKED",
        "claim_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_targets": [],
        "positive_pressure_evidence_targets": ["p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"],
        "unified_mechanism_operator_table": {
            "targets": rows,
            "operator_counts": {},
            "positive_folding_evidence_targets": [],
            "claim_allowed": False,
        },
    }


def _v28() -> dict:
    return {
        "test_status": "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST_PASSED_CLAIM_DISABLED",
        "claim_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "state_condition_contrast_preserved": True,
        "mixed_state_pollution": False,
        "single_fold_claim_made": False,
        "fold_switch_claim_made": False,
    }


def test_v29_locks_operator_panel_without_md_or_folding_claim() -> None:
    mod = _load_module()
    cert = mod.build_v29(_v25(), _v28())
    assert cert["summary_status"] == "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_LOCKED"
    assert cert["operator_panel_target_count"] == 6
    assert cert["universal_core_operator_present"] is True
    assert cert["xcl1_state_condition_operator_preserved"] is True
    assert cert["positive_folding_evidence_targets"] == []
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["new_MD_allowed"] is False
    assert cert["new_MD_recommended"] is False


def test_v29_blocks_md_until_external_constraints_are_locked() -> None:
    mod = _load_module()
    cert = mod.build_v29(_v25(), _v28())
    dec = cert["md_readiness_decision"]
    assert dec["decision_status"] == "V29_NO_NEW_MD_READY_EXTERNAL_CONSTRAINTS_REQUIRED"
    assert dec["selected_next_panel"] == "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT"
    assert dec["new_MD_allowed"] is False
    assert dec["new_MD_recommended"] is False
    required = set(dec["required_before_any_MD"])
    assert "state_specific_external_couplings_or_constraints_if_available" in required
    assert "external_couplings_if_available" in required


def test_v29_writes_required_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_v29(_v25(), _v28())
    paths = mod.write_outputs(tmp_path, cert)
    assert paths["certificate"].exists()
    assert paths["operator_panel"].exists()
    assert paths["md_decision"].exists()
    assert paths["readiness_rows"].exists()
    assert paths["readiness_csv"].exists()
    written = json.loads(paths["certificate"].read_text())
    assert written["artifacts"]["md_decision"] == str(paths["md_decision"])
    assert json.loads(paths["md_decision"].read_text())["new_MD_allowed"] is False
