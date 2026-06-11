from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v26_xcl1_state_separation_operator_test_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v26_xcl1_state_separation", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _v25() -> dict:
    return {
        "sprint_status": "V25_FAST_MECHANISM_EVIDENCE_SPRINT_LOCKED",
        "claim_allowed": False,
        "new_md_executed": False,
        "next_mechanism_test_decision": {
            "selected_V26_target": "XCL1_lymphotactin",
            "selected_V26_test": "XCL1_STATE_SEPARATION_OPERATOR_TEST",
        },
        "xcl1_state_specific_readout": {
            "state_A_detected": True,
            "state_B_detected": True,
            "state_specific_role_evidence_found": True,
            "mixed_state_pollution": False,
            "single_fold_forcing": False,
            "mixed_state_contact_pooling_used": False,
            "fold_switch_claim_made": False,
            "available_evidence": ["state_A_and_state_B_labels"],
            "missing_evidence": ["condition_labels_if_available"],
        },
    }


def _v20() -> dict:
    return {
        "test_status": "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED",
        "claim_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "state_A_role_evidence_found": True,
        "state_B_role_evidence_found": True,
        "state_specific_role_evidence_found": True,
        "mixed_state_pollution": False,
        "single_fold_claim_made": False,
        "mixed_state_contact_pooling_used": False,
        "fold_switch_claim_made": False,
        "available_evidence": ["state_specific_structural_or_contact_evidence"],
        "missing_evidence": ["state_specific_external_couplings_or_constraints_if_available"],
    }


def test_v26_passes_clean_xcl1_state_separation_operator() -> None:
    mod = _load_module()
    cert = mod.build_v26(_v25(), _v20())
    assert cert["test_status"] == "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST_PASSED_CLAIM_DISABLED"
    assert cert["state_separation_operator_passed"] is True
    assert cert["positive_pressure_evidence_found"] is True
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    op = cert["operator_readout"]
    assert op["mixed_state_pollution"] is False
    assert op["single_fold_forcing"] is False
    assert op["mixed_state_contact_pooling_used"] is False
    assert op["selection_threshold_used"] is False
    assert cert["next_mechanism_decision"]["selected_next_panel"] == "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION"
    assert cert["next_mechanism_decision"]["new_MD_allowed"] is False


def test_v26_blocks_if_v25_did_not_select_xcl1() -> None:
    mod = _load_module()
    v25 = _v25()
    v25["next_mechanism_test_decision"]["selected_V26_target"] = "KcsA"
    cert = mod.build_v26(v25, _v20())
    assert cert["test_status"] == "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST_CLEAN_ABSTAIN_OR_BLOCKED_CLAIM_DISABLED"
    assert "V25_selected_XCL1_state_separation_operator" in cert["failed_checks"]
    assert cert["claim_allowed"] is False


def test_v26_writes_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_v26(_v25(), _v20())
    paths = mod.write_outputs(tmp_path, cert)
    assert paths["certificate"].exists()
    assert paths["operator_readout"].exists()
    assert paths["decision"].exists()
    assert paths["report"].exists()
    written = json.loads(paths["certificate"].read_text())
    assert written["artifacts"]["operator_readout"] == str(paths["operator_readout"])
    readout = json.loads(paths["operator_readout"].read_text())
    assert readout["mechanism_operator"] == "state_separation"
