from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v27_xcl1_condition_and_coupling_evidence_acquisition_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v27_xcl1_condition_coupling", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _v26() -> dict:
    return {
        "test_status": "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST_PASSED_CLAIM_DISABLED",
        "state_separation_operator_passed": True,
        "claim_allowed": False,
        "new_md_executed": False,
        "operator_readout": {
            "state_A_detected": True,
            "state_B_detected": True,
            "state_specific_role_evidence_found": True,
            "state_separation_guard_passed": True,
            "mixed_state_pollution": False,
            "mixed_state_contact_pooling_used": False,
            "single_fold_forcing": False,
            "fold_switch_claim_made": False,
            "available_evidence": ["monomer_dimer_context", "leakage_guard_preventing_mixed_state_fake_core"],
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
        "monomer_dimer_context_present": True,
        "mixed_state_pollution": False,
        "single_fold_claim_made": False,
        "fold_switch_claim_made": False,
        "available_evidence": ["state_A_and_state_B_context_structures", "monomer_dimer_context"],
        "missing_evidence": ["state_specific_external_couplings_or_constraints_if_available"],
    }


def _preflight() -> dict:
    return {
        "data_preflight_status": "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT",
        "target_results": [
            {
                "target_id": "XCL1_lymphotactin",
                "target_material_status": "ready_for_zero_md_role_readout",
            }
        ],
    }


def test_v27_passes_condition_context_without_synthesizing_couplings() -> None:
    mod = _load_module()
    cert = mod.build_v27(_v26(), _v20(), _preflight())
    assert cert["test_status"] == "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION_PASSED_CLAIM_DISABLED"
    assert cert["condition_label_context_locked"] is True
    assert cert["state_A_context_label_locked"] is True
    assert cert["state_B_context_label_locked"] is True
    assert cert["mixed_state_leakage_guard_preserved"] is True
    assert cert["positive_pressure_evidence_found"] is True
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["external_couplings_required_for_V27"] is False
    assert cert["external_couplings_required_for_future_MD_or_contact_test"] is True
    assert cert["next_decision"]["selected_next_panel"] == "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST"
    assert cert["next_decision"]["new_MD_allowed"] is False


def test_v27_blocks_if_v26_state_separation_not_clean() -> None:
    mod = _load_module()
    v26 = _v26()
    v26["operator_readout"]["mixed_state_pollution"] = True
    cert = mod.build_v27(v26, _v20(), _preflight())
    assert cert["test_status"] == "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION_CLEAN_ABSTAIN_OR_BLOCKED_CLAIM_DISABLED"
    assert cert["condition_label_context_locked"] is False
    assert "mixed_state_leakage_guard_preserved" in cert["failed_checks"]
    assert cert["claim_allowed"] is False


def test_v27_writes_required_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_v27(_v26(), _v20(), _preflight())
    paths = mod.write_outputs(tmp_path, cert)
    expected = {
        "certificate",
        "condition_label_manifest",
        "condition_preflight",
        "state_specific_coupling_availability",
        "next_decision",
        "report",
    }
    assert expected.issubset(paths)
    for p in paths.values():
        assert p.exists()
    written = json.loads(paths["certificate"].read_text())
    assert written["artifacts"]["condition_label_manifest"] == str(paths["condition_label_manifest"])
    manifest = json.loads(paths["condition_label_manifest"].read_text())
    assert manifest["state_labels_locked"] is True
    coupling = json.loads(paths["state_specific_coupling_availability"].read_text())
    assert coupling["external_couplings_required_for_V27"] is False
