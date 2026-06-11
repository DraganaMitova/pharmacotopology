from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v28_xcl1_state_condition_evidence_contrast_test_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v28_xcl1_contrast", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _v27() -> dict:
    return {
        "test_status": "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION_PASSED_CLAIM_DISABLED",
        "claim_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "condition_label_context_locked": True,
        "state_A_context_label_locked": True,
        "state_B_context_label_locked": True,
        "state_labels_locked": True,
        "monomer_dimer_context_locked": True,
        "mixed_state_leakage_guard_preserved": True,
        "state_specific_couplings_present": False,
        "forbidden_misclassification_violations": [],
        "condition_label_manifest": {
            "state_A_label": "state_A_chemokine_like_or_monomer_context",
            "state_B_label": "state_B_alternative_beta_sandwich_or_dimer_context",
        },
        "state_specific_coupling_availability": {
            "state_specific_couplings_present": False,
            "external_coupling_or_constraint_files": [],
        },
        "next_decision": {
            "selected_next_panel": "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST",
            "claim_allowed": False,
        },
    }


def _v26() -> dict:
    return {
        "test_status": "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST_PASSED_CLAIM_DISABLED",
        "state_separation_operator_passed": True,
        "claim_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "forbidden_misclassification_violations": [],
    }


def _v20() -> dict:
    return {
        "test_status": "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED",
        "state_A_role_evidence_found": True,
        "state_B_role_evidence_found": True,
        "state_specific_role_evidence_found": True,
        "claim_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "mixed_state_pollution": False,
        "single_fold_claim_made": False,
        "fold_switch_claim_made": False,
        "forbidden_misclassification_violations": [],
        "state_specific_readouts": {
            "state_A": {
                "pdb_id": "2HDM",
                "state_role_bucket": "state_A_chemokine_monomer_support_context",
                "expected_context": "chemokine_like_or_monomer_context",
                "state_role_evidence_found": True,
                "chain_ca_counts": {"A": 1480},
                "model_count": 20,
            },
            "state_B": {
                "pdb_id": "2JP1",
                "state_role_bucket": "state_B_beta_sandwich_dimer_or_alternative_state_support_context",
                "expected_context": "alternative_state_or_dimer_context",
                "state_role_evidence_found": True,
                "chain_ca_counts": {"A": 1200, "B": 1200},
                "model_count": 20,
                "chain_interface_readout": {"interface_pairs_present": ["A-B"]},
            },
        },
    }


def test_v28_passes_state_condition_contrast_without_folding_claim() -> None:
    mod = _load_module()
    cert = mod.build_v28(_v27(), _v26(), _v20())
    assert cert["test_status"] == "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST_PASSED_CLAIM_DISABLED"
    assert cert["state_A_condition_evidence_found"] is True
    assert cert["state_B_condition_evidence_found"] is True
    assert cert["state_condition_contrast_preserved"] is True
    assert cert["positive_pressure_evidence_found"] is True
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False


def test_v28_preserves_mixed_state_and_single_fold_guards() -> None:
    mod = _load_module()
    cert = mod.build_v28(_v27(), _v26(), _v20())
    assert cert["mixed_state_pollution"] is False
    assert cert["mixed_state_contact_pooling_used"] is False
    assert cert["single_fold_forcing"] is False
    assert cert["single_fold_claim_made"] is False
    assert cert["fold_switch_claim_made"] is False
    assert cert["forbidden_misclassification_violations"] == []
    assert cert["state_condition_contrast_readout"]["guard"]["selection_threshold_used"] is False


def test_v28_does_not_require_couplings_for_condition_contrast() -> None:
    mod = _load_module()
    cert = mod.build_v28(_v27(), _v26(), _v20())
    assert cert["state_specific_couplings_present"] is False
    assert cert["external_couplings_required_for_V28"] is False
    assert cert["external_couplings_required_for_future_MD_or_contact_test"] is True
    assert "state_specific_external_couplings_or_constraints_if_available" in cert["missing_evidence"]


def test_v28_blocks_if_state_b_label_missing() -> None:
    mod = _load_module()
    v27 = _v27()
    v27["state_B_context_label_locked"] = False
    cert = mod.build_v28(v27, _v26(), _v20())
    assert cert["test_status"] == "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST_BLOCKED_OR_CLEAN_ABSTAIN_CLAIM_DISABLED"
    assert "V27_condition_state_label_context_locked" in cert["failed_checks"] or "state_B_condition_evidence_found" in cert["failed_checks"]
    assert cert["claim_allowed"] is False


def test_v28_writer_outputs_certificate_and_decision(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_v28(_v27(), _v26(), _v20())
    paths = mod.write_outputs(tmp_path, cert)
    assert paths["certificate"].exists()
    assert paths["state_condition_contrast_readout"].exists()
    assert paths["next_decision"].exists()
    assert paths["report"].exists()
    written = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert written["test_status"] == "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST_PASSED_CLAIM_DISABLED"
    assert written["positive_folding_evidence_found"] is False
