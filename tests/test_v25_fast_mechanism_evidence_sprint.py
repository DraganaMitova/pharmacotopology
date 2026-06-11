from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v25_fast_mechanism_evidence_sprint_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v25_fast_mechanism", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _v23() -> dict:
    return {
        "test_status": "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST_PASSED_CLAIM_DISABLED",
        "positive_pressure_evidence_found": True,
        "positive_external_annotation_evidence_found": True,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "isolated_TAD_clean_abstain_preserved": True,
        "partner_bound_interface_evidence_preserved": True,
        "missing_evidence": ["external_couplings_if_available"],
    }


def _v24(couplings: bool = False) -> dict:
    return {
        "test_status": "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION_PASSED_CLAIM_DISABLED",
        "positive_pressure_evidence_found": True,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "external_couplings_present": couplings,
        "external_coupling_or_msa_files": ["external_msa/kcsa/example.sto"] if couplings else [],
        "pore_filter_annotation_present": True,
        "transmembrane_role_present": True,
        "membrane_environment_context_present": True,
        "oligomer_or_chain_interface_context_present": True,
        "biological_assembly_context_locked": True,
        "soluble_core_misclassification_avoided": True,
        "missing_evidence": ["external_couplings_if_available"] if not couplings else [],
    }


def _v20() -> dict:
    return {
        "test_status": "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED",
        "positive_pressure_evidence_found": True,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "state_A_role_evidence_found": True,
        "state_B_role_evidence_found": True,
        "state_specific_role_evidence_found": True,
        "mixed_state_pollution": False,
        "single_fold_claim_made": False,
        "mixed_state_contact_pooling_used": False,
        "available_evidence": ["state_A_and_state_B_labels"],
        "missing_evidence": ["state_specific_external_couplings_or_constraints_if_available"],
        "forbidden_misclassification_violations": [],
    }


def _v15() -> dict:
    return {
        "protein_rows": [
            {"target_id": "4AKE", "artifact_status": "present_machine_readable_4ake_role_artifact", "selected_pairs": ["124-135"]},
            {"target_id": "1UBQ", "artifact_status": "present", "selected_pairs": ["23-48"]},
            {"target_id": "1CLL", "artifact_status": "present", "selected_pairs": ["97-133"]},
        ]
    }


def test_v25_builds_fast_sprint_and_selects_xcl1_when_kcsa_couplings_missing() -> None:
    mod = _load_module()
    cert = mod.build_v25(_v15(), _v23(), _v24(couplings=False), _v20())
    assert cert["sprint_status"] == "V25_FAST_MECHANISM_EVIDENCE_SPRINT_LOCKED"
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["positive_folding_evidence_targets"] == []
    assert set(cert["positive_pressure_evidence_targets"]) == {"p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}
    assert cert["kcsa_coupling_interface_readout"]["kcsa_coupling_available"] is False
    assert cert["kcsa_coupling_interface_readout"]["interface_evidence_present"] is True
    assert cert["xcl1_state_specific_readout"]["mixed_state_pollution"] is False
    decision = cert["next_mechanism_test_decision"]
    assert decision["selected_V26_target"] == "XCL1_lymphotactin"
    assert decision["new_MD_allowed"] is False


def test_v25_prefers_kcsa_if_coupling_interface_evidence_is_available() -> None:
    mod = _load_module()
    cert = mod.build_v25(_v15(), _v23(), _v24(couplings=True), _v20())
    assert cert["kcsa_coupling_interface_readout"]["kcsa_coupling_available"] is True
    assert cert["kcsa_coupling_interface_readout"]["pore_filter_coupling_support"] is True
    assert cert["next_mechanism_test_decision"]["selected_V26_target"] == "KcsA"
    assert cert["next_mechanism_test_decision"]["new_MD_allowed"] is False


def test_v25_writes_required_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_v25(_v15(), _v23(), _v24(couplings=False), _v20())
    paths = mod.write_outputs(tmp_path, cert)
    assert paths["certificate"].exists()
    assert paths["kcsa"].exists()
    assert paths["xcl1"].exists()
    assert paths["operator_table"].exists()
    assert paths["decision"].exists()
    written = json.loads(paths["certificate"].read_text())
    assert written["artifacts"]["kcsa"] == str(paths["kcsa"])
    table = json.loads(paths["operator_table"].read_text())
    assert {row["target"] for row in table["targets"]} == {"4AKE", "1UBQ", "1CLL", "p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}
