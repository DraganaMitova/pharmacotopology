from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v20_xcl1_state_specific_evidence_readout_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v20_xcl1", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _atom(serial: int, chain: str, res: int, resname: str, x: float, y: float, z: float) -> str:
    return (
        f"ATOM  {serial:5d}  CA  {resname:>3s} {chain}{res:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
    )


def _state_a_pdb(tmp_path: Path) -> Path:
    pdb = tmp_path / "2HDM.pdb"
    lines = [
        "HEADER    CHEMOKINE-LIKE XCL1 STATE A\n",
        "TITLE     XCL1 LYMPHOTACTIN CHEMOKINE MONOMER CONTEXT\n",
        "MODEL        1\n",
    ]
    serial = 1
    for i in range(1, 9):
        lines.append(_atom(serial, "A", i, "ALA", float(i), 0.0, 0.0)); serial += 1
    lines.append("ENDMDL\n")
    pdb.write_text("".join(lines), encoding="utf-8")
    return pdb


def _state_b_pdb(tmp_path: Path) -> Path:
    pdb = tmp_path / "2JP1.pdb"
    lines = [
        "HEADER    ALTERNATIVE XCL1 STATE B DIMER\n",
        "TITLE     XCL1 LYMPHOTACTIN BETA SANDWICH DIMER CONTEXT\n",
        "MODEL        1\n",
    ]
    serial = 1
    for i in range(1, 9):
        lines.append(_atom(serial, "A", i, "ALA", float(i), 0.0, 0.0)); serial += 1
    for i in range(1, 9):
        lines.append(_atom(serial, "B", i, "ALA", float(i), 4.0, 0.0)); serial += 1
    lines.append("ENDMDL\n")
    pdb.write_text("".join(lines), encoding="utf-8")
    return pdb


def _preflight(a: Path, b: Path) -> dict:
    return {
        "data_preflight_status": "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT",
        "target_results": [
            {
                "target_id": "XCL1_lymphotactin",
                "target_material_status": "ready_for_zero_md_role_readout",
                "expected_role_class": "metamorphic_switch_object",
                "material_results": [
                    {
                        "material_id": "XCL1_state_A_chemokine_like_structure",
                        "material_status": "present_and_provenanced",
                        "path": str(a),
                    },
                    {
                        "material_id": "XCL1_state_B_alternative_conformation_structure",
                        "material_status": "present_and_provenanced",
                        "path": str(b),
                    },
                ],
            }
        ],
    }


def _zero_md(pass_role: bool = True) -> dict:
    return {
        "role_transfer_status": "V16_ZERO_MD_ROLE_TRANSFER_READOUT_COMPLETED_CLAIM_DISABLED",
        "role_classification_passed_targets": ["XCL1_lymphotactin"] if pass_role else [],
        "pressure_role_transfer_passed_targets": ["XCL1_lymphotactin"] if pass_role else [],
        "target_results": [
            {
                "target_id": "XCL1_lymphotactin",
                "expected_role_class": "metamorphic_switch_object",
                "role_classification_passed": pass_role,
                "positive_role_context_found": pass_role,
                "forbidden_misclassification_violations": [],
            }
        ],
    }


def _gap_lock() -> dict:
    return {
        "gap_lock_status": "V16_PRESSURE_EVIDENCE_GAP_LOCKED",
        "target_results": [
            {
                "target_id": "XCL1_lymphotactin",
                "missing_for_evidence_test": ["state_specific_external_couplings_or_constraints_if_available"],
            }
        ],
    }


def test_v20_xcl1_passes_state_specific_readout_without_mixed_core(tmp_path: Path) -> None:
    mod = _load_module()
    a = _state_a_pdb(tmp_path)
    b = _state_b_pdb(tmp_path)
    cert = mod.build_evidence(_preflight(a, b), _zero_md(), _gap_lock(), a, b)
    assert cert["test_status"] == "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED"
    assert cert["positive_pressure_evidence_found"] is True
    assert cert["state_specific_role_evidence_found"] is True
    assert cert["state_A_role_evidence_found"] is True
    assert cert["state_B_role_evidence_found"] is True
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["mixed_state_pollution"] is False
    assert cert["mixed_state_contact_pooling_used"] is False
    assert cert["single_fold_claim_made"] is False
    assert cert["fold_switch_claim_made"] is False
    assert cert["forbidden_misclassification_violations"] == []
    assert "state_A_chemokine_monomer_support" in cert["role_buckets_assigned"]
    assert "state_B_beta_sandwich_dimer_support" in cert["role_buckets_assigned"]
    assert cert["state_separation_guard"]["selection_threshold_used"] is False


def test_v20_xcl1_clean_abstains_if_state_b_lacks_dimer_context(tmp_path: Path) -> None:
    mod = _load_module()
    a = _state_a_pdb(tmp_path)
    b = tmp_path / "2JP1_single_chain.pdb"
    lines = ["HEADER    ALTERNATIVE XCL1 STATE B SINGLE CHAIN\n"]
    serial = 1
    for i in range(1, 9):
        lines.append(_atom(serial, "A", i, "ALA", float(i), 0.0, 0.0)); serial += 1
    b.write_text("".join(lines), encoding="utf-8")
    cert = mod.build_evidence(_preflight(a, b), _zero_md(), _gap_lock(), a, b)
    assert cert["test_status"] == "V20_XCL1_CLEAN_ABSTAIN_STATE_SPECIFIC_EVIDENCE_INSUFFICIENT_CLAIM_DISABLED"
    assert cert["positive_pressure_evidence_found"] is False
    assert cert["clean_abstain_valid"] is True
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False
    assert "monomer_dimer_context" in cert["missing_evidence"]


def test_v20_xcl1_forbidden_zero_md_violation_blocks_evidence(tmp_path: Path) -> None:
    mod = _load_module()
    a = _state_a_pdb(tmp_path)
    b = _state_b_pdb(tmp_path)
    zero = _zero_md()
    zero["target_results"][0]["forbidden_misclassification_violations"] = ["mixing_two_states_into_one_false_core"]
    cert = mod.build_evidence(_preflight(a, b), zero, _gap_lock(), a, b)
    assert cert["test_status"] == "V20_XCL1_CLEAN_ABSTAIN_STATE_SPECIFIC_EVIDENCE_INSUFFICIENT_CLAIM_DISABLED"
    assert "mixing_two_states_into_one_false_core" in cert["forbidden_misclassification_violations"]
    assert cert["positive_pressure_evidence_found"] is False
    assert cert["positive_folding_evidence_found"] is False


def test_v20_xcl1_writer_outputs_certificate_and_state_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    a = _state_a_pdb(tmp_path)
    b = _state_b_pdb(tmp_path)
    cert = mod.build_evidence(_preflight(a, b), _zero_md(), _gap_lock(), a, b)
    out = tmp_path / "out"
    mod.write_outputs(out, cert)
    for name in [
        "v20_xcl1_state_specific_evidence_readout_certificate.json",
        "v20_xcl1_state_A_readout.json",
        "v20_xcl1_state_B_readout.json",
        "v20_xcl1_state_separation_guard.json",
        "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_REPORT.md",
    ]:
        assert (out / name).exists()
    written = json.loads((out / "v20_xcl1_state_specific_evidence_readout_certificate.json").read_text(encoding="utf-8"))
    assert written["claim_allowed"] is False
    assert written["positive_folding_evidence_found"] is False
