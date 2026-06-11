from __future__ import annotations

import importlib.util
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_v17_pressure_evidence_sprint_lock_v0.py"
spec = importlib.util.spec_from_file_location("run_v17_pressure_evidence_sprint_lock_v0", SCRIPT_PATH)
assert spec is not None and spec.loader is not None
v17 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v17)


def _write_mock_p53_pdb(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    atom_id = 1
    for chain, count, y in [("A", 85, 0.0), ("B", 13, 3.0)]:
        for resi in range(1, count + 1):
            x = float(resi % 10)
            lines.append(
                f"ATOM  {atom_id:5d}  CA  ALA {chain:1s}{resi:4d}    "
                f"{x:8.3f}{y:8.3f}{0.0:8.3f}  1.00 20.00           C\n"
            )
            atom_id += 1
    path.write_text("".join(lines) + "END\n", encoding="utf-8")


def _certs(tmp_path: Path) -> tuple[dict, dict, dict]:
    p53_pdb = tmp_path / "1YCR.pdb"
    _write_mock_p53_pdb(p53_pdb)
    preflight = {
        "kind": "V16_TRANSFER_DATA_PREFLIGHT_v0",
        "data_preflight_status": "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT",
        "claim_allowed": False,
        "new_md_executed": False,
        "target_results": [
            {
                "target_id": "p53_TAD_MDM2",
                "pressure_class": "disorder_partner_induced_binding",
                "target_material_status": "ready_for_zero_md_role_readout",
                "material_results": [
                    {
                        "material_id": "p53_TAD_MDM2_complex_structure",
                        "material_status": "present_and_provenanced",
                        "path": str(p53_pdb),
                        "pdb_summary": {"chain_ca_counts": {"A": 85, "B": 13}},
                    }
                ],
            },
            {
                "target_id": "KcsA",
                "pressure_class": "membrane_pore_oligomer_environment",
                "target_material_status": "ready_for_zero_md_role_readout",
                "material_results": [
                    {
                        "material_id": "KcsA_Fab_complex_structure",
                        "material_status": "present_and_provenanced",
                        "path": str(tmp_path / "1K4C.pdb"),
                        "pdb_summary": {"chain_ca_counts": {"A": 219, "B": 212, "C": 103}},
                    }
                ],
            },
            {
                "target_id": "XCL1_lymphotactin",
                "pressure_class": "metamorphic_fold_switching",
                "target_material_status": "ready_for_zero_md_role_readout",
                "material_results": [
                    {
                        "material_id": "XCL1_state_A_chemokine_like_structure",
                        "material_status": "present_and_provenanced",
                        "path": str(tmp_path / "2HDM.pdb"),
                        "pdb_summary": {"chain_ca_counts": {"A": 1480}},
                    },
                    {
                        "material_id": "XCL1_state_B_alternative_conformation_structure",
                        "material_status": "present_and_provenanced",
                        "path": str(tmp_path / "2JP1.pdb"),
                        "pdb_summary": {"chain_ca_counts": {"A": 1200, "B": 1200}},
                    },
                ],
            },
        ],
    }
    zero = {
        "kind": "V16_ZERO_MD_ROLE_TRANSFER_READOUT_v0",
        "role_transfer_status": "V16_ZERO_MD_ROLE_TRANSFER_READOUT_COMPLETED_CLAIM_DISABLED",
        "role_classification_passed_targets": ["p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"],
        "pressure_role_transfer_passed_targets": ["p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"],
        "positive_folding_evidence_targets": [],
        "claim_allowed": False,
        "new_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "target_specific_threshold_tuning_allowed": False,
        "target_rows": [
            {"target_id": "p53_TAD_MDM2", "role_classification_passed": True, "positive_role_context_found": True, "forbidden_misclassification_violations": []},
            {"target_id": "KcsA", "role_classification_passed": True, "positive_role_context_found": True, "forbidden_misclassification_violations": []},
            {"target_id": "XCL1_lymphotactin", "role_classification_passed": True, "positive_role_context_found": True, "forbidden_misclassification_violations": []},
        ],
    }
    gap = {
        "kind": "V16_PRESSURE_EVIDENCE_GAP_LOCK_v0",
        "gap_lock_status": "V16_PRESSURE_EVIDENCE_GAP_LOCKED",
        "claim_allowed": False,
        "new_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "target_specific_threshold_tuning_allowed": False,
        "target_gap_rows": [
            {"target_id": "p53_TAD_MDM2", "missing_for_evidence_test": ["isolated_TAD_disorder_or_clean_abstain_context"]},
            {"target_id": "KcsA", "missing_for_evidence_test": ["membrane_topology_annotation"]},
            {"target_id": "XCL1_lymphotactin", "missing_for_evidence_test": ["leakage_guard_preventing_mixed_state_fake_core"]},
        ],
    }
    return preflight, zero, gap


def test_v17_sprint_selects_single_p53_v18_target_without_claim(tmp_path: Path) -> None:
    cert = v17.build_sprint(*_certs(tmp_path))
    assert cert["sprint_lock_status"] == "V17_PRESSURE_EVIDENCE_SPRINT_LOCKED"
    assert cert["selected_V18_target"] == "p53_TAD_MDM2"
    assert cert["selected_V18_test"] == "p53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST"
    assert cert["ready_for_V18_targets"] == ["p53_TAD_MDM2"]
    assert cert["positive_folding_evidence_targets"] == []
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["fixed_residue_cutoff_used"] is False
    assert cert["native_metrics_used_for_selection"] is False


def test_v17_sprint_explains_why_kcsa_and_xcl1_are_deferred(tmp_path: Path) -> None:
    cert = v17.build_sprint(*_certs(tmp_path))
    rows = {row["target_id"]: row for row in cert["target_rows"]}
    assert rows["KcsA"]["ready_for_V18"] is False
    assert "membrane_topology_annotation" in rows["KcsA"]["missing_for_first_v18_test"]
    assert rows["XCL1_lymphotactin"]["ready_for_V18"] is False
    assert "leakage_guard_preventing_mixed_state_fake_core" in rows["XCL1_lymphotactin"]["missing_for_first_v18_test"]
    assert rows["p53_TAD_MDM2"]["ready_for_V18"] is True
    assert "interface_or_contact_evidence" in rows["p53_TAD_MDM2"]["available_evidence"]


def test_v17_sprint_writes_four_requested_artifacts(tmp_path: Path) -> None:
    cert = v17.build_sprint(*_certs(tmp_path))
    out_dir = tmp_path / "out"
    v17.write_outputs(out_dir, cert)
    assert (out_dir / "v17_pressure_evidence_manifest.json").exists()
    assert (out_dir / "v17_pressure_evidence_preflight.json").exists()
    assert (out_dir / "v17_zero_md_evidence_readout.json").exists()
    assert (out_dir / "v17_next_target_decision.json").exists()
    decision = json.loads((out_dir / "v17_next_target_decision.json").read_text())
    assert decision["selected_V18_target"] == "p53_TAD_MDM2"
