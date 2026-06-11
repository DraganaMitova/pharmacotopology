from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v18_p53_tad_mdm2_partner_induced_evidence_test_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v18_p53", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _atom(serial: int, chain: str, res: int, x: float, y: float, z: float) -> str:
    return (
        f"ATOM  {serial:5d}  CA  ALA {chain}{res:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
    )


def test_v18_p53_partner_induced_evidence_passes_without_folding_claim(tmp_path: Path) -> None:
    mod = _load_module()
    pdb = tmp_path / "1YCR.pdb"
    lines = []
    # Large partner chain A and small p53-like chain B, close enough for interface probes.
    for i in range(10):
        lines.append(_atom(i + 1, "A", i + 1, float(i), 0.0, 0.0))
    for i in range(8):
        lines.append(_atom(100 + i, "B", i + 1, float(i), 3.0, 0.0))
    pdb.write_text("".join(lines), encoding="utf-8")

    v17 = {
        "sprint_lock_status": "V17_PRESSURE_EVIDENCE_SPRINT_LOCKED",
        "selected_V18_target": "p53_TAD_MDM2",
        "selected_V18_test": "p53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST",
        "claim_allowed": False,
    }
    preflight = {
        "data_preflight_status": "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT",
        "target_results": [
            {
                "target_id": "p53_TAD_MDM2",
                "target_material_status": "ready_for_zero_md_role_readout",
                "material_results": [
                    {
                        "material_id": "p53_TAD_MDM2_complex_structure",
                        "material_status": "present_and_provenanced",
                    }
                ],
            }
        ],
    }
    cert = mod.build_evidence(v17, preflight, pdb)
    assert cert["test_status"] == "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_PASSED_CLAIM_DISABLED"
    assert cert["positive_pressure_evidence_found"] is True
    assert cert["partner_induced_role_evidence_found"] is True
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["fixed_residue_cutoff_used"] is False
    assert cert["native_metrics_used_for_selection"] is False
    assert cert["forbidden_misclassification_violations"] == []
    assert cert["isolated_TAD_autonomous_fold_status"] == "clean_abstain_no_isolated_TAD_material_no_autonomous_fold_claim"
    assert cert["interface_contact_readout"]["interface_contact_evidence_present"] is True
    assert cert["interface_contact_readout"]["selection_threshold_used"] is False


def test_v18_p53_writer_emits_certificate_report_and_interface_readout(tmp_path: Path) -> None:
    mod = _load_module()
    cert = {
        "kind": "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_v0",
        "test_status": "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_PASSED_CLAIM_DISABLED",
        "positive_pressure_evidence_found": True,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "interface_contact_readout": {"chain_pair": "A-B", "multi_radius_contact_counts": {"6.0": 1}},
        "partner_bound_helix_proxy_readout": {"small_chain": "B"},
        "isolated_TAD_autonomous_fold_status": "clean_abstain_no_isolated_TAD_material_no_autonomous_fold_claim",
        "forbidden_misclassification_violations": [],
    }
    out = tmp_path / "out"
    mod.write_outputs(out, cert)
    cpath = out / "v18_p53_tad_mdm2_partner_induced_evidence_certificate.json"
    rpath = out / "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_REPORT.md"
    ipath = out / "v18_p53_tad_mdm2_interface_contact_readout.json"
    assert cpath.exists()
    assert rpath.exists()
    assert ipath.exists()
    written = json.loads(cpath.read_text(encoding="utf-8"))
    assert written["claim_allowed"] is False
    assert written["positive_folding_evidence_found"] is False
