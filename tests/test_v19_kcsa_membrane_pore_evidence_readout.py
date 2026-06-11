from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v19_kcsa_membrane_pore_evidence_readout_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v19_kcsa", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _atom(serial: int, chain: str, res: int, resname: str, x: float, y: float, z: float) -> str:
    return (
        f"ATOM  {serial:5d}  CA  {resname:>3s} {chain}{res:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
    )


def _het_k(serial: int, chain: str, res: int, x: float, y: float, z: float) -> str:
    return (
        f"HETATM{serial:5d}  K     K {chain}{res:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           K\n"
    )


def _preflight(pdb: Path) -> dict:
    return {
        "data_preflight_status": "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT",
        "target_results": [
            {
                "target_id": "KcsA",
                "target_material_status": "ready_for_zero_md_role_readout",
                "expected_role_class": "membrane_pore_oligomer_object",
                "material_results": [
                    {
                        "material_id": "KcsA_Fab_complex_structure",
                        "material_status": "present_and_provenanced",
                        "path": str(pdb),
                    }
                ],
            }
        ],
    }


def _zero_md() -> dict:
    return {
        "role_transfer_status": "V16_ZERO_MD_ROLE_TRANSFER_READOUT_COMPLETED_CLAIM_DISABLED",
        "target_results": [
            {
                "target_id": "KcsA",
                "expected_role_class": "membrane_pore_oligomer_object",
                "role_classification_passed": True,
                "forbidden_misclassification_violations": [],
            }
        ],
    }


def _kcsa_like_pdb(tmp_path: Path) -> Path:
    pdb = tmp_path / "1K4C.pdb"
    lines = [
        "HEADER    MEMBRANE PROTEIN                         1K4C\n",
        "TITLE     POTASSIUM CHANNEL KCSA SELECTIVITY FILTER PORE CONTEXT\n",
        "KEYWDS    POTASSIUM CHANNEL, MEMBRANE, SELECTIVITY FILTER, PORE\n",
        "HELIX    1   1 ALA A    1  GLY A   15  1                                  15\n",
    ]
    serial = 1
    # Fab/context-like chains A and B plus channel-like chain C. Coordinates are close
    # enough to create report-only chain-interface contacts.
    for i in range(1, 9):
        lines.append(_atom(serial, "A", i, "ALA", float(i), 0.0, 0.0)); serial += 1
    for i in range(1, 9):
        lines.append(_atom(serial, "B", i, "ALA", float(i), 3.0, 0.0)); serial += 1
    for idx, resname in enumerate(["THR", "VAL", "GLY", "TYR", "GLY", "ASP", "ALA", "LEU"], start=1):
        lines.append(_atom(serial, "C", idx, resname, float(idx), 6.0, 0.0)); serial += 1
    lines.append(_het_k(serial, "C", 900, 3.0, 6.2, 0.0))
    pdb.write_text("".join(lines), encoding="utf-8")
    return pdb


def test_v19_kcsa_membrane_pore_evidence_passes_without_folding_claim(tmp_path: Path) -> None:
    mod = _load_module()
    pdb = _kcsa_like_pdb(tmp_path)
    cert = mod.build_evidence(_preflight(pdb), _zero_md(), pdb)
    assert cert["test_status"] == "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED"
    assert cert["positive_pressure_evidence_found"] is True
    assert cert["membrane_pore_role_evidence_found"] is True
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["membrane_md_executed"] is False
    assert cert["fixed_residue_cutoff_used"] is False
    assert cert["native_metrics_used_for_selection"] is False
    assert cert["pore_filter_annotation_present"] is True
    assert cert["transmembrane_role_present"] is True
    assert cert["oligomer_or_chain_interface_context_present"] is True
    assert cert["soluble_core_misclassification_avoided"] is True
    assert cert["forbidden_misclassification_violations"] == []
    assert "pore_selectivity_filter_core" in cert["role_buckets_assigned"]
    assert cert["pore_filter_readout"]["sequence_motif_probe"]["filter_motif_detected"] is True
    assert cert["chain_interface_readout"]["selection_threshold_used"] is False


def test_v19_kcsa_clean_abstains_when_pore_context_missing(tmp_path: Path) -> None:
    mod = _load_module()
    pdb = tmp_path / "not_kcsa.pdb"
    lines = ["HEADER    SOLUBLE TEST\n"]
    serial = 1
    for i in range(1, 6):
        lines.append(_atom(serial, "A", i, "ALA", float(i), 0.0, 0.0)); serial += 1
    pdb.write_text("".join(lines), encoding="utf-8")
    cert = mod.build_evidence(_preflight(pdb), _zero_md(), pdb)
    assert cert["test_status"] == "V19_KcsA_CLEAN_ABSTAIN_MISSING_MEMBRANE_PORE_OR_ASSEMBLY_EVIDENCE_CLAIM_DISABLED"
    assert cert["positive_pressure_evidence_found"] is False
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False
    assert cert["clean_abstain_valid"] is True
    assert "pore_selectivity_filter_annotation_or_diagnostic" in cert["missing_evidence"]


def test_v19_writer_emits_certificate_report_and_readouts(tmp_path: Path) -> None:
    mod = _load_module()
    cert = {
        "kind": "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_v0",
        "test_status": "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED",
        "positive_pressure_evidence_found": True,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "pore_filter_readout": {"sequence_motif_probe": {"filter_motif_detected": True}},
        "chain_interface_readout": {"chain_ca_counts": {"A": 8, "B": 8, "C": 8}},
        "forbidden_misclassification_violations": [],
    }
    out = tmp_path / "out"
    mod.write_outputs(out, cert)
    cpath = out / "v19_kcsa_membrane_pore_evidence_readout_certificate.json"
    rpath = out / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_REPORT.md"
    ppath = out / "v19_kcsa_pore_filter_readout.json"
    ipath = out / "v19_kcsa_chain_interface_readout.json"
    assert cpath.exists()
    assert rpath.exists()
    assert ppath.exists()
    assert ipath.exists()
    written = json.loads(cpath.read_text(encoding="utf-8"))
    assert written["claim_allowed"] is False
    assert written["positive_folding_evidence_found"] is False


def test_v19_reads_v16_zero_md_target_rows_schema_for_soluble_core_guard(tmp_path: Path) -> None:
    mod = _load_module()
    pdb = _kcsa_like_pdb(tmp_path)
    zero_md_target_rows_schema = {
        "role_transfer_status": "V16_ZERO_MD_ROLE_TRANSFER_READOUT_COMPLETED_CLAIM_DISABLED",
        "role_classification_passed_targets": ["KcsA"],
        "pressure_role_transfer_passed_targets": ["KcsA"],
        "target_rows": [
            {
                "target_id": "KcsA",
                "expected_role_class": "membrane_pore_oligomer_object",
                "positive_role_context_found": True,
                "forbidden_misclassification_violations": [],
                "target_transfer_status": "membrane_pore_context_detected_soluble_core_misclassification_avoided_no_tetramer_claim",
            }
        ],
    }
    cert = mod.build_evidence(_preflight(pdb), zero_md_target_rows_schema, pdb)
    assert cert["test_status"] == "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED"
    assert cert["positive_pressure_evidence_found"] is True
    assert cert["membrane_pore_role_evidence_found"] is True
    assert cert["soluble_core_misclassification_avoided"] is True
    assert "leakage_guard_against_soluble_core_misread" in cert["available_evidence"]
    assert "leakage_guard_against_soluble_core_misread" not in cert["missing_evidence"]
