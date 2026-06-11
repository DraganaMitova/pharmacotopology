from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v22_external_evidence_and_annotation_acquisition_panel_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v22_acquisition", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _v21_summary() -> dict:
    return {
        "summary_status": "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_LOCKED",
        "claim_allowed": False,
        "new_md_executed": False,
        "positive_pressure_evidence_targets": ["p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"],
        "positive_folding_evidence_targets": [],
        "target_rows": [
            {
                "target_id": "p53_TAD_MDM2",
                "positive_pressure_evidence_found": True,
                "available_evidence": [
                    "partner_bound_complex_context",
                    "interface_or_contact_evidence",
                    "leakage_guard_autonomous_fold_vs_partner_induced_fold",
                ],
                "missing_evidence": ["external_couplings_if_available"],
                "forbidden_misclassification_violations": [],
            },
            {
                "target_id": "KcsA",
                "positive_pressure_evidence_found": True,
                "available_evidence": ["KcsA_complex_coordinate_context", "leakage_guard_against_soluble_core_misread"],
                "missing_evidence": ["biological_tetramer_assembly_annotation_for_tetramer_claim"],
                "forbidden_misclassification_violations": [],
            },
            {
                "target_id": "XCL1_lymphotactin",
                "positive_pressure_evidence_found": True,
                "available_evidence": ["state_A_and_state_B_context_structures", "leakage_guard_preventing_mixed_state_fake_core"],
                "missing_evidence": ["condition_labels_if_available"],
                "forbidden_misclassification_violations": [],
            },
        ],
    }


def _write_required_structures(tmp_path: Path) -> None:
    paths = [
        tmp_path / "data" / "v16_pressure_targets" / "p53_TAD_MDM2" / "1YCR.pdb",
        tmp_path / "data" / "v16_pressure_targets" / "KcsA" / "1K4C.pdb",
        tmp_path / "data" / "v16_pressure_targets" / "XCL1_lymphotactin" / "2HDM.pdb",
        tmp_path / "data" / "v16_pressure_targets" / "XCL1_lymphotactin" / "2JP1.pdb",
    ]
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("HEADER TEST\nATOM      1  CA  ALA A   1       0.0   0.0   0.0\n", encoding="utf-8")


def _patch_structure_files(mod, tmp_path: Path) -> None:
    mod.LOCAL_STRUCTURE_FILES = {
        "p53_TAD_MDM2": [("p53_TAD_MDM2_complex_structure", tmp_path / "data" / "v16_pressure_targets" / "p53_TAD_MDM2" / "1YCR.pdb")],
        "KcsA": [("KcsA_Fab_complex_structure", tmp_path / "data" / "v16_pressure_targets" / "KcsA" / "1K4C.pdb")],
        "XCL1_lymphotactin": [
            ("XCL1_state_A_chemokine_like_structure", tmp_path / "data" / "v16_pressure_targets" / "XCL1_lymphotactin" / "2HDM.pdb"),
            ("XCL1_state_B_alternative_conformation_structure", tmp_path / "data" / "v16_pressure_targets" / "XCL1_lymphotactin" / "2JP1.pdb"),
        ],
    }


def test_v22_selects_p53_for_v23_without_folding_claim(tmp_path: Path) -> None:
    mod = _load_module()
    _write_required_structures(tmp_path)
    _patch_structure_files(mod, tmp_path)
    cert = mod.build_v22(_v21_summary(), tmp_path)
    assert cert["panel_status"] == "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL_LOCKED"
    assert cert["selected_V23_target"] == "p53_TAD_MDM2"
    assert cert["selected_V23_test"] == "p53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST"
    assert cert["ready_for_V23_targets"] == ["p53_TAD_MDM2"]
    assert cert["positive_folding_evidence_targets"] == []
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["fixed_residue_cutoff_used"] is False
    assert cert["native_metrics_used_for_selection"] is False


def test_v22_reports_couplings_missing_without_blocking_p53_first_v23(tmp_path: Path) -> None:
    mod = _load_module()
    _write_required_structures(tmp_path)
    _patch_structure_files(mod, tmp_path)
    cert = mod.build_v22(_v21_summary(), tmp_path)
    coupling_rows = {row["target_id"]: row for row in cert["coupling_availability_scan"]["targets"]}
    assert coupling_rows["p53_TAD_MDM2"]["coupling_status"] == "missing_or_not_yet_acquired"
    assert coupling_rows["p53_TAD_MDM2"]["coupling_required_for_selected_V23"] is False
    assert cert["selected_V23_target"] == "p53_TAD_MDM2"


def test_v22_deferred_targets_require_local_annotations_or_couplings(tmp_path: Path) -> None:
    mod = _load_module()
    _write_required_structures(tmp_path)
    _patch_structure_files(mod, tmp_path)
    cert = mod.build_v22(_v21_summary(), tmp_path)
    preflight = {row["target_id"]: row for row in cert["annotation_preflight"]["targets"]}
    assert preflight["KcsA"]["ready_for_V23"] is False
    assert "local_biological_assembly_or_OPM_membrane_annotation_file" in preflight["KcsA"]["missing_evidence"]
    assert preflight["XCL1_lymphotactin"]["ready_for_V23"] is False
    assert "local_state_condition_or_monomer_dimer_annotation_file" in preflight["XCL1_lymphotactin"]["missing_evidence"]


def test_v22_writer_outputs_required_four_core_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    _write_required_structures(tmp_path)
    _patch_structure_files(mod, tmp_path)
    cert = mod.build_v22(_v21_summary(), tmp_path)
    out = tmp_path / "out"
    paths = mod.write_outputs(out, cert)
    for key in ["manifest", "annotation_preflight", "coupling_availability_scan", "next_target_decision", "certificate", "report"]:
        assert paths[key].exists()
    written = json.loads((out / "v22_external_evidence_and_annotation_acquisition_panel_certificate.json").read_text(encoding="utf-8"))
    assert written["panel_status"] == "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL_LOCKED"
    assert written["selected_V23_target"] == "p53_TAD_MDM2"
