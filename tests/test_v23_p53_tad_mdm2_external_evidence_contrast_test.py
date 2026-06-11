from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v23_p53_tad_mdm2_external_evidence_contrast_test_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v23_p53_external", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _v22() -> dict:
    return {
        "panel_status": "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL_LOCKED",
        "selected_V23_target": "p53_TAD_MDM2",
        "selected_V23_test": "p53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST",
        "claim_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_targets": [],
        "manifest": {
            "targets": [
                {
                    "target_id": "p53_TAD_MDM2",
                    "external_annotation_source_references": [
                        {
                            "source_id": "DisProt_DP00086",
                            "kind": "disorder_annotation_reference",
                            "url": "https://disprot.org/DP00086",
                            "use_boundary": "reference_lock_for_isolated_TAD_disorder_context_not_raw_API_claim",
                        },
                        {
                            "source_id": "RCSB_1YCR",
                            "kind": "partner_bound_interface_structure_reference",
                            "url": "https://www.rcsb.org/structure/1YCR",
                            "use_boundary": "bound_p53_TAD_MDM2_complex_context_only_no_autonomous_fold_claim",
                        },
                    ],
                }
            ]
        },
        "annotation_preflight": {
            "targets": [
                {
                    "target_id": "p53_TAD_MDM2",
                    "external_annotation_source_references_locked": True,
                    "raw_local_external_annotation_files": [],
                    "ready_for_V23": True,
                }
            ]
        },
        "coupling_availability_scan": {
            "targets": [
                {
                    "target_id": "p53_TAD_MDM2",
                    "coupling_status": "missing_or_not_yet_acquired",
                    "matching_files": [],
                }
            ]
        },
    }


def _v18() -> dict:
    return {
        "test_status": "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_PASSED_CLAIM_DISABLED",
        "positive_pressure_evidence_found": True,
        "partner_induced_role_evidence_found": True,
        "positive_folding_evidence_found": False,
        "new_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "interface_contact_readout": {
            "chain_pair": "A-B",
            "chain_ca_counts": {"A": 85, "B": 13},
            "min_ca_distance": 4.884,
            "multi_radius_contact_counts": {"6.0": 7, "8.0": 28, "10.0": 96},
            "contact_probe_policy": "multi_radius_report_only_no_single_fixed_selection_threshold",
            "interface_contact_evidence_present": True,
        },
        "partner_bound_helix_proxy_readout": {
            "small_chain": "B",
            "small_chain_ca_count": 13,
            "helix_proxy_policy": "report_only_partner_bound_segment_geometry_no_autonomous_fold_claim",
        },
    }


def _v18b() -> dict:
    return {
        "test_status": "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT_PASSED_CLAIM_DISABLED",
        "isolated_TAD_autonomous_core_selected": False,
        "isolated_TAD_fold_claim_made": False,
        "bound_complex_partner_induced_evidence_preserved": True,
        "new_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
    }


def test_v23_passes_external_contrast_without_folding_claim() -> None:
    mod = _load_module()
    cert = mod.build_v23(_v22(), _v18(), _v18b())
    assert cert["test_status"] == "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST_PASSED_CLAIM_DISABLED"
    assert cert["positive_pressure_evidence_found"] is True
    assert cert["positive_external_annotation_evidence_found"] is True
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["isolated_TAD_clean_abstain_preserved"] is True
    assert cert["partner_bound_interface_evidence_preserved"] is True
    assert cert["forbidden_misclassification_violations"] == []


def test_v23_does_not_require_couplings_for_first_external_contrast() -> None:
    mod = _load_module()
    cert = mod.build_v23(_v22(), _v18(), _v18b())
    assert cert["external_couplings_or_msa_files_present"] is False
    assert cert["external_couplings_required_for_this_test"] is False
    assert "external_couplings_if_available" in cert["missing_evidence"]
    assert cert["failed_checks"] == []


def test_v23_blocks_if_external_disorder_reference_missing() -> None:
    mod = _load_module()
    v22 = _v22()
    v22["manifest"]["targets"][0]["external_annotation_source_references"] = [
        {
            "source_id": "RCSB_1YCR",
            "kind": "partner_bound_interface_structure_reference",
            "url": "https://www.rcsb.org/structure/1YCR",
            "use_boundary": "bound_p53_TAD_MDM2_complex_context_only_no_autonomous_fold_claim",
        }
    ]
    cert = mod.build_v23(v22, _v18(), _v18b())
    assert cert["test_status"] == "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST_BLOCKED_CLAIM_DISABLED"
    assert "isolated_TAD_disorder_reference_locked" in cert["failed_checks"]
    assert cert["claim_allowed"] is False


def test_v23_writer_outputs_certificate_and_report(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_v23(_v22(), _v18(), _v18b())
    paths = mod.write_outputs(tmp_path, cert)
    assert paths["certificate"].exists()
    assert paths["report"].exists()
    written = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert written["test_status"] == "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST_PASSED_CLAIM_DISABLED"
    assert written["positive_folding_evidence_found"] is False
