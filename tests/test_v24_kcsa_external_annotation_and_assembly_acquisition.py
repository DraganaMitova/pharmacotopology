from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v24_kcsa_external_annotation_and_assembly_acquisition_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v24_kcsa", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _v19_readout() -> dict:
    return {
        "test_status": "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED",
        "positive_pressure_evidence_found": True,
        "positive_folding_evidence_found": False,
        "membrane_pore_role_evidence_found": True,
        "claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "pore_filter_annotation_present": True,
        "transmembrane_role_present": True,
        "oligomer_or_chain_interface_context_present": True,
        "membrane_environment_context_present": True,
        "soluble_core_misclassification_avoided": True,
        "whole_channel_fold_claim_made": False,
        "tetramer_claim_made": False,
        "forbidden_misclassification_violations": [],
        "pore_filter_readout": {
            "sequence_motif_probe": {
                "strong_TVGYG_motif_detected": True,
                "selection_threshold_used": False,
            },
            "potassium_ion_probe": {
                "ion_filter_diagnostic_shell_present": True,
                "potassium_ion_count": 7,
                "selection_threshold_used": False,
            },
        },
        "transmembrane_helix_readout": {
            "transmembrane_helix_scaffold_context_present": True,
            "selection_threshold_used": False,
        },
        "chain_interface_readout": {
            "interface_pairs_present": ["A-B", "A-C", "B-C"],
            "chain_ca_counts": {"A": 219, "B": 212, "C": 103},
            "selection_threshold_used": False,
        },
    }


def _v19_lock() -> dict:
    return {
        "lock_status": "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCKED",
        "positive_pressure_evidence_found": True,
        "membrane_pore_role_evidence_found": True,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "soluble_core_misclassification_avoided": True,
        "whole_channel_fold_claim_made": False,
        "tetramer_claim_made": False,
        "forbidden_misclassification_violations": [],
    }


def _v22() -> dict:
    return {
        "manifest": {
            "targets": [
                {
                    "target_id": "KcsA",
                    "external_annotation_source_references": [
                        {
                            "source_id": "RCSB_1K4C",
                            "kind": "KcsA_Fab_ion_filter_structure_reference",
                            "url": "https://www.rcsb.org/structure/1K4C",
                            "use_boundary": "pore_filter_ion_context_only_no_whole_channel_fold_claim",
                        },
                        {
                            "source_id": "RCSB_1BL8",
                            "kind": "KcsA_integral_membrane_potassium_channel_reference",
                            "url": "https://www.rcsb.org/structure/1BL8",
                            "use_boundary": "membrane_pore_context_only_no_tetramer_or_MD_claim",
                        },
                        {
                            "source_id": "OPM_KcsA",
                            "kind": "transmembrane_topology_reference",
                            "url": "https://opm.phar.umich.edu/proteins/59",
                            "use_boundary": "topology_annotation_reference_only_no_membrane_MD_claim",
                        },
                    ],
                }
            ]
        },
        "coupling_availability_scan": {
            "targets": [
                {"target_id": "KcsA", "matching_files": [], "coupling_status": "missing_or_not_yet_acquired"}
            ]
        },
    }


def test_v24_passes_kcsa_external_annotation_layer_without_folding_claim() -> None:
    mod = _load_module()
    cert = mod.build_v24(_v19_readout(), _v19_lock(), _v22())
    assert cert["test_status"] == "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION_PASSED_CLAIM_DISABLED"
    assert cert["positive_pressure_evidence_found"] is True
    assert cert["positive_folding_evidence_found"] is False
    assert cert["membrane_pore_role_evidence_found"] is True
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["membrane_md_executed"] is False
    assert cert["pore_filter_annotation_present"] is True
    assert cert["TVGYG_motif_detected"] is True
    assert cert["K_plus_ion_diagnostic_shell_present"] is True
    assert cert["transmembrane_role_present"] is True
    assert cert["soluble_core_misclassification_avoided"] is True


def test_v24_locks_assembly_context_but_forbids_tetramer_claim() -> None:
    mod = _load_module()
    cert = mod.build_v24(_v19_readout(), _v19_lock(), _v22())
    assert cert["tetramer_context_annotation_present"] is True
    assert cert["biological_assembly_context_locked"] is True
    assert cert["tetramer_claim_allowed"] is False
    assert cert["tetramer_claim_made"] is False
    assert cert["whole_channel_fold_claim_allowed"] is False
    assert cert["whole_channel_fold_claim_made"] is False
    assert "expanded_biological_assembly_coordinate_or_author_defined_assembly_model_for_tetramer_claim" in cert["missing_evidence"]


def test_v24_does_not_require_couplings_for_acquisition_but_marks_future_need() -> None:
    mod = _load_module()
    cert = mod.build_v24(_v19_readout(), _v19_lock(), _v22())
    assert cert["external_couplings_present"] is False
    assert cert["external_couplings_required_for_V24"] is False
    assert cert["external_couplings_required_for_future_MD_or_contact_test"] is True
    assert cert["external_coupling_availability"]["coupling_policy"] == "do_not_synthesize_couplings_missing_is_not_v24_failure"


def test_v24_blocks_if_v19_lock_missing() -> None:
    mod = _load_module()
    v19_lock = _v19_lock()
    v19_lock["lock_status"] = "V19_LOCK_BLOCKED"
    cert = mod.build_v24(_v19_readout(), v19_lock, _v22())
    assert cert["test_status"] == "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION_BLOCKED_CLAIM_DISABLED"
    assert "source_v19_membrane_pore_lock_present" in cert["v24_failed_checks"]
    assert cert["claim_allowed"] is False


def test_v24_writes_all_requested_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_v24(_v19_readout(), _v19_lock(), _v22())
    paths = mod.write_outputs(tmp_path, cert)
    expected = {
        "certificate",
        "manifest",
        "assembly_preflight",
        "pore_filter_annotation",
        "membrane_topology_annotation",
        "external_coupling_availability",
        "next_decision",
        "report",
    }
    assert expected.issubset(paths)
    for key in expected:
        assert paths[key].exists()
    written = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert written["test_status"] == "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION_PASSED_CLAIM_DISABLED"
    assert written["positive_folding_evidence_found"] is False
