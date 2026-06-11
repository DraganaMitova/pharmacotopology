from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v21_pressure_evidence_panel_summary_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v21_summary", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _v18_lock() -> dict:
    return {
        "lock_status": "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCKED",
        "positive_pressure_evidence_found": True,
        "partner_induced_role_evidence_found": True,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "available_evidence": ["partner_bound_complex_context", "interface_or_contact_evidence"],
        "missing_evidence": ["external_couplings_if_available"],
        "forbidden_misclassification_violations": [],
    }


def _v18b() -> dict:
    return {
        "test_status": "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT_PASSED_CLAIM_DISABLED",
        "bound_complex_partner_induced_evidence_preserved": True,
        "positive_pressure_evidence_found": True,
        "positive_folding_evidence_found": False,
        "isolated_TAD_autonomous_core_selected": False,
        "isolated_TAD_fold_claim_made": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "selected_core_or_clean_abstain": "isolated_clean_abstain_bound_partner_induced_evidence_preserved",
        "available_evidence": ["isolated_context_clean_abstain_guard"],
        "missing_evidence": ["isolated_TAD_disorder_prior_or_external_annotation"],
        "forbidden_misclassification_violations": [],
    }


def _v19_lock() -> dict:
    return {
        "lock_status": "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCKED",
        "positive_pressure_evidence_found": True,
        "membrane_pore_role_evidence_found": True,
        "positive_folding_evidence_found": False,
        "soluble_core_misclassification_avoided": True,
        "whole_channel_fold_claim_made": False,
        "tetramer_claim_made": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "available_evidence": ["pore_selectivity_filter_annotation_or_diagnostic"],
        "missing_evidence": ["biological_tetramer_assembly_annotation_for_tetramer_claim"],
        "forbidden_misclassification_violations": [],
    }


def _v20() -> dict:
    return {
        "test_status": "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED",
        "positive_pressure_evidence_found": True,
        "state_specific_role_evidence_found": True,
        "state_A_role_evidence_found": True,
        "state_B_role_evidence_found": True,
        "positive_folding_evidence_found": False,
        "mixed_state_pollution": False,
        "single_fold_claim_made": False,
        "fold_switch_claim_made": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "selected_core_or_clean_abstain": "state_specific_role_separation_evidence_found_no_mixed_state_fake_core",
        "available_evidence": ["state_A_and_state_B_context_structures", "leakage_guard_preventing_mixed_state_fake_core"],
        "missing_evidence": ["condition_labels_if_available"],
        "forbidden_misclassification_violations": [],
    }


def test_v21_summary_locks_three_pressure_evidence_wins_without_folding_claim() -> None:
    mod = _load_module()
    cert = mod.build_summary(_v18_lock(), _v18b(), _v19_lock(), _v20())
    assert cert["summary_status"] == "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_LOCKED"
    assert cert["pressure_evidence_count"] == 3
    assert cert["positive_pressure_evidence_targets"] == ["p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"]
    assert cert["positive_folding_evidence_targets"] == []
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["fixed_residue_cutoff_used"] is False
    assert cert["native_metrics_used_for_selection"] is False
    assert cert["summary_failed_checks"] == []
    assert cert["md_readiness_decision"] == "no_new_MD_recommended_yet_external_evidence_and_annotations_first"


def test_v21_summary_blocks_if_xcl1_mixed_state_pollution_occurs() -> None:
    mod = _load_module()
    v20 = _v20()
    v20["mixed_state_pollution"] = True
    cert = mod.build_summary(_v18_lock(), _v18b(), _v19_lock(), v20)
    assert cert["summary_status"] == "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_BLOCKED"
    assert "XCL1_lymphotactin_pressure_evidence_passed" in cert["summary_failed_checks"]
    assert cert["positive_folding_evidence_targets"] == []
    assert cert["claim_allowed"] is False


def test_v21_summary_blocks_if_kcsa_md_claim_is_made() -> None:
    mod = _load_module()
    v19 = _v19_lock()
    v19["membrane_md_executed"] = True
    cert = mod.build_summary(_v18_lock(), _v18b(), v19, _v20())
    assert cert["summary_status"] == "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_BLOCKED"
    assert "KcsA_pressure_evidence_passed" in cert["summary_failed_checks"]
    assert cert["new_md_executed"] is False


def test_v21_writer_outputs_certificate_report_and_table(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_summary(_v18_lock(), _v18b(), _v19_lock(), _v20())
    out = tmp_path / "out"
    mod.write_outputs(out, cert)
    assert (out / "v21_pressure_evidence_panel_summary_certificate.json").exists()
    assert (out / "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_REPORT.md").exists()
    assert (out / "v21_pressure_evidence_panel_table.csv").exists()
    written = json.loads((out / "v21_pressure_evidence_panel_summary_certificate.json").read_text(encoding="utf-8"))
    assert written["summary_status"] == "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_LOCKED"
    assert written["positive_folding_evidence_targets"] == []
