#!/usr/bin/env python3
from __future__ import annotations
import os
from pathlib import Path

ROOT = Path(os.environ.get("REPO_ROOT", Path.cwd())).resolve()
BASE = ROOT / "first_contact_clean_pharmacotopology_layer_run"
needed = [
    "V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCKED/v15_dynamic_role_grammar_panel_locked_certificate.json",
    "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCK/v16_target_manifest_and_role_expectation_lock_certificate.json",
    "V17_PRESSURE_EVIDENCE_SPRINT_LOCK/v17_pressure_evidence_sprint_lock_certificate.json",
    "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCK/v18_p53_partner_induced_evidence_lock_certificate.json",
    "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT/v18b_p53_isolated_tad_abstain_context_certificate.json",
    "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCK/v19_kcsa_membrane_pore_evidence_lock_certificate.json",
    "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT/v20_xcl1_state_specific_evidence_readout_certificate.json",
    "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY/v21_pressure_evidence_panel_summary_certificate.json",
    "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL/v22_external_evidence_and_annotation_acquisition_panel_certificate.json",
    "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST/v23_p53_tad_mdm2_external_evidence_contrast_certificate.json",
    "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION/v24_kcsa_external_annotation_and_assembly_acquisition_certificate.json",
    "V25_FAST_MECHANISM_EVIDENCE_SPRINT/v25_fast_mechanism_evidence_sprint_certificate.json",
    "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST/v26_xcl1_state_separation_operator_test_certificate.json",
    "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION/v27_xcl1_condition_and_coupling_evidence_acquisition_certificate.json",
    "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST/v28_xcl1_state_condition_evidence_contrast_certificate.json",
    "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_DECISION/v29_mechanism_operator_panel_summary_and_md_readiness_certificate.json",
    "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT/v30_external_constraint_and_coupling_acquisition_sprint_certificate.json",
]
missing=[]
for rel in needed:
    p=BASE/rel
    ok=p.exists()
    print(("OK      " if ok else "MISSING ")+str(p))
    if not ok: missing.append(str(p))
print(f"missing_count = {len(missing)}")
