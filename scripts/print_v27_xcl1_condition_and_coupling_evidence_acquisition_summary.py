#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION" / "v27_xcl1_condition_and_coupling_evidence_acquisition_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V27 XCL1 CONDITION AND COUPLING EVIDENCE ACQUISITION ===")
    for key in [
        "run_mode",
        "test_status",
        "target_id",
        "role_class",
        "mechanism_operator",
        "source_v26_status",
        "condition_label_context_locked",
        "state_A_context_label_locked",
        "state_B_context_label_locked",
        "state_labels_locked",
        "monomer_dimer_context_locked",
        "mixed_state_leakage_guard_preserved",
        "state_specific_couplings_present",
        "state_specific_external_constraints_present",
        "external_couplings_required_for_V27",
        "external_couplings_required_for_future_MD_or_contact_test",
        "positive_pressure_evidence_found",
        "positive_folding_evidence_found",
        "claim_allowed",
        "new_md_executed",
        "fixed_threshold_policy",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "mixed_state_contact_pooling_used",
        "mixed_state_pollution",
        "single_fold_forcing",
        "single_fold_claim_made",
        "fold_switch_claim_made",
        "forbidden_misclassification_violations",
        "failed_checks",
        "available_evidence",
        "missing_evidence",
    ]:
        print(f"{key}: {data.get(key)}")

    print("\nCondition manifest:")
    manifest = data.get("condition_label_manifest", {})
    for key in [
        "state_A_label",
        "state_B_label",
        "state_label_source_policy",
        "condition_label_context_locked",
        "raw_local_condition_annotation_files_present",
        "raw_local_condition_annotation_files",
    ]:
        print(f"  {key}: {manifest.get(key)}")

    print("\nCoupling availability:")
    coupling = data.get("state_specific_coupling_availability", {})
    for key in [
        "state_specific_couplings_present",
        "state_specific_external_constraints_present",
        "external_coupling_or_constraint_files",
        "external_couplings_required_for_V27",
        "external_couplings_required_for_future_MD_or_contact_test",
        "coupling_policy",
    ]:
        print(f"  {key}: {coupling.get(key)}")

    print("\nNext decision:")
    dec = data.get("next_decision", {})
    for key in [
        "decision_status",
        "selected_next_panel",
        "selected_V28_target",
        "selected_V28_test",
        "reason",
        "new_MD_allowed",
        "new_MD_recommended",
        "new_MD_allowed_policy",
        "parallel_MD_paths_allowed",
        "claim_allowed",
        "required_before_any_MD",
    ]:
        print(f"  {key}: {dec.get(key)}")


if __name__ == "__main__":
    main()
