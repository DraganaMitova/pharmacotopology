#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST" / "v28_xcl1_state_condition_evidence_contrast_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V28 XCL1 STATE-CONDITION EVIDENCE CONTRAST TEST ===")
    for key in [
        "run_mode",
        "test_status",
        "target_id",
        "role_class",
        "mechanism_operator_tested",
        "source_v27_status",
        "source_v26_status",
        "source_v20_status",
        "condition_label_context_locked",
        "state_A_condition_evidence_found",
        "state_B_condition_evidence_found",
        "state_condition_contrast_preserved",
        "state_specific_buckets_preserved",
        "state_specific_couplings_present",
        "external_couplings_required_for_V28",
        "external_couplings_required_for_future_MD_or_contact_test",
        "positive_pressure_evidence_found",
        "positive_folding_evidence_found",
        "claim_allowed",
        "new_md_executed",
        "fixed_threshold_policy",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "mixed_state_pollution",
        "mixed_state_contact_pooling_used",
        "mixed_state_fake_core_selected",
        "single_fold_forcing",
        "single_fold_claim_made",
        "fold_switch_claim_made",
        "forbidden_misclassification_violations",
        "failed_checks",
        "available_evidence",
        "missing_evidence",
    ]:
        print(f"{key}: {data.get(key)}")

    print("\nState-condition contrast readout:")
    contrast = data.get("state_condition_contrast_readout", {})
    for state in ["state_A", "state_B"]:
        row = contrast.get(state, {})
        print(f"  {state}:")
        for key in ["label", "source_pdb", "condition_role_evidence_found", "role_bucket", "expected_context", "chain_ca_counts", "model_count", "interface_pairs_present"]:
            print(f"    {key}: {row.get(key)}")
    guard = contrast.get("guard", {})
    print("\nGuard:")
    for key in ["guard_policy", "state_A_state_B_contact_pooling_used", "mixed_state_fake_core_selected", "single_native_state_forced", "single_fold_claim_made", "fold_switch_claim_made", "selection_threshold_used", "state_specific_buckets_preserved"]:
        print(f"  {key}: {guard.get(key)}")

    print("\nNext decision:")
    dec = data.get("next_decision", {})
    for key in ["decision_status", "selected_next_panel", "reason", "new_MD_allowed", "new_MD_recommended", "new_MD_allowed_policy", "parallel_MD_paths_allowed", "claim_allowed", "required_before_any_MD"]:
        print(f"  {key}: {dec.get(key)}")


if __name__ == "__main__":
    main()
