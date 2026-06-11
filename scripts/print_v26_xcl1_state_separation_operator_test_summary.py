#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST" / "v26_xcl1_state_separation_operator_test_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V26 XCL1 STATE-SEPARATION OPERATOR TEST ===")
    for key in [
        "run_mode",
        "test_status",
        "target_id",
        "role_class",
        "mechanism_operator_tested",
        "source_v25_status",
        "source_v20_status",
        "selected_by_v25",
        "state_separation_operator_passed",
        "mechanism_operator_evidence_found",
        "positive_pressure_evidence_found",
        "positive_folding_evidence_found",
        "claim_allowed",
        "new_md_executed",
        "fixed_threshold_policy",
        "target_specific_threshold_tuning_allowed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "failed_checks",
    ]:
        print(f"{key}: {data.get(key)}")
    op = data.get("operator_readout", {})
    print("\nOperator readout:")
    for key in [
        "state_A_detected",
        "state_B_detected",
        "state_specific_role_evidence_found",
        "state_specific_buckets_preserved",
        "mixed_state_fake_core_selected",
        "mixed_state_pollution",
        "mixed_state_contact_pooling_used",
        "single_fold_forcing",
        "single_fold_claim_made",
        "fold_switch_claim_made",
        "state_separation_guard_passed",
        "selection_threshold_used",
        "available_evidence",
        "missing_evidence",
        "operator_policy",
    ]:
        print(f"  {key}: {op.get(key)}")
    print("\nNext decision:")
    dec = data.get("next_mechanism_decision", {})
    for key in [
        "decision_status",
        "selected_next_panel",
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
