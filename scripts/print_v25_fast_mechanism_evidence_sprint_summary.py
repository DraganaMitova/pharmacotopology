#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V25_FAST_MECHANISM_EVIDENCE_SPRINT" / "v25_fast_mechanism_evidence_sprint_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V25 FAST MECHANISM EVIDENCE SPRINT ===")
    for key in [
        "run_mode",
        "sprint_status",
        "claim_allowed",
        "new_md_executed",
        "membrane_md_executed",
        "fixed_threshold_policy",
        "target_specific_threshold_tuning_allowed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "positive_pressure_evidence_targets",
        "positive_folding_evidence_targets",
        "fast_sprint_failed_checks",
    ]:
        print(f"{key}: {data.get(key)}")

    print("\nKcsA coupling/interface readout:")
    kcsa = data.get("kcsa_coupling_interface_readout", {})
    for key in [
        "readout_status",
        "positive_pressure_evidence_found",
        "kcsa_coupling_available",
        "interface_evidence_present",
        "pore_filter_coupling_support",
        "tetramer_claim_allowed",
        "whole_channel_fold_claim_allowed",
        "positive_folding_evidence_found",
        "claim_allowed",
        "new_md_executed",
        "membrane_md_executed",
        "missing_evidence",
    ]:
        print(f"  {key}: {kcsa.get(key)}")

    print("\nXCL1 state-specific readout:")
    xcl1 = data.get("xcl1_state_specific_readout", {})
    for key in [
        "readout_status",
        "state_A_detected",
        "state_B_detected",
        "state_specific_role_evidence_found",
        "mixed_state_pollution",
        "single_fold_forcing",
        "positive_folding_evidence_found",
        "claim_allowed",
        "new_md_executed",
        "missing_evidence",
    ]:
        print(f"  {key}: {xcl1.get(key)}")

    print("\nMechanism operators:")
    table = data.get("unified_mechanism_operator_table", {})
    print("  operator_counts:", table.get("operator_counts"))
    for row in table.get("targets", []):
        print(f"  {row.get('target')}:")
        print(f"    pressure_type: {row.get('pressure_type')}")
        print(f"    role_detected: {row.get('role_detected')}")
        print(f"    evidence_type: {row.get('evidence_type')}")
        print(f"    selected_signal: {row.get('selected_signal')}")
        print(f"    abstain_reason_if_any: {row.get('abstain_reason_if_any')}")
        print(f"    missing_next_evidence: {row.get('missing_next_evidence')}")
        print(f"    mechanism_operator_used: {row.get('mechanism_operator_used')}")
        print(f"    claim_allowed: {row.get('claim_allowed')}")

    print("\nNext V26 decision:")
    dec = data.get("next_mechanism_test_decision", {})
    for key in [
        "decision_status",
        "selected_V26_target",
        "selected_V26_test",
        "reason",
        "new_MD_allowed",
        "new_MD_recommended",
        "new_MD_allowed_policy",
        "parallel_MD_paths_allowed",
        "claim_allowed",
    ]:
        print(f"  {key}: {dec.get(key)}")


if __name__ == "__main__":
    main()
