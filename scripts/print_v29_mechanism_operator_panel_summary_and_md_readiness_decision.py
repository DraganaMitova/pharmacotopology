#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_DECISION" / "v29_mechanism_operator_panel_summary_and_md_readiness_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V29 MECHANISM OPERATOR PANEL SUMMARY AND MD READINESS DECISION ===")
    for key in [
        "run_mode",
        "summary_status",
        "source_v25_status",
        "source_v28_status",
        "claim_allowed",
        "new_md_executed",
        "membrane_md_executed",
        "fixed_threshold_policy",
        "target_specific_threshold_tuning_allowed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "positive_pressure_evidence_targets",
        "positive_folding_evidence_targets",
        "operator_panel_target_count",
        "repeated_mechanism_operators",
        "core_operators_across_all_targets",
        "universal_core_operator_present",
        "xcl1_state_condition_operator_preserved",
        "new_MD_allowed",
        "new_MD_recommended",
        "md_readiness_decision_status",
        "selected_next_panel",
        "summary_failed_checks",
    ]:
        print(f"{key}: {data.get(key)}")

    print("\nOperator panel:")
    panel = data.get("mechanism_operator_panel_summary", {})
    print("  operator_counts:", panel.get("operator_counts"))
    print("  operator_interpretation:", panel.get("operator_interpretation"))
    for row in panel.get("targets", []):
        print(f"  {row.get('target')}:")
        print(f"    pressure_type: {row.get('pressure_type')}")
        print(f"    role_detected: {row.get('role_detected')}")
        print(f"    evidence_type: {row.get('evidence_type')}")
        print(f"    selected_signal: {row.get('selected_signal')}")
        print(f"    missing_next_evidence: {row.get('missing_next_evidence')}")
        print(f"    mechanism_operator_used: {row.get('mechanism_operator_used')}")
        print(f"    claim_allowed: {row.get('claim_allowed')}")

    print("\nMD readiness rows:")
    for row in data.get("md_readiness_rows", []):
        print(f"  {row.get('target')}: {row.get('md_readiness')}")
        print(f"    new_MD_allowed: {row.get('new_MD_allowed')}")
        print(f"    new_MD_recommended: {row.get('new_MD_recommended')}")
        print(f"    blockers: {row.get('md_blockers')}")

    print("\nMD readiness decision:")
    dec = data.get("md_readiness_decision", {})
    for key in [
        "decision_status",
        "md_ready_targets",
        "selected_next_panel",
        "selected_next_focus",
        "reason",
        "new_MD_allowed",
        "new_MD_recommended",
        "membrane_MD_recommended",
        "parallel_MD_paths_allowed",
        "required_before_any_MD",
        "claim_allowed",
    ]:
        print(f"  {key}: {dec.get(key)}")


if __name__ == "__main__":
    main()
