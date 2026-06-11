#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL" / "v22_external_evidence_and_annotation_acquisition_panel_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V22 EXTERNAL EVIDENCE AND ANNOTATION ACQUISITION PANEL ===")
    for key in [
        "run_mode",
        "panel_status",
        "source_v21_summary_status",
        "claim_allowed",
        "evidence_claim_allowed",
        "new_md_executed",
        "fixed_threshold_policy",
        "target_specific_threshold_tuning_allowed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "positive_pressure_evidence_targets_from_v21",
        "positive_folding_evidence_targets",
        "ready_for_V23_targets",
        "selected_V23_target",
        "selected_V23_test",
        "blocked_or_deferred_targets",
        "annotation_source_refs_locked_targets",
        "targets_with_local_raw_annotation_files",
        "targets_with_couplings_or_msa",
        "panel_failed_checks",
    ]:
        print(f"{key}: {data.get(key)}")
    print("locked_interpretation:", data.get("locked_interpretation"))
    print("\nNext target decision:")
    decision = data.get("next_target_decision", {})
    for key in ["decision_status", "selected_V23_target", "selected_V23_test", "reason", "v23_first_run_mode", "parallel_md_paths_allowed"]:
        print(f"  {key}: {decision.get(key)}")
    print("\nAnnotation preflight rows:")
    for row in data.get("annotation_preflight", {}).get("targets", []):
        print(f"\n{row.get('target_id')}:")
        for key in [
            "role_class",
            "external_annotation_source_references_locked",
            "prior_pressure_evidence_locked",
            "prior_forbidden_misclassification_guard_clean",
            "ready_for_V23",
            "claim_allowed",
        ]:
            print(f"  {key}: {row.get(key)}")
        print("  raw_local_external_annotation_files:", row.get("raw_local_external_annotation_files"))
        print("  available_evidence:", row.get("available_evidence"))
        print("  missing_evidence:", row.get("missing_evidence"))
    print("\nCoupling availability:")
    for row in data.get("coupling_availability_scan", {}).get("targets", []):
        print(f"  {row.get('target_id')}: {row.get('coupling_status')} files={row.get('matching_files')}")


if __name__ == "__main__":
    main()
