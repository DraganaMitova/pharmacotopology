#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY" / "v21_pressure_evidence_panel_summary_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V21 PRESSURE EVIDENCE PANEL SUMMARY ===")
    for key in [
        "run_mode",
        "summary_status",
        "positive_pressure_evidence_targets",
        "positive_folding_evidence_targets",
        "pressure_classes_with_evidence",
        "pressure_evidence_count",
        "claim_allowed",
        "evidence_claim_allowed",
        "new_md_executed",
        "membrane_md_executed",
        "fixed_threshold_policy",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "md_readiness_decision",
        "next_recommended_panel",
        "summary_failed_checks",
    ]:
        print(f"{key}: {data.get(key)}")
    print("locked_interpretation:", data.get("locked_interpretation"))
    print("global_missing_evidence:", data.get("global_missing_evidence"))
    print("\nTarget rows:")
    for row in data.get("target_rows", []):
        print(f"\n{row.get('target_id')}:")
        for key in [
            "pressure_class",
            "role_class",
            "pressure_evidence_status",
            "positive_pressure_evidence_found",
            "positive_folding_evidence_found",
            "claim_allowed",
            "new_md_executed",
            "selected_core_or_clean_abstain",
            "evidence_summary",
            "md_readiness",
        ]:
            print(f"  {key}: {row.get(key)}")
        print("  available_evidence:", row.get("available_evidence"))
        print("  missing_evidence:", row.get("missing_evidence"))
        print("  forbidden_misclassification_violations:", row.get("forbidden_misclassification_violations"))


if __name__ == "__main__":
    main()
