#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERT = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "V16_PRESSURE_EVIDENCE_GAP_LOCK"
    / "v16_pressure_evidence_gap_lock_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Print V16 pressure evidence gap lock summary.")
    parser.add_argument("--cert", default=str(DEFAULT_CERT))
    args = parser.parse_args()
    path = Path(args.cert)
    if not path.exists():
        raise SystemExit(f"missing V16 pressure evidence gap lock certificate: {path}")
    o = json.loads(path.read_text(encoding="utf-8"))

    print("=== V16 PRESSURE EVIDENCE GAP LOCK ===")
    for key in [
        "run_mode",
        "gap_lock_status",
        "source_zero_md_status",
        "source_v16_lock_status",
        "global_missing_layer",
        "next_panel",
        "claim_allowed",
        "new_md_executed",
        "fixed_threshold_policy",
        "threshold_policy",
        "target_specific_threshold_tuning_allowed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
    ]:
        print(f"{key}: {o.get(key)}")
    print("gap_lock_failed_checks:", o.get("gap_lock_failed_checks"))
    print("role_classification_passed_targets:", o.get("role_classification_passed_targets"))
    print("pressure_role_transfer_passed_targets:", o.get("pressure_role_transfer_passed_targets"))
    print("positive_folding_evidence_targets:", o.get("positive_folding_evidence_targets"))
    print("evidence_claim_allowed_targets:", o.get("evidence_claim_allowed_targets"))
    print()
    print("Locked interpretation:")
    print(o.get("locked_interpretation"))
    print()
    print("Target gaps:")
    for row in o.get("target_gap_rows", []):
        print(f"\n{row.get('target_id')}:")
        print("  role_class:", row.get("role_class"))
        print("  role_transfer_status:", row.get("role_transfer_status"))
        print("  role_classification_passed:", row.get("role_classification_passed"))
        print("  positive_folding_evidence_found:", row.get("positive_folding_evidence_found"))
        print("  forbidden_misclassification_avoided:", row.get("forbidden_misclassification_avoided"))
        print("  forbidden_misclassification_guard:", row.get("forbidden_misclassification_guard"))
        print("  missing_for_evidence_test:")
        for item in row.get("missing_for_evidence_test", []):
            print("    -", item)
        print("  evidence_claim_forbidden_until:")
        for item in row.get("evidence_claim_forbidden_until", []):
            print("    -", item)
        print("  next_evidence_test:", row.get("next_evidence_test"))
        print("  v17_readiness:", row.get("v17_readiness"))


if __name__ == "__main__":
    main()
