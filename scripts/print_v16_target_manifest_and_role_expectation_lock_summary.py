#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCK"
    / "v16_target_manifest_and_role_expectation_lock_certificate.json"
)


def main() -> None:
    if not CERT.exists():
        raise SystemExit(f"missing V16 lock certificate: {CERT}\nrun scripts/run_v16_target_manifest_and_role_expectation_lock.sh first")
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V16 TARGET MANIFEST AND ROLE-EXPECTATION LOCK ===")
    for key in [
        "run_mode",
        "lock_status",
        "source_v15_lock_status",
        "transfer_panel_claim",
        "claim_allowed",
        "data_preflight_executed",
        "new_md_executed",
        "fixed_threshold_policy",
        "threshold_policy",
        "native_metrics_not_used_for_selection",
        "fixed_residue_cutoff_used",
    ]:
        print(f"{key}: {cert.get(key)}")
    print("lock_failed_checks:", cert.get("lock_failed_checks"))
    print("\nLocked interpretation:")
    print(cert.get("locked_interpretation"))
    print("\nTargets:")
    for target in cert.get("locked_targets", []):
        print(f"\n{target.get('target_id')}:")
        print("  pressure_class:", target.get("pressure_class"))
        print("  expected_role_class:", target.get("expected_role_class"))
        print("  clean_abstain_allowed:", target.get("clean_abstain_allowed"))
        print("  claim_allowed:", target.get("claim_allowed"))
        print("  allowed_evidence_roles:")
        for role in target.get("allowed_evidence_roles", []):
            print("    -", role)
        print("  forbidden_misclassification:")
        for item in target.get("forbidden_misclassification", []):
            print("    -", item)


if __name__ == "__main__":
    main()
