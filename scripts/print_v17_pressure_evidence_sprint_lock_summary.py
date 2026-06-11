#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERT = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "V17_PRESSURE_EVIDENCE_SPRINT_LOCK"
    / "v17_pressure_evidence_sprint_lock_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Print V17 pressure evidence sprint lock summary.")
    parser.add_argument("--cert", default=str(DEFAULT_CERT))
    args = parser.parse_args()
    path = Path(args.cert)
    if not path.exists():
        raise SystemExit(f"missing V17 sprint certificate: {path}")
    o = json.loads(path.read_text(encoding="utf-8"))
    print("=== V17 PRESSURE EVIDENCE SPRINT LOCK ===")
    for key in [
        "run_mode",
        "sprint_lock_status",
        "source_gap_lock_status",
        "source_zero_md_status",
        "source_preflight_status",
        "claim_allowed",
        "new_md_executed",
        "fixed_threshold_policy",
        "threshold_policy",
        "target_specific_threshold_tuning_allowed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "positive_folding_evidence_targets",
        "ready_for_V18_targets",
        "blocked_or_deferred_targets",
        "selected_V18_target",
        "selected_V18_test",
    ]:
        print(f"{key}: {o.get(key)}")
    print("sprint_failed_checks:", o.get("sprint_failed_checks"))
    print()
    print("Next target decision:")
    d = o.get("next_target_decision", {})
    if isinstance(d, dict):
        for key in ["decision_status", "selected_V18_target", "selected_V18_test", "reason", "v18_first_run_mode", "parallel_md_paths_allowed"]:
            print(f"  {key}: {d.get(key)}")
    print()
    print("Target sprint rows:")
    for row in o.get("target_rows", []):
        print(f"\n{row.get('target_id')}:")
        print("  role_class:", row.get("role_class"))
        print("  role_classification_passed:", row.get("role_classification_passed"))
        print("  zero_MD_readout_status:", row.get("zero_MD_readout_status"))
        print("  ready_for_V18:", row.get("ready_for_V18"))
        print("  next_evidence_test:", row.get("next_evidence_test"))
        print("  available_evidence:")
        for item in row.get("available_evidence", []):
            print("    -", item)
        print("  missing_for_first_v18_test:")
        for item in row.get("missing_for_first_v18_test", []):
            print("    -", item)
        print("  missing_evidence:")
        for item in row.get("missing_evidence", []):
            print("    -", item)
        obs = row.get("material_observations")
        if obs:
            print("  material_observations:", obs)


if __name__ == "__main__":
    main()
