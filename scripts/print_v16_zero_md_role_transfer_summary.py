#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERT = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "V16_ZERO_MD_ROLE_TRANSFER_READOUT"
    / "v16_zero_md_role_transfer_readout_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Print V16 zero-MD role transfer summary.")
    parser.add_argument("--cert", default=str(DEFAULT_CERT))
    args = parser.parse_args()
    path = Path(args.cert)
    if not path.exists():
        raise SystemExit(f"missing V16 zero-MD role transfer certificate: {path}")
    o = json.loads(path.read_text(encoding="utf-8"))

    print("=== V16 ZERO-MD ROLE TRANSFER READOUT ===")
    for key in [
        "run_mode",
        "role_transfer_status",
        "source_preflight_status",
        "source_v16_lock_status",
        "claim_allowed",
        "new_md_executed",
        "fixed_threshold_policy",
        "target_specific_threshold_tuning_allowed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "selection_policy",
    ]:
        print(f"{key}: {o.get(key)}")
    print("positive_role_targets:", o.get("positive_role_targets"))
    print("clean_abstain_targets:", o.get("clean_abstain_targets"))
    print("blocked_targets:", o.get("blocked_targets"))
    print("forbidden_misclassification_violations:", o.get("forbidden_misclassification_violations"))
    print()
    print("Targets:")
    for row in o.get("target_rows", []):
        print(f"\n{row.get('target_id')}:")
        print("  pressure_class:", row.get("pressure_class"))
        print("  expected_role_class:", row.get("expected_role_class"))
        print("  target_transfer_status:", row.get("target_transfer_status"))
        print("  selected_core_or_clean_abstain:", row.get("selected_core_or_clean_abstain"))
        print("  positive_role_context_found:", row.get("positive_role_context_found"))
        print("  role_buckets_assigned:", row.get("role_buckets_assigned"))
        print("  monitor_only_roles:", row.get("monitor_only_roles"))
        print("  clean_abstain_roles:", row.get("clean_abstain_roles"))
        print("  forbidden_misclassification_violations:", row.get("forbidden_misclassification_violations"))
        print("  material_observations:", row.get("material_observations"))
        print("  limitations:", row.get("limitations"))


if __name__ == "__main__":
    main()
