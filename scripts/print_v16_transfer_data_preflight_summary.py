#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERT = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "V16_TRANSFER_DATA_PREFLIGHT"
    / "v16_transfer_data_preflight_certificate.json"
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Print V16 transfer data preflight summary.")
    parser.add_argument("--cert", default=str(DEFAULT_CERT))
    args = parser.parse_args()
    path = Path(args.cert)
    if not path.exists():
        raise SystemExit(f"missing V16 transfer data preflight certificate: {path}")
    o = json.loads(path.read_text(encoding="utf-8"))

    print("=== V16 TRANSFER DATA PREFLIGHT ===")
    for key in [
        "run_mode",
        "data_preflight_status",
        "source_v16_lock_status",
        "claim_allowed",
        "data_preflight_executed",
        "new_md_executed",
        "download_attempted_by_this_script",
        "fixed_threshold_policy",
        "native_metrics_not_used_for_selection",
        "fixed_residue_cutoff_used",
    ]:
        print(f"{key}: {o.get(key)}")
    print("ready_targets:", o.get("ready_targets"))
    print("blocked_targets:", o.get("blocked_targets"))
    print("missing_material_by_target:", o.get("missing_material_by_target"))
    print()
    print("Targets:")
    for t in o.get("target_results", []):
        print(f"\n{t.get('target_id')}:")
        print("  pressure_class:", t.get("pressure_class"))
        print("  expected_role_class:", t.get("expected_role_class"))
        print("  target_material_status:", t.get("target_material_status"))
        print("  target_role_preflight_status:", t.get("target_role_preflight_status"))
        print("  clean_abstain_allowed:", t.get("clean_abstain_allowed"))
        print("  claim_allowed:", t.get("claim_allowed"))
        print("  missing_required_material_ids:", t.get("missing_required_material_ids"))
        print("  optional_or_later_material:", t.get("optional_or_later_material"))
        print("  forbidden_misclassification:", t.get("forbidden_misclassification"))
        for mat in t.get("material_results", []):
            pdb = mat.get("pdb_summary", {})
            print(f"  material {mat.get('material_id')} ({mat.get('pdb_id')}):")
            print("    status:", mat.get("material_status"))
            print("    failed_checks:", mat.get("failed_checks"))
            print("    ca_atom_count:", pdb.get("ca_atom_count"))
            print("    chain_ca_counts:", pdb.get("chain_ca_counts"))
            print("    provenance_safe:", mat.get("provenance_summary", {}).get("usage_boundary_safe"))


if __name__ == "__main__":
    main()
