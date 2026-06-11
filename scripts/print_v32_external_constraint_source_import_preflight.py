#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

ROOT = Path(os.environ.get("REPO_ROOT", Path.cwd())).resolve()
CERT = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT" / "v32_external_constraint_source_import_preflight_certificate.json"


def _short(value):
    if isinstance(value, list) and len(value) > 8:
        return value[:8] + [f"...({len(value)} total)"]
    return value


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V32 EXTERNAL CONSTRAINT SOURCE IMPORT PREFLIGHT ===")
    for key in [
        "run_mode", "preflight_status", "source_v31_status", "import_manifest_present", "import_manifest_path",
        "claim_allowed", "new_md_executed", "membrane_md_executed", "new_MD_allowed", "new_MD_recommended",
        "fixed_threshold_policy", "target_specific_threshold_tuning_allowed", "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection", "selected_V31_targets", "selected_V33_target", "selected_V33_panel",
        "md_ready_targets", "provenance_clean", "preflight_failed_checks",
    ]:
        print(f"{key}: {data.get(key)}")

    print("\nImport requirements:")
    for target, reqs in data.get("import_requirements", {}).items():
        print(f"\n{target}:")
        for req in reqs:
            print(f"  - {req}")

    print("\nTarget import rows:")
    for row in data.get("target_rows", []):
        print(f"\n{row.get('target')}:")
        for key in [
            "import_row_count", "target_import_status", "target_ready_for_V33", "missing_for_V33",
            "valid_real_external_constraint_rows", "role_context_rows", "constraint_derivation_only_rows",
            "invalid_or_excluded_rows", "ready_for_MD", "new_MD_allowed", "claim_allowed",
        ]:
            print(f"  {key}: {_short(row.get(key))}")

    decision = data.get("next_constraint_backed_operator_readout_decision", {})
    print("\nNext decision:")
    for key in ["decision_status", "selected_V33_target", "selected_V33_panel", "next_action", "reason", "new_MD_allowed", "new_MD_recommended", "claim_allowed"]:
        print(f"  {key}: {decision.get(key)}")


if __name__ == "__main__":
    main()
