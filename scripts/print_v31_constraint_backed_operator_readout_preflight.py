#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

ROOT = Path(os.environ.get("REPO_ROOT", Path.cwd())).resolve()
CERT = ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V31_CONSTRAINT_BACKED_OPERATOR_READOUT_PREFLIGHT" / "v31_constraint_backed_operator_readout_preflight_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V31 CONSTRAINT-BACKED OPERATOR READOUT PREFLIGHT ===")
    for key in [
        "run_mode", "preflight_status", "source_v30_status", "claim_allowed", "new_md_executed",
        "membrane_md_executed", "new_MD_allowed", "new_MD_recommended", "fixed_threshold_policy",
        "target_specific_threshold_tuning_allowed", "fixed_residue_cutoff_used", "native_metrics_used_for_selection",
        "selected_V31_targets", "real_external_constraint_targets_in_selected_scope", "real_external_constraint_targets_all",
        "selected_V32_target", "selected_V32_panel", "md_ready_targets", "provenance_clean", "preflight_failed_checks",
    ]:
        print(f"{key}: {data.get(key)}")
    print("\nAllowed-use matrix:")
    for source_class, allowed in data.get("allowed_use_matrix", {}).items():
        print(f"  {source_class}: {allowed}")
    print("\nTarget preflight rows:")
    for row in data.get("target_rows", []):
        print(f"\n{row.get('target')}:")
        for key in [
            "in_selected_V31_scope", "v30_acquisition_status", "candidate_file_count", "preflight_status",
            "constraint_backed_operator_readout_allowed", "real_external_constraint_or_coupling_files",
            "constraint_derivation_only_files", "annotation_context_files", "generated_internal_report_files", "excluded_files",
            "allowed_use_counts", "evidence_source_class_counts", "ready_for_MD", "new_MD_allowed", "claim_allowed",
        ]:
            value = row.get(key)
            if isinstance(value, list) and len(value) > 8:
                value = value[:8] + [f"...({len(row.get(key, []))} total)"]
            print(f"  {key}: {value}")
    decision = data.get("next_operator_readout_decision", {})
    print("\nNext operator readout decision:")
    for key in ["decision_status", "selected_V32_target", "selected_V32_panel", "next_action", "reason", "new_MD_allowed", "new_MD_recommended", "claim_allowed"]:
        print(f"  {key}: {decision.get(key)}")


if __name__ == "__main__":
    main()
