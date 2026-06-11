#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT" / "v30_external_constraint_and_coupling_acquisition_sprint_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V30 EXTERNAL CONSTRAINT AND COUPLING ACQUISITION SPRINT ===")
    for key in [
        "run_mode", "sprint_status", "source_v29_status", "claim_allowed", "new_md_executed", "membrane_md_executed",
        "new_MD_allowed", "new_MD_recommended", "fixed_threshold_policy", "target_specific_threshold_tuning_allowed",
        "fixed_residue_cutoff_used", "native_metrics_used_for_selection", "positive_pressure_evidence_targets",
        "positive_folding_evidence_targets", "constraint_or_coupling_ready_targets", "annotation_context_locked_targets",
        "md_ready_targets", "selected_next_panel", "selected_V31_targets", "sprint_failed_checks",
    ]:
        print(f"{key}: {data.get(key)}")

    print("\nTarget acquisition rows:")
    for row in data.get("target_rows", []):
        print(f"\n{row.get('target')}:")
        print(f"  pressure_type: {row.get('pressure_type')}")
        print(f"  selected_signal: {row.get('selected_signal')}")
        print(f"  acquisition_status: {row.get('acquisition_status')}")
        print(f"  local_constraints_or_couplings_present: {row.get('local_constraints_or_couplings_present')}")
        print(f"  local_annotation_context_present_or_locked: {row.get('local_annotation_context_present_or_locked')}")
        print(f"  candidate_coupling_or_constraint_files: {row.get('candidate_coupling_or_constraint_files')}")
        print(f"  candidate_annotation_files: {row.get('candidate_annotation_files')}")
        print(f"  missing_external_constraints_or_couplings: {row.get('missing_external_constraints_or_couplings')}")
        print(f"  recommended_next_action: {row.get('recommended_next_action')}")
        print(f"  ready_for_MD: {row.get('ready_for_MD')}")
        print(f"  new_MD_allowed: {row.get('new_MD_allowed')}")
        print(f"  claim_allowed: {row.get('claim_allowed')}")

    print("\nNext acquisition decision:")
    dec = data.get("next_acquisition_decision", {})
    for key in [
        "decision_status", "selected_next_panel", "selected_V31_targets", "primary_focus", "reason",
        "coupling_ready_targets", "md_ready_targets", "new_MD_allowed", "new_MD_recommended", "membrane_MD_recommended",
        "parallel_MD_paths_allowed", "required_before_any_MD", "claim_allowed",
    ]:
        print(f"  {key}: {dec.get(key)}")


if __name__ == "__main__":
    main()
