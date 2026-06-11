#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V35_NONCOORDINATE_EVOLUTIONARY_HOLDOUT" / "v35_noncoordinate_evolutionary_holdout_certificate.json"


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V35 certificate: {CERT}")
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V35 NON-COORDINATE EVOLUTIONARY HOLDOUT ===")
    print(f"run_mode: {data.get('run_mode')}")
    print(f"control_status: {data.get('control_status')}")
    print(f"source_status: {data.get('source_status')}")
    print(f"selected_target: {data.get('selected_target')}")
    print(f"evidence_row_count: {data.get('evidence_row_count')}")
    print(f"noncoordinate_external_source_count: {data.get('noncoordinate_external_source_count')}")
    print(f"coordinate_derived_source_count: {data.get('coordinate_derived_source_count')}")
    print(f"internal_runtime_source_count: {data.get('internal_runtime_source_count')}")
    print(f"operator_candidate_found: {data.get('operator_candidate_found')}")
    print(f"controls: {data.get('passed_control_count')} / {data.get('control_count')}")
    print(f"claim_allowed: {data.get('claim_allowed')}")
    print(f"new_MD_allowed: {data.get('new_MD_allowed')}")
    print(f"new_md_executed: {data.get('new_md_executed')}")
    print(f"positive_folding_evidence_found: {data.get('positive_folding_evidence_found')}")
    print(f"folding_problem_solved: {data.get('folding_problem_solved')}")
    print(f"coordinate_truth_used_before_selection: {data.get('coordinate_truth_used_before_selection')}")
    print(f"native_metrics_used_for_selection: {data.get('native_metrics_used_for_selection')}")
    print("")
    print("Selected evidence sources:")
    sources = data.get("selected_source_paths") or []
    if sources:
        for source in sources:
            print(f"  - {source}")
    else:
        print("  - None")
    print("")
    print("Failed checks:")
    failures = data.get("failed_checks") or []
    if failures:
        for failure in failures:
            print(f"  - {failure}")
    else:
        print("  - None")
    print("")
    print("Controls:")
    for row in data.get("controls", []):
        print(f"  - {row.get('control_id')}: {'PASS' if row.get('passed') else 'FAIL'}")
        print(f"    observed_status: {row.get('observed_status')}")
        print(f"    expected_status: {row.get('expected_status')}")
    print("")
    print("Next decision:")
    print(json.dumps(data.get("next_decision", {}), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
