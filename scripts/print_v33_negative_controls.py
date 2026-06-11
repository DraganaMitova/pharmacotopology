#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V33_NEGATIVE_CONTROLS" / "v33_negative_controls_certificate.json"


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V33 negative controls certificate: {CERT}")
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V33 NEGATIVE CONTROLS ===")
    print(f"run_mode: {data.get('run_mode')}")
    print(f"negative_control_status: {data.get('negative_control_status')}")
    print(f"passed_control_count: {data.get('passed_control_count')} / {data.get('control_count')}")
    print(f"claim_allowed: {data.get('claim_allowed')}")
    print(f"new_md_executed: {data.get('new_md_executed')}")
    print(f"membrane_md_executed: {data.get('membrane_md_executed')}")
    print(f"new_MD_allowed: {data.get('new_MD_allowed')}")
    print(f"positive_folding_evidence_found: {data.get('positive_folding_evidence_found')}")
    print(f"folding_problem_solved: {data.get('folding_problem_solved')}")
    print("")
    print("Controls:")
    for row in data.get("controls", []):
        print(f"  - {row.get('control_id')}: {'PASS' if row.get('passed') else 'FAIL'}")
        print(f"    observed: {row.get('observed_status')}")
        print(f"    expected: {row.get('expected_status')}")
        if row.get("observed_missing_operator_evidence"):
            print(f"    missing_operator_evidence: {row.get('observed_missing_operator_evidence')}")
        if row.get("observed_failed_checks"):
            print(f"    failed_checks: {row.get('observed_failed_checks')}")
        if row.get("observed_selected_V33_target") is not None:
            print(f"    selected_V33_target: {row.get('observed_selected_V33_target')}")
        print(f"    reason: {row.get('reason')}")
    print("")
    print("Next decision:")
    print(json.dumps(data.get("next_decision", {}), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
