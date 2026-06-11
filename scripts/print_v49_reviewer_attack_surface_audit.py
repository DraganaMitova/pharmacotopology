#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V49_REVIEWER_ATTACK_SURFACE_HARDENING" / "v49_reviewer_attack_surface_audit_certificate.json"


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V49 certificate: {CERT}")
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V49 REVIEWER ATTACK-SURFACE HARDENING ===")
    print(f"V49 status: {cert.get('status')}")
    print(f"evidence_taxonomy_present: {cert.get('evidence_taxonomy_present')}")
    print(f"spatial_proxy_boundary_present: {cert.get('spatial_proxy_boundary_present')}")
    print(f"operator_scoring_rubric_present: {cert.get('operator_scoring_rubric_present')}")
    print(f"latest_completed_cycle_documented: {cert.get('latest_completed_cycle_documented')}")
    print(f"negative_control_inventory_present: {cert.get('negative_control_inventory_present')}")
    print(f"v46_full_cycle_visible: {cert.get('v46_full_cycle_visible')}")
    print(f"readme_claim_boundary_preserved: {cert.get('readme_claim_boundary_preserved')}")
    print(f"controls passed: {cert.get('passed_control_count')} / {cert.get('control_count')}")
    print(f"negative/null controls visible: {cert.get('negative_null_control_count')}")
    print(f"folding_problem_solved: {cert.get('folding_problem_solved')}")
    print(f"claim_allowed: {cert.get('claim_allowed')}")
    print(f"failed_checks: {cert.get('failed_checks')}")
    print(f"certificate: {cert.get('artifacts', {}).get('certificate')}")
    print(f"report: {cert.get('artifacts', {}).get('report')}")
    print(f"V46 certificate: {cert.get('artifacts', {}).get('v46_certificate')}")
    print(f"V46 report: {cert.get('artifacts', {}).get('v46_report')}")
    print("")
    print("Plain English interpretation:")
    print(cert.get("plain_english_interpretation"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
