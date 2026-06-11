#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V37_MECHANISM_QUESTION_PROBES" / "v37_mechanism_question_probes_certificate.json"


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V37 certificate: {CERT}")
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V37 MECHANISM QUESTION PROBES ===")
    print(f"V37 status: {data.get('control_status')}")
    print(f"input_v36_status: {data.get('input_v36_status')}")
    print(f"assigned_target_count: {data.get('assigned_target_count')} / {data.get('target_count')}")
    print(f"partial_target_count: {data.get('partial_target_count')}")
    print(f"blocked_target_count: {data.get('blocked_target_count')}")
    print(f"target_name_masking_passed: {data.get('target_name_masking_passed')}")
    print(f"swapped_dossier_detection_passed: {data.get('swapped_dossier_detection_passed')}")
    print(f"coordinate_derived_source_count: {data.get('coordinate_derived_source_count')}")
    print(f"internal_runtime_source_count: {data.get('internal_runtime_source_count')}")
    print(f"placeholder_source_count: {data.get('placeholder_source_count')}")
    print(f"controls: {data.get('passed_control_count')} / {data.get('control_count')}")
    print(f"claim_allowed: {data.get('claim_allowed')}")
    print(f"new_MD_allowed: {data.get('new_MD_allowed')}")
    print(f"new_md_executed: {data.get('new_md_executed')}")
    print(f"positive_folding_evidence_found: {data.get('positive_folding_evidence_found')}")
    print(f"folding_problem_solved: {data.get('folding_problem_solved')}")
    print("")
    print("Mechanism classes:")
    for target, mechanism in data.get("mechanism_classes_assigned", {}).items():
        print(f"  - {target}: {mechanism}")
    print("")
    print("Target details:")
    for row in data.get("target_results", []):
        print(f"  - {row.get('target')}: {row.get('target_status')}")
        print(f"    required operators: {row.get('required_operator_buckets_present')}")
        print(f"    forbidden operators: {row.get('forbidden_operator_buckets_triggered')}")
        print("    scientific questions:")
        for answer in row.get("scientific_question_answers", []):
            print(f"      - {answer.get('answer')}: {answer.get('question')}")
    print("")
    print(f"certificate: {data.get('artifacts', {}).get('certificate')}")
    print(f"report: {data.get('artifacts', {}).get('report')}")
    print(f"answer_table: {data.get('artifacts', {}).get('answer_table')}")
    print("")
    print("Plain English interpretation:")
    print(data.get("locked_interpretation"))
    print("")
    print("Next decision:")
    print(json.dumps(data.get("next_decision", {}), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
