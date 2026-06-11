#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V38_BLIND_MECHANISM_GENERALIZATION_PANEL" / "v38_blind_mechanism_generalization_panel_certificate.json"


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V38 certificate: {CERT}")
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V38 BLIND MECHANISM GENERALIZATION PANEL ===")
    print(f"V38 status: {data.get('control_status')}")
    print(f"known positive accuracy: {data.get('correct_known_positive_count')} / {data.get('known_positive_count')}")
    print(f"false positive hard-grammar count: {data.get('false_positive_hard_grammar_count')}")
    print(f"false negative known-positive count: {data.get('false_negative_known_positive_count')}")
    print(f"controls: {data.get('passed_control_count')} / {data.get('control_count')}")
    print(f"coordinate_derived_source_count: {data.get('coordinate_derived_source_count')}")
    print(f"internal_runtime_source_count: {data.get('internal_runtime_source_count')}")
    print(f"placeholder_source_count: {data.get('placeholder_source_count')}")
    print(f"target names masked: {data.get('masked_panel_used')}")
    print(f"answer key used for assignment: {data.get('answer_key_used_for_assignment')}")
    print(f"claim_allowed: {data.get('claim_allowed')}")
    print(f"folding_problem_solved: {data.get('folding_problem_solved')}")
    print("")
    print("Masked target table:")
    for row in data.get("masked_assignment_results", []):
        print(f"  - {row.get('masked_id')}: {row.get('assigned_mechanism_class')} strength={row.get('confidence_strength')} abstain={row.get('abstain_reason')}")
    print("")
    print("Unmasked answer-key evaluation:")
    for row in data.get("answer_key_evaluation", []):
        print(f"  - {row.get('masked_id')} {row.get('target')}: expected={row.get('expected_mechanism_class')} assigned={row.get('assigned_mechanism_class')}")
    print("")
    print(f"certificate: {data.get('artifacts', {}).get('certificate')}")
    print(f"report: {data.get('artifacts', {}).get('report')}")
    print("")
    print("Plain English interpretation:")
    print(data.get("locked_interpretation"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
