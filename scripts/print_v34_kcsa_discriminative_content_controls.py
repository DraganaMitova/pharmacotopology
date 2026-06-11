#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS" / "v34_kcsa_discriminative_content_controls_certificate.json"


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V34 certificate: {CERT}")
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V34 KCSA DISCRIMINATIVE CONTENT CONTROLS ===")
    print(f"run_mode: {data.get('run_mode')}")
    print(f"control_status: {data.get('control_status')}")
    print(f"passed_control_count: {data.get('passed_control_count')} / {data.get('control_count')}")
    print(f"claim_allowed: {data.get('claim_allowed')}")
    print(f"new_md_executed: {data.get('new_md_executed')}")
    print(f"membrane_md_executed: {data.get('membrane_md_executed')}")
    print(f"new_MD_allowed: {data.get('new_MD_allowed')}")
    print(f"positive_folding_evidence_found: {data.get('positive_folding_evidence_found')}")
    print(f"folding_problem_solved: {data.get('folding_problem_solved')}")
    print("")
    print("Baseline content signature:")
    summary = data.get("baseline_content_signature", {})
    print(f"  pore_row_count: {summary.get('pore_row_count')}")
    print(f"  interface_row_count: {summary.get('interface_row_count')}")
    print(f"  pore_residue_numbers: {summary.get('pore_residue_numbers')}")
    print(f"  pore_chains: {summary.get('pore_chains')}")
    print(f"  pore_motif_labels: {summary.get('pore_motif_labels')}")
    print(f"  pore_ion_names: {summary.get('pore_ion_names')}")
    print(f"  interface_chain_pair_count: {summary.get('interface_chain_pair_count')}")
    print("")
    print("Controls:")
    for row in data.get("controls", []):
        print(f"  - {row.get('control_id')}: {'PASS' if row.get('passed') else 'FAIL'}")
        print(f"    observed_valid: {row.get('observed_valid')}")
        print(f"    expected_valid: {row.get('expected_valid')}")
        print(f"    observed_failures: {row.get('observed_failures')}")
        print(f"    reason: {row.get('reason')}")
    print("")
    print("Next decision:")
    print(json.dumps(data.get("next_decision", {}), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
