#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V43_PROTEIN_FOLDING_SOLVED_FLAG_TRIAL" / "v43_protein_folding_solved_flag_trial_certificate.json"

V43_PATHS = [
    "scripts/build_v43_solved_flag_trial_panel_v0.py",
    "scripts/run_v43_protein_folding_solved_flag_trial_v0.py",
    "scripts/print_v43_protein_folding_solved_flag_trial.py",
    "tests/test_v43_protein_folding_solved_flag_trial.py",
    "docs/V43_PROTEIN_FOLDING_SOLVED_FLAG_TRIAL_PROTOCOL.md",
    "data/solved_flag_trial/V43",
    "first_contact_clean_pharmacotopology_layer_run/V43_PROTEIN_FOLDING_SOLVED_FLAG_TRIAL",
]


def _git(args: list[str]) -> str:
    result = subprocess.run(["git", *args], cwd=REPO_ROOT, check=False, text=True, capture_output=True)
    return result.stdout.strip()


def _v42_commit_hash() -> str:
    return _git(["log", "--grep", "Add V42 de novo universal mechanism language challenge", "-n", "1", "--format=%H"])


def _v43_committed() -> bool:
    return _git(["status", "--short", "--", *V43_PATHS]) == "" and CERT.exists()


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V43 certificate: {CERT}")
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V43 PROTEIN FOLDING SOLVED FLAG TRIAL ===")
    print(f"V42 committed: {bool(_v42_commit_hash())}")
    print(f"V42 commit hash: {_v42_commit_hash()}")
    print(f"V43 status: {cert.get('control_status')}")
    print(f"protein_folding_solved_candidate: {cert.get('protein_folding_solved_candidate')}")
    print(f"folding_problem_solved: {cert.get('folding_problem_solved')}")
    print(f"claim_allowed: {cert.get('claim_allowed')}")
    print(f"panel target count: {cert.get('panel_target_count')}")
    print(f"holdout validated target count: {cert.get('holdout_validated_target_count')}")
    print(f"mechanism accuracy: {cert.get('mechanism_class_accuracy'):.3f}")
    print(f"hard precision / recall: {cert.get('hard_class_precision'):.3f} / {cert.get('hard_class_recall'):.3f}")
    print(f"operator support: {cert.get('operator_region_support_rate'):.3f}")
    print(f"contact precision / enrichment: {cert.get('top_L_long_range_contact_precision'):.3f} / {cert.get('contact_enrichment_vs_random'):.3f}")
    print(f"low-resolution structure/ensemble support: {cert.get('low_resolution_structure_or_ensemble_support_rate'):.3f}")
    print(f"perturbation support: {cert.get('perturbation_support_rate'):.3f}")
    print("")
    print("Baseline scores:")
    for name, score in cert.get("baseline_scores", {}).items():
        print(f"  - {name}: {score.get('mechanism_class_accuracy'):.3f}")
    print("")
    print("Threshold pass/fail table:")
    for row in cert.get("thresholds", []):
        print(f"  - {row.get('threshold')}: {row.get('passed')} observed={row.get('observed')} required={row.get('required')}")
    print("")
    print(f"failed targets: {cert.get('targets_failed')}")
    print(f"abstained targets: {cert.get('targets_abstained')}")
    print(f"controls: {cert.get('passed_control_count')} / {cert.get('control_count')}")
    print("tests passed: pending focused pytest command")
    print(f"leakage counts: coordinate={cert.get('coordinate_derived_source_count_before_prediction')} internal={cert.get('internal_runtime_source_count_for_prediction')} holdout={cert.get('holdout_leakage_detected')}")
    print("")
    print(f"certificate: {cert.get('artifacts', {}).get('certificate')}")
    print(f"report: {cert.get('artifacts', {}).get('report')}")
    print(f"V43 committed: {_v43_committed()}")
    print("")
    print("Plain English interpretation:")
    print(cert.get("locked_interpretation"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
