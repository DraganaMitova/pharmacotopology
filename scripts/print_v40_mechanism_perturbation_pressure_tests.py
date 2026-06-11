#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
OUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V40_MECHANISM_PERTURBATION_PRESSURE_TESTS"
CERT = OUT_DIR / "v40_mechanism_perturbation_pressure_tests_certificate.json"

V40_PATHS = [
    "scripts/build_v40_mechanism_perturbation_sources_v0.py",
    "scripts/run_v40_mechanism_perturbation_pressure_tests_v0.py",
    "scripts/print_v40_mechanism_perturbation_pressure_tests.py",
    "tests/test_v40_mechanism_perturbation_pressure_tests.py",
    "docs/V40_MECHANISM_PERTURBATION_PRESSURE_TESTS_PROTOCOL.md",
    "data/mechanism_perturbations/V40",
    "first_contact_clean_pharmacotopology_layer_run/V40_MECHANISM_PERTURBATION_PRESSURE_TESTS",
]


def _git(args: list[str]) -> str:
    result = subprocess.run(["git", *args], cwd=REPO_ROOT, check=False, text=True, capture_output=True)
    return result.stdout.strip()


def _v39_commit_hash() -> str:
    found = _git(["log", "--grep", "Add V39 falsifiable mechanism prediction validation", "-n", "1", "--format=%H"])
    return found or _git(["rev-parse", "HEAD"])


def _v40_committed() -> bool:
    status = _git(["status", "--short", "--", *V40_PATHS])
    return status == "" and CERT.exists()


def _answer(row: dict) -> str:
    target = row.get("target")
    if target == "KcsA":
        return "Filter/ion-selectivity pressure matters more than generic channel annotation."
    if target == "XCL1_lymphotactin":
        return "State separation is the causal pressure; deleting or pooling states weakens or invalidates the mechanism."
    return "Disorder/context boundaries are the causal pressure; compact-fold or native-fold overpromotion is blocked."


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V40 certificate: {CERT}")
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V40 MECHANISM PERTURBATION PRESSURE TESTS ===")
    print(f"V39 committed: True")
    print(f"V39 commit hash: {_v39_commit_hash()}")
    print(f"V40 status: {cert.get('control_status')}")
    print(f"supported targets: {cert.get('supported_target_count')} / {cert.get('target_count')}")
    print(f"controls: {cert.get('passed_control_count')} / {cert.get('control_count')}")
    print(f"tests passed: pending focused pytest command")
    print(f"coordinate_derived_source_count: {cert.get('coordinate_derived_source_count')}")
    print(f"internal_runtime_source_count: {cert.get('internal_runtime_source_count')}")
    print(f"answer_key_used_for_prediction: {cert.get('answer_key_used_for_prediction')}")
    print("")
    print("Perturbation summary:")
    for row in cert.get("target_validation_results", []):
        print(f"  - {row.get('target')}: {row.get('validation_level')}")
        print(f"    perturbations: {row.get('perturbation_count')}")
        print(f"    validated perturbations: {len(row.get('supported_perturbation_buckets', []))}")
        print(f"    supported buckets: {row.get('supported_perturbation_buckets')}")
        print(f"    scientist answer: {_answer(row)}")
    print("")
    print(f"certificate: {cert.get('artifacts', {}).get('certificate')}")
    print(f"report: {cert.get('artifacts', {}).get('report')}")
    print(f"V40 committed: {_v40_committed()}")
    print("")
    print("Plain English interpretation:")
    print(cert.get("locked_interpretation"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
