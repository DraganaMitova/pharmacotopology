#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V41_PROTEIN_FOLDING_CLAIM_GATE" / "v41_protein_folding_claim_gate_certificate.json"

V41_PATHS = [
    "scripts/run_v41_protein_folding_claim_gate_v0.py",
    "scripts/print_v41_protein_folding_claim_gate.py",
    "tests/test_v41_protein_folding_claim_gate.py",
    "docs/V41_PROTEIN_FOLDING_CLAIM_GATE_PROTOCOL.md",
    "docs/CLAIM_BOUNDARY_CURRENT.md",
    "data/claim_gate/V41",
    "first_contact_clean_pharmacotopology_layer_run/V41_PROTEIN_FOLDING_CLAIM_GATE",
]


def _git(args: list[str]) -> str:
    result = subprocess.run(["git", *args], cwd=REPO_ROOT, check=False, text=True, capture_output=True)
    return result.stdout.strip()


def _v40_commit_hash() -> str:
    return _git(["log", "--grep", "Add V40 mechanism perturbation pressure tests", "-n", "1", "--format=%H"])


def _v41_committed() -> bool:
    return _git(["status", "--short", "--", *V41_PATHS]) == "" and CERT.exists()


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V41 certificate: {CERT}")
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V41 PROTEIN FOLDING CLAIM GATE ===")
    print(f"V40 committed: {bool(_v40_commit_hash())}")
    print(f"V40 commit hash: {_v40_commit_hash()}")
    print(f"V41 status: {cert.get('control_status')}")
    print(f"max allowed claim level: {cert.get('max_allowed_claim_level')}")
    print(f"protein_folding_solved_claim_allowed: {cert.get('protein_folding_solved_claim_allowed')}")
    print(f"mechanism_claim_allowed: {cert.get('mechanism_claim_allowed')}")
    print(f"controls: {cert.get('passed_control_count')} / {cert.get('control_count')}")
    print(f"tests passed: pending focused pytest command")
    print("")
    print("Allowed claim text:")
    print(cert.get("allowed_public_claim"))
    print("")
    print("Forbidden claims:")
    for claim in cert.get("forbidden_public_claims", []):
        print(f"  - {claim}")
    print("")
    print("Missing evidence for C5:")
    for item in cert.get("evidence_missing_for_c5", []):
        print(f"  - {item}")
    print("")
    print("Scientist attack surface:")
    for item in cert.get("scientist_attack_surface", []):
        print(f"  - {item}")
    print("")
    print(f"certificate: {cert.get('artifacts', {}).get('certificate')}")
    print(f"report: {cert.get('artifacts', {}).get('report')}")
    print(f"V41 committed: {_v41_committed()}")
    print("")
    print("Plain English interpretation:")
    print(cert.get("locked_interpretation"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
