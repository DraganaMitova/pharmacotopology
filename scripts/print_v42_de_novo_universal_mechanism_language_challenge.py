#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V42_DE_NOVO_UNIVERSAL_MECHANISM_LANGUAGE_CHALLENGE" / "v42_de_novo_universal_mechanism_language_challenge_certificate.json"
PANEL = REPO_ROOT / "data" / "de_novo_mechanism_language" / "V42" / "panel" / "panel_manifest.json"

V42_PATHS = [
    "scripts/build_v42_de_novo_universal_mechanism_challenge_panel_v0.py",
    "scripts/run_v42_de_novo_universal_mechanism_language_challenge_v0.py",
    "scripts/print_v42_de_novo_universal_mechanism_language_challenge.py",
    "tests/test_v42_de_novo_universal_mechanism_language_challenge.py",
    "docs/V42_DE_NOVO_UNIVERSAL_MECHANISM_LANGUAGE_CHALLENGE_PROTOCOL.md",
    "data/de_novo_mechanism_language/V42",
    "first_contact_clean_pharmacotopology_layer_run/V42_DE_NOVO_UNIVERSAL_MECHANISM_LANGUAGE_CHALLENGE",
]


def _git(args: list[str]) -> str:
    result = subprocess.run(["git", *args], cwd=REPO_ROOT, check=False, text=True, capture_output=True)
    return result.stdout.strip()


def _v42_committed() -> bool:
    return _git(["status", "--short", "--", *V42_PATHS]) == "" and CERT.exists()


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V42 certificate: {CERT}")
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    panel = json.loads(PANEL.read_text(encoding="utf-8")) if PANEL.exists() else {"panel_groups": {}}
    print("=== V42 DE NOVO UNIVERSAL MECHANISM LANGUAGE CHALLENGE ===")
    print(f"V42 status: {cert.get('control_status')}")
    print(f"panel target count: {cert.get('panel_target_count')}")
    print(f"sealed prediction count: {cert.get('sealed_prediction_count')}")
    print(f"mechanism-class accuracy: {cert.get('mechanism_class_accuracy'):.3f}")
    print(f"hard-class precision / recall: {cert.get('hard_class_precision'):.3f} / {cert.get('hard_class_recall'):.3f}")
    print(f"operator-region support rate: {cert.get('operator_region_support_rate'):.3f}")
    print(f"perturbation support rate: {cert.get('perturbation_prediction_support_rate'):.3f}")
    print(f"controls: {cert.get('passed_control_count')} / {cert.get('control_count')}")
    print("tests passed: pending focused pytest command")
    print("")
    print("Panel composition:")
    for group, count in sorted(panel.get("panel_groups", {}).items()):
        print(f"  - {group}: {count}")
    print("")
    print("Sealed prediction / holdout summary:")
    for row in cert.get("target_scoring_results", []):
        print(f"  - {row.get('target_id')} {row.get('target_name')}: {row.get('predicted_mechanism_class')} -> {row.get('validation_level')}")
    print("")
    print("Baseline comparison:")
    for name, score in cert.get("baseline_scores", {}).items():
        print(f"  - {name}: {score.get('mechanism_class_accuracy'):.3f}")
    print(f"beats baselines: random={cert.get('beats_random_baseline')} keyword={cert.get('beats_keyword_baseline')} majority={cert.get('beats_majority_baseline')}")
    print("")
    print(f"failed targets: {cert.get('targets_failed')}")
    print(f"abstained targets: {cert.get('targets_abstained')}")
    print(f"leakage counts: coordinate={cert.get('coordinate_derived_source_count_for_prediction')} internal={cert.get('internal_runtime_source_count_for_prediction')} holdout={cert.get('holdout_leakage_detected')}")
    print("")
    print(f"certificate: {cert.get('artifacts', {}).get('certificate')}")
    print(f"report: {cert.get('artifacts', {}).get('report')}")
    print(f"V42 committed: {_v42_committed()}")
    print("")
    print("Plain English interpretation:")
    print(cert.get("locked_interpretation"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
