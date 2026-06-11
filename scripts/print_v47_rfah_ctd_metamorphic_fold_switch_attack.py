#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_ATTACK" / "v47_rfah_ctd_metamorphic_fold_switch_attack_certificate.json"

V47_PATHS = [
    "scripts/build_v47_rfah_ctd_metamorphic_sources_v0.py",
    "scripts/run_v47_rfah_ctd_metamorphic_fold_switch_attack_v0.py",
    "scripts/print_v47_rfah_ctd_metamorphic_fold_switch_attack.py",
    "tests/test_v47_rfah_ctd_metamorphic_fold_switch_attack.py",
    "docs/V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_PROTOCOL.md",
    "data/live_unsolved_targets/V47/RfaH_CTD",
    "first_contact_clean_pharmacotopology_layer_run/V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_ATTACK",
]


def _git(args: list[str]) -> str:
    result = subprocess.run(["git", *args], cwd=REPO_ROOT, check=False, text=True, capture_output=True)
    return result.stdout.strip()


def _v47_committed() -> bool:
    return _git(["status", "--short", "--", *V47_PATHS]) == "" and CERT.exists()


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V47 certificate: {CERT}")
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V47 RFAH-CTD METAMORPHIC FOLD-SWITCH ATTACK ===")
    print(f"V46 commit hash: {cert.get('v46_commit_hash')}")
    print(f"V47 status: {cert.get('control_status')}")
    print(f"live_fold_switch_solution_packet: {cert.get('live_fold_switch_solution_packet')}")
    print(f"protein_folding_solved_candidate_strengthened: {cert.get('protein_folding_solved_candidate_strengthened')}")
    print(f"target and region: {cert.get('target')} / {cert.get('sequence_region_scope')}")
    print(f"mechanism class: {cert.get('mechanism_class')}")
    print(f"operator bucket count: {cert.get('operator_bucket_count')}")
    print("operator buckets:")
    for bucket in cert.get("operator_buckets", []):
        print(f"  - {bucket}")
    print(f"perturbation prediction count: {cert.get('perturbation_prediction_count')}")
    print("perturbation predictions:")
    for row in cert.get("perturbation_predictions", []):
        print(f"  - {row.get('prediction_id')}: {row.get('perturbation')} -> {row.get('predicted_effect')}")
    print("")
    print("validation support per prediction:")
    for row in cert.get("validation_support_per_prediction", []):
        print(f"  - {row.get('prediction_id')}: {row.get('support_level')}")
    print(f"validated prediction count: {cert.get('validated_prediction_count')}")
    print(f"contradicted prediction count: {cert.get('contradicted_prediction_count')}")
    print(f"contradicted predictions: {cert.get('contradicted_predictions')}")
    print(f"operator support rate: {cert.get('operator_support_rate'):.3f}")
    print(f"perturbation support rate: {cert.get('perturbation_support_rate'):.3f}")
    print(f"state separation support rate: {cert.get('state_separation_support_rate'):.3f}")
    print(f"partner context support rate: {cert.get('partner_context_support_rate'):.3f}")
    print("")
    print("proposed experiments:")
    for row in cert.get("proposed_experimental_tests", []):
        print(f"  - {row}")
    print("")
    print("falsification criteria:")
    for row in cert.get("falsification_criteria", []):
        print(f"  - {row}")
    print("")
    print(f"leakage counts: coordinate={cert.get('coordinate_derived_source_count_before_prediction')} internal={cert.get('internal_runtime_source_count_for_prediction')} holdout={cert.get('holdout_leakage_detected')} native_metrics={cert.get('native_metrics_used_before_prediction')} coordinate_truth={cert.get('coordinate_truth_used_before_prediction')}")
    print(f"controls passed: {cert.get('passed_control_count')} / {cert.get('control_count')}")
    print("tests passed: pending focused pytest command")
    print(f"certificate: {cert.get('artifacts', {}).get('certificate')}")
    print(f"report: {cert.get('artifacts', {}).get('report')}")
    print(f"V46 committed: {cert.get('v46_committed')}")
    print(f"V47 committed: {_v47_committed()}")
    print("")
    print("Plain English interpretation:")
    print(cert.get("plain_english_interpretation"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
