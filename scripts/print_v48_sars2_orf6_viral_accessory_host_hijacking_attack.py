#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_ATTACK" / "v48_sars2_orf6_viral_accessory_host_hijacking_attack_certificate.json"

V48_PATHS = [
    "scripts/build_v48_sars2_orf6_viral_accessory_sources_v0.py",
    "scripts/run_v48_sars2_orf6_viral_accessory_host_hijacking_attack_v0.py",
    "scripts/print_v48_sars2_orf6_viral_accessory_host_hijacking_attack.py",
    "tests/test_v48_sars2_orf6_viral_accessory_host_hijacking_attack.py",
    "docs/V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_PROTOCOL.md",
    "data/live_unsolved_targets/V48/SARS2_ORF6",
    "first_contact_clean_pharmacotopology_layer_run/V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_ATTACK",
]


def _git(args: list[str]) -> str:
    result = subprocess.run(["git", *args], cwd=REPO_ROOT, check=False, text=True, capture_output=True)
    return result.stdout.strip()


def _v48_committed() -> bool:
    return _git(["status", "--short", "--", *V48_PATHS]) == "" and CERT.exists()


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V48 certificate: {CERT}")
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V48 SARS2 ORF6 VIRAL ACCESSORY HOST-HIJACKING ATTACK ===")
    print(f"V47 commit hash: {cert.get('v47_commit_hash')}")
    print(f"V48 status: {cert.get('control_status')}")
    print(f"live_viral_host_hijacking_solution_packet: {cert.get('live_viral_host_hijacking_solution_packet')}")
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
    print(f"host interaction support rate: {cert.get('host_interaction_support_rate'):.3f}")
    print(f"functional consequence support rate: {cert.get('functional_consequence_support_rate'):.3f}")
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
    print(f"V47 committed: {cert.get('v47_committed')}")
    print(f"V48 committed: {_v48_committed()}")
    print("")
    print("Plain English interpretation:")
    print(cert.get("plain_english_interpretation"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
