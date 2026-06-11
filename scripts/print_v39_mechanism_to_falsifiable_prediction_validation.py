#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V39_MECHANISM_TO_FALSIFIABLE_PREDICTION_VALIDATION" / "v39_mechanism_to_falsifiable_prediction_validation_certificate.json"


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V39 certificate: {CERT}")
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V39 MECHANISM TO FALSIFIABLE PREDICTION VALIDATION ===")
    print(f"V39 status: {cert.get('control_status')}")
    print(f"validated targets: {cert.get('validated_target_count')} / {cert.get('target_count')}")
    print(f"controls: {cert.get('passed_control_count')} / {cert.get('control_count')}")
    print(f"coordinate_derived_source_count: {cert.get('coordinate_derived_source_count')}")
    print(f"internal_runtime_source_count: {cert.get('internal_runtime_source_count')}")
    print(f"answer_key_used_for_prediction: {cert.get('answer_key_used_for_prediction')}")
    print("")
    print("Prediction / holdout summary:")
    for row in cert.get("target_validation_results", []):
        print(f"  - {row.get('target')}: {row.get('validation_level')}")
        print(f"    predictions: {row.get('prediction_count')}")
        print(f"    holdout sources: {row.get('holdout_source_count')}")
        print(f"    supported buckets: {row.get('supported_prediction_buckets')}")
    print("")
    print(f"certificate: {cert.get('artifacts', {}).get('certificate')}")
    print(f"report: {cert.get('artifacts', {}).get('report')}")
    print("")
    print("Plain English interpretation:")
    print(cert.get("locked_interpretation"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
