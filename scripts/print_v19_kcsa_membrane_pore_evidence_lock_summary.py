#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCK" / "v19_kcsa_membrane_pore_evidence_lock_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V19 KcsA MEMBRANE/PORE EVIDENCE LOCK ===")
    for key in [
        "run_mode",
        "lock_status",
        "source_v19_status",
        "target_id",
        "role_class",
        "locked_claim",
        "positive_pressure_evidence_found",
        "membrane_pore_role_evidence_found",
        "positive_folding_evidence_found",
        "soluble_core_misclassification_avoided",
        "whole_channel_fold_claim_made",
        "tetramer_claim_made",
        "claim_allowed",
        "new_md_executed",
        "membrane_md_executed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
    ]:
        print(f"{key}: {data.get(key)}")
    print("role_buckets_assigned:", data.get("role_buckets_assigned"))
    print("available_evidence:", data.get("available_evidence"))
    print("missing_evidence:", data.get("missing_evidence"))
    print("forbidden_misclassification_violations:", data.get("forbidden_misclassification_violations"))
    print("lock_failed_checks:", data.get("lock_failed_checks"))
    print("locked_interpretation:", data.get("locked_interpretation"))
    print("next_step:", data.get("next_step"))


if __name__ == "__main__":
    main()
