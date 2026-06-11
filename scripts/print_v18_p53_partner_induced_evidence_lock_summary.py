#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCK" / "v18_p53_partner_induced_evidence_lock_certificate.json"

def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V18 P53 PARTNER-INDUCED EVIDENCE LOCK ===")
    for key in [
        "run_mode",
        "lock_status",
        "source_v18_status",
        "target_id",
        "role_class",
        "locked_claim",
        "positive_pressure_evidence_found",
        "partner_induced_role_evidence_found",
        "positive_folding_evidence_found",
        "autonomous_TAD_fold_claim",
        "claim_allowed",
        "new_md_executed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "isolated_TAD_autonomous_fold_status",
        "selected_core_or_clean_abstain",
    ]:
        print(f"{key}: {data.get(key)}")
    print("forbidden_misclassification_violations:", data.get("forbidden_misclassification_violations"))
    print("lock_failed_checks:", data.get("lock_failed_checks"))
    print("available_evidence:", data.get("available_evidence"))
    print("missing_evidence:", data.get("missing_evidence"))
    iface = data.get("interface_contact_readout", {})
    print("\nInterface lock summary:")
    print("  chain_pair:", iface.get("chain_pair"))
    print("  min_ca_distance:", iface.get("min_ca_distance"))
    print("  multi_radius_contact_counts:", iface.get("multi_radius_contact_counts"))
    print("  contact_probe_policy:", iface.get("contact_probe_policy"))
    print("next_step:", data.get("next_step"))

if __name__ == "__main__":
    main()
