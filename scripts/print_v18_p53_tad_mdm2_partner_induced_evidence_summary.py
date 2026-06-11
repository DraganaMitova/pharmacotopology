#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST" / "v18_p53_tad_mdm2_partner_induced_evidence_certificate.json"

def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V18 P53 TAD/MDM2 PARTNER-INDUCED EVIDENCE TEST ===")
    for key in [
        "run_mode",
        "test_status",
        "target_id",
        "role_class",
        "positive_pressure_evidence_found",
        "positive_folding_evidence_found",
        "partner_induced_role_evidence_found",
        "claim_allowed",
        "evidence_claim_allowed",
        "new_md_executed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "selected_core_or_clean_abstain",
        "isolated_TAD_autonomous_fold_status",
    ]:
        print(f"{key}: {data.get(key)}")
    print("forbidden_misclassification_violations:", data.get("forbidden_misclassification_violations"))
    print("available_evidence:", data.get("available_evidence"))
    print("missing_evidence:", data.get("missing_evidence"))
    iface = data.get("interface_contact_readout", {})
    print("\nInterface readout:")
    print("  chain_pair:", iface.get("chain_pair"))
    print("  chain_ca_counts:", iface.get("chain_ca_counts"))
    print("  min_ca_distance:", iface.get("min_ca_distance"))
    print("  multi_radius_contact_counts:", iface.get("multi_radius_contact_counts"))
    print("  contact_probe_policy:", iface.get("contact_probe_policy"))
    print("  interface_contact_evidence_present:", iface.get("interface_contact_evidence_present"))
    helix = data.get("partner_bound_helix_proxy_readout", {})
    print("\nPartner-bound helix proxy:")
    print("  small_chain:", helix.get("small_chain"))
    print("  small_chain_ca_count:", helix.get("small_chain_ca_count"))
    print("  ca_i_to_i_plus_3_distance_summary:", helix.get("ca_i_to_i_plus_3_distance_summary"))
    print("  ca_i_to_i_plus_4_distance_summary:", helix.get("ca_i_to_i_plus_4_distance_summary"))
    print("  helix_proxy_policy:", helix.get("helix_proxy_policy"))

if __name__ == "__main__":
    main()
