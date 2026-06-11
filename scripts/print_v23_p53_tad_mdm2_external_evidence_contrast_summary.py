#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST" / "v23_p53_tad_mdm2_external_evidence_contrast_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V23 P53 TAD/MDM2 EXTERNAL EVIDENCE CONTRAST TEST ===")
    for key in [
        "run_mode", "test_status", "target_id", "role_class", "pressure_class",
        "positive_pressure_evidence_found", "positive_external_annotation_evidence_found",
        "positive_folding_evidence_found", "claim_allowed", "evidence_claim_allowed",
        "new_md_executed", "fixed_residue_cutoff_used", "native_metrics_used_for_selection",
        "selected_core_or_clean_abstain", "contrast_status",
        "external_disorder_reference_context_found", "partner_bound_interface_reference_context_found",
        "isolated_TAD_clean_abstain_preserved", "isolated_TAD_autonomous_core_selected",
        "isolated_TAD_fold_claim_made", "partner_bound_interface_evidence_preserved",
        "external_couplings_or_msa_files_present", "external_couplings_required_for_this_test",
        "forbidden_misclassification_violations", "failed_checks",
    ]:
        print(f"{key}: {data.get(key)}")
    print("available_evidence:", data.get("available_evidence"))
    print("missing_evidence:", data.get("missing_evidence"))
    print("\nState labels:")
    for k, v in data.get("state_labels", {}).items():
        print(f"  {k}: {v}")
    print("\nInterface readout:")
    interface = data.get("interface_contact_readout", {})
    for key in ["chain_pair", "chain_ca_counts", "min_ca_distance", "multi_radius_contact_counts", "contact_probe_policy", "interface_contact_evidence_present"]:
        print(f"  {key}: {interface.get(key)}")
    print("\nPartner-bound helix proxy:")
    helix = data.get("partner_bound_helix_proxy_readout", {})
    for key in ["small_chain", "small_chain_ca_count", "ca_i_to_i_plus_3_distance_summary", "ca_i_to_i_plus_4_distance_summary", "helix_proxy_policy"]:
        print(f"  {key}: {helix.get(key)}")
    print("\nExternal annotation source references:")
    for ref in data.get("external_annotation_source_references", []):
        print(f"  {ref.get('source_id')}: {ref.get('kind')} [{ref.get('use_boundary')}]")
    print("locked_interpretation:", data.get("locked_interpretation"))


if __name__ == "__main__":
    main()
