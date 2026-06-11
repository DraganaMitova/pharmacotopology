#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT" / "v19_kcsa_membrane_pore_evidence_readout_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V19 KcsA MEMBRANE/PORE EVIDENCE READOUT ===")
    for key in [
        "run_mode",
        "test_status",
        "target_id",
        "role_class",
        "positive_pressure_evidence_found",
        "positive_folding_evidence_found",
        "membrane_pore_role_evidence_found",
        "claim_allowed",
        "evidence_claim_allowed",
        "new_md_executed",
        "membrane_md_executed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "selected_core_or_clean_abstain",
        "pore_filter_annotation_present",
        "transmembrane_role_present",
        "oligomer_or_chain_interface_context_present",
        "membrane_environment_context_present",
        "soluble_core_misclassification_avoided",
        "whole_channel_fold_claim_made",
        "tetramer_claim_made",
    ]:
        print(f"{key}: {data.get(key)}")
    print("role_buckets_assigned:", data.get("role_buckets_assigned"))
    print("available_evidence:", data.get("available_evidence"))
    print("missing_evidence:", data.get("missing_evidence"))
    print("forbidden_misclassification_violations:", data.get("forbidden_misclassification_violations"))
    print("material_observations:", data.get("material_observations"))

    pore = data.get("pore_filter_readout", {})
    print("\nPore/filter readout:")
    print("  sequence_motif_probe:", pore.get("sequence_motif_probe"))
    print("  annotation_probe:", pore.get("annotation_probe"))
    print("  potassium_ion_probe:", pore.get("potassium_ion_probe"))
    print("  policy:", pore.get("policy"))

    helix = data.get("transmembrane_helix_readout", {})
    print("\nTransmembrane/helix readout:")
    print("  helix_record_count:", helix.get("helix_record_count"))
    print("  helix_records_by_chain:", helix.get("helix_records_by_chain"))
    print("  transmembrane_helix_scaffold_context_present:", helix.get("transmembrane_helix_scaffold_context_present"))
    print("  helix_policy:", helix.get("helix_policy"))

    iface = data.get("chain_interface_readout", {})
    print("\nChain/interface readout:")
    print("  chain_ca_counts:", iface.get("chain_ca_counts"))
    print("  interface_pairs_present:", iface.get("interface_pairs_present"))
    print("  min_ca_distance_by_chain_pair:", iface.get("min_ca_distance_by_chain_pair"))
    print("  contact_probe_policy:", iface.get("contact_probe_policy"))


if __name__ == "__main__":
    main()
