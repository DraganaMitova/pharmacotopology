#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT" / "v20_xcl1_state_specific_evidence_readout_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V20 XCL1 STATE-SPECIFIC EVIDENCE READOUT ===")
    for key in [
        "run_mode",
        "test_status",
        "target_id",
        "role_class",
        "positive_pressure_evidence_found",
        "state_specific_role_evidence_found",
        "state_A_role_evidence_found",
        "state_B_role_evidence_found",
        "positive_folding_evidence_found",
        "claim_allowed",
        "evidence_claim_allowed",
        "new_md_executed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "selected_core_or_clean_abstain",
        "state_labels_present",
        "monomer_dimer_context_present",
        "mixed_state_fake_core_forbidden",
        "single_fold_forcing_forbidden",
        "mixed_state_contact_pooling_used",
        "mixed_state_pollution",
        "single_fold_claim_made",
        "fold_switch_claim_made",
    ]:
        print(f"{key}: {data.get(key)}")
    print("role_buckets_assigned:", data.get("role_buckets_assigned"))
    print("available_evidence:", data.get("available_evidence"))
    print("missing_evidence:", data.get("missing_evidence"))
    print("forbidden_misclassification_violations:", data.get("forbidden_misclassification_violations"))

    readouts = data.get("state_specific_readouts", {})
    for key in ["state_A", "state_B"]:
        row = readouts.get(key, {}) if isinstance(readouts, dict) else {}
        print(f"\n{key} readout:")
        print("  pdb_id:", row.get("pdb_id"))
        print("  expected_context:", row.get("expected_context"))
        print("  state_role_bucket:", row.get("state_role_bucket"))
        print("  state_role_evidence_found:", row.get("state_role_evidence_found"))
        print("  model_count:", row.get("model_count"))
        print("  chain_ca_counts:", row.get("chain_ca_counts"))
        iface = row.get("chain_interface_readout", {})
        print("  interface_pairs_present:", iface.get("interface_pairs_present"))
        print("  contact_probe_policy:", iface.get("contact_probe_policy"))

    print("\nState separation guard:", data.get("state_separation_guard"))
    print("next_recommended_panel:", data.get("next_recommended_panel"))


if __name__ == "__main__":
    main()
