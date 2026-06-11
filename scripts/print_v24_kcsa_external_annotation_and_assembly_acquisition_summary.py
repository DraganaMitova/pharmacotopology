#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION" / "v24_kcsa_external_annotation_and_assembly_acquisition_certificate.json"


def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V24 KcsA EXTERNAL ANNOTATION AND ASSEMBLY ACQUISITION ===")
    for key in [
        "run_mode",
        "test_status",
        "target_id",
        "target_role",
        "positive_pressure_evidence_found",
        "positive_folding_evidence_found",
        "membrane_pore_role_evidence_found",
        "claim_allowed",
        "new_md_executed",
        "membrane_md_executed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "pore_filter_annotation_present",
        "TVGYG_motif_detected",
        "K_plus_ion_diagnostic_shell_present",
        "transmembrane_role_present",
        "membrane_environment_context_present",
        "oligomer_or_chain_interface_context_present",
        "tetramer_context_annotation_present",
        "biological_assembly_context_locked",
        "tetramer_claim_allowed",
        "whole_channel_fold_claim_allowed",
        "tetramer_claim_made",
        "whole_channel_fold_claim_made",
        "soluble_core_misclassification_avoided",
        "external_couplings_present",
        "external_couplings_required_for_V24",
        "external_couplings_required_for_future_MD_or_contact_test",
    ]:
        print(f"{key}: {data.get(key)}")
    print("available_evidence:", data.get("available_evidence"))
    print("missing_evidence:", data.get("missing_evidence"))
    print("forbidden_misclassification_violations:", data.get("forbidden_misclassification_violations"))
    print("v24_failed_checks:", data.get("v24_failed_checks"))

    print("\nPore/filter annotation:")
    pore = data.get("pore_filter_annotation", {})
    for key in ["selectivity_filter_role", "TVGYG_motif_detected", "K_plus_ion_diagnostic_shell_present", "potassium_ion_count", "multi_radius_report_only", "selection_threshold_used"]:
        print(f"  {key}: {pore.get(key)}")

    print("\nMembrane topology annotation:")
    topo = data.get("membrane_topology_annotation", {})
    for key in ["external_membrane_topology_reference_locked", "transmembrane_helix_scaffold", "membrane_environment_context_present", "soluble_core_misclassification_avoided"]:
        print(f"  {key}: {topo.get(key)}")

    print("\nAssembly preflight:")
    assm = data.get("assembly_preflight", {})
    for key in ["tetramer_context_annotation_present", "biological_assembly_context_locked", "biological_assembly_source_reference_locked", "oligomer_or_chain_interface_context_present", "tetramer_claim_allowed", "whole_channel_fold_claim_allowed", "missing_for_tetramer_claim"]:
        print(f"  {key}: {assm.get(key)}")

    print("\nExternal coupling availability:")
    coup = data.get("external_coupling_availability", {})
    for key in ["external_couplings_present", "external_coupling_or_msa_files", "external_couplings_required_for_V24", "external_couplings_required_for_future_MD_or_contact_test", "coupling_policy"]:
        print(f"  {key}: {coup.get(key)}")

    print("\nNext decision:")
    dec = data.get("next_decision", {})
    for key in ["decision_status", "selected_next_panel", "reason", "new_md_recommended", "membrane_md_recommended", "parallel_md_paths_allowed"]:
        print(f"  {key}: {dec.get(key)}")


if __name__ == "__main__":
    main()
