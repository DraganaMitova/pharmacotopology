# V25 Fast Mechanism Evidence Sprint

Status: `V25_FAST_MECHANISM_EVIDENCE_SPRINT_LOCKED`
Positive pressure evidence targets: `['p53_TAD_MDM2', 'KcsA', 'XCL1_lymphotactin']`
Positive folding evidence targets: `[]`
Claim allowed: `False`
New MD executed: `False`

## Locked interpretation
V25 shifts from target-by-target debugging to mechanism sprint mode. It preserves KcsA coupling/interface readiness, XCL1 state-specific separation, and an all-target mechanism-operator table across 4AKE, 1UBQ, 1CLL, p53/MDM2, KcsA, and XCL1. This is evidence for a repeated role/context/operator grammar, not a universal folding claim.

## KcsA coupling/interface
`{'kind': 'V25A_KcsA_COUPLING_INTERFACE_READOUT_v0', 'target_id': 'KcsA', 'run_mode': 'zero_md_coupling_interface_readout_no_simulation_no_threshold_tuning', 'readout_status': 'V25A_KcsA_INTERFACE_CONTEXT_FOUND_COUPLING_MISSING_CLAIM_DISABLED', 'role_class': 'membrane_pore_oligomer_object', 'positive_pressure_evidence_found': True, 'positive_folding_evidence_found': False, 'kcsa_coupling_available': False, 'external_coupling_or_msa_files': [], 'interface_evidence_present': True, 'pore_filter_coupling_support': False, 'pore_filter_role_preserved': True, 'transmembrane_scaffold_preserved': True, 'membrane_context_preserved': True, 'assembly_context_preserved': True, 'tetramer_claim_allowed': False, 'whole_channel_fold_claim_allowed': False, 'tetramer_claim_made': False, 'whole_channel_fold_claim_made': False, 'claim_allowed': False, 'new_md_executed': False, 'membrane_md_executed': False, 'fixed_residue_cutoff_used': False, 'native_metrics_used_for_selection': False, 'selection_threshold_used': False, 'coupling_policy': 'do_not_synthesize_couplings_missing_is_reported_not_backfilled', 'available_evidence': ['pore_filter_role', 'TM_scaffold_context', 'interface_context', 'assembly_context_reference', 'soluble_core_guard'], 'missing_evidence': ['external_couplings_if_available', 'pore_filter_coupling_support'], 'source_v24_status': 'V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION_PASSED_CLAIM_DISABLED'}`

## XCL1 state-specific
`{'kind': 'V25B_XCL1_STATE_SPECIFIC_READOUT_v0', 'target_id': 'XCL1_lymphotactin', 'run_mode': 'zero_md_state_separation_preservation_no_simulation_no_threshold_tuning', 'readout_status': 'V25B_XCL1_STATE_SPECIFIC_EVIDENCE_PRESERVED_CLAIM_DISABLED', 'role_class': 'metamorphic_switch_object', 'positive_pressure_evidence_found': True, 'positive_folding_evidence_found': False, 'state_A_detected': True, 'state_B_detected': True, 'state_specific_role_evidence_found': True, 'mixed_state_pollution': False, 'single_fold_forcing': False, 'single_fold_claim_made': False, 'fold_switch_claim_made': False, 'mixed_state_contact_pooling_used': False, 'claim_allowed': False, 'new_md_executed': False, 'fixed_residue_cutoff_used': False, 'native_metrics_used_for_selection': False, 'selection_threshold_used': False, 'state_separation_policy': 'state_specific_buckets_no_cross_state_pooling_no_single_native_assumption', 'available_evidence': ['leakage_guard_preventing_mixed_state_fake_core', 'monomer_dimer_context', 'state_A_and_state_B_context_structures', 'state_A_and_state_B_labels', 'state_specific_structural_or_contact_evidence'], 'missing_evidence': ['condition_labels_if_available', 'state_specific_external_couplings_or_constraints_if_available'], 'forbidden_misclassification_violations': [], 'source_v20_status': 'V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED'}`

## Next V26 decision
Selected target: `XCL1_lymphotactin`
Selected test: `XCL1_STATE_SEPARATION_OPERATOR_TEST`
Reason: XCL1 has clean state-specific separation with no mixed-state pollution; strongest conceptual mechanism target, but still no MD until condition/coupling evidence is available
