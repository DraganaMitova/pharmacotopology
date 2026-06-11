# V16 Zero-MD Role Transfer Readout

This readout does not run MD, does not tune thresholds, and does not claim solved folding.

Status: `V16_ZERO_MD_ROLE_TRANSFER_READOUT_COMPLETED_CLAIM_DISABLED`
Claim allowed: `False`
New MD executed: `False`
Role-classification passed targets: `['p53_TAD_MDM2', 'KcsA', 'XCL1_lymphotactin']`
Positive folding-evidence targets: `[]`
Clean abstain targets: `[]`

## Target rows
### p53_TAD_MDM2
- Status: `partner_induced_complex_role_context_detected_no_autonomous_fold_claim`
- Selected/abstain: `partner_induced_interface_or_helix_signal_found`
- Role buckets: `['partner_induced_interface_or_helix_signal_found', 'binding_stabilized_interface_core_context']`
- Monitor-only roles: `[]`
- Clean abstain roles: `['clean_abstain_no_autonomous_isolated_TAD_core_claim']`
- Forbidden violations: `[]`
- Limitations: `['isolated_TAD_material_not_required_at_this_stage', 'disorder_prior_and_external_couplings_later', 'no_universal_p53_fold_claim']`

### KcsA
- Status: `membrane_pore_context_detected_soluble_core_misclassification_avoided_no_tetramer_claim`
- Selected/abstain: `membrane_pore_roles_detected_without_soluble_core_misclassification`
- Role buckets: `['membrane_pore_oligomer_context', 'transmembrane_helix_scaffold_context']`
- Monitor-only roles: `['pore_selectivity_filter_core_later_annotation_required', 'tetramer_interface_support_later_biological_assembly_required', 'ion_filter_diagnostic_shell_later_annotation_required']`
- Clean abstain roles: `[]`
- Forbidden violations: `[]`
- Limitations: `['biological_assembly_or_oligomer_annotation_later', 'membrane_topology_annotation_later', 'pore_filter_residue_context_later', 'no_heavy_membrane_MD_executed', 'no_whole_fold_claim']`

### XCL1_lymphotactin
- Status: `metamorphic_two_state_context_detected_single_fold_forcing_avoided_no_switch_claim`
- Selected/abstain: `multiple_state_roles_supported_or_clean_state_specific_abstain`
- Role buckets: `['state_A_chemokine_monomer_support_context', 'state_B_alternative_or_beta_sandwich_state_support_context', 'pressure_condition_switch_context']`
- Monitor-only roles: `['shared_local_support_later_coupling_readout_required']`
- Clean abstain roles: `[]`
- Forbidden violations: `[]`
- Limitations: `['condition_or_oligomerization_context_later', 'state_specific_external_couplings_later', 'no_fold_switch_claim_without_state_specific_support']`
