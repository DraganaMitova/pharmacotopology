# V16 Pressure Evidence Gap Lock

This is a gap lock, not a new evidence claim.

Status: `V16_PRESSURE_EVIDENCE_GAP_LOCKED`
Claim allowed: `False`
New MD executed: `False`
Role-classification passed targets: `['KcsA', 'XCL1_lymphotactin', 'p53_TAD_MDM2']`
Positive folding-evidence targets: `[]`

## Locked interpretation

V16 passed zero-MD pressure-role transfer as role classification, not folding evidence. The gap lock freezes exactly what is missing before each pressure target can become a real folding/contact evidence test. No MD, threshold tuning, native-metric selection, fixed residue cutoff, or claim upgrade is allowed here.

## Target gaps
### KcsA
- Role class: `membrane_pore_oligomer_object`
- Role transfer status: `membrane_pore_context_detected_soluble_core_misclassification_avoided_no_tetramer_claim`
- Role classification passed: `True`
- Positive folding evidence found: `False`
- Forbidden misclassification avoided: `True`
- Next evidence test: `membrane_pore_role_evidence_readout`
- Missing for evidence test:
  - `membrane_topology_annotation`
  - `pore_selectivity_filter_annotation`
  - `oligomer_or_tetramer_assembly_context`
  - `chain_and_interface_identity`
  - `transmembrane_helix_roles`
  - `external_couplings_if_available`
  - `leakage_guard_against_soluble_core_misread`
- Evidence claim forbidden until:
  - `membrane_topology_annotation_present`
  - `assembly_or_interface_context_present`
  - `pore_filter_residue_context_present`

### XCL1_lymphotactin
- Role class: `metamorphic_switch_object`
- Role transfer status: `metamorphic_two_state_context_detected_single_fold_forcing_avoided_no_switch_claim`
- Role classification passed: `True`
- Positive folding evidence found: `False`
- Forbidden misclassification avoided: `True`
- Next evidence test: `state_specific_role_separation_readout`
- Missing for evidence test:
  - `state_A_and_state_B_labels`
  - `monomer_dimer_context`
  - `condition_labels_if_available`
  - `state_specific_structural_or_contact_evidence`
  - `state_specific_external_couplings_or_constraints_if_available`
  - `leakage_guard_preventing_mixed_state_fake_core`
- Evidence claim forbidden until:
  - `state_specific_labels_present`
  - `state_specific_evidence_present`
  - `mixed_state_leakage_guard_present`

### p53_TAD_MDM2
- Role class: `disorder_partner_induced_object`
- Role transfer status: `partner_induced_complex_role_context_detected_no_autonomous_fold_claim`
- Role classification passed: `True`
- Positive folding evidence found: `False`
- Forbidden misclassification avoided: `True`
- Next evidence test: `partner_induced_interface_or_helix_readout`
- Missing for evidence test:
  - `isolated_TAD_disorder_or_clean_abstain_context`
  - `partner_bound_complex_context`
  - `interface_or_contact_evidence`
  - `state_labels_isolated_vs_MDM2_bound`
  - `external_couplings_if_available`
  - `leakage_guard_autonomous_fold_vs_partner_induced_fold`
- Evidence claim forbidden until:
  - `isolated_TAD_abstain_context_present`
  - `interface_contact_evidence_present`
  - `state_label_guard_present`
