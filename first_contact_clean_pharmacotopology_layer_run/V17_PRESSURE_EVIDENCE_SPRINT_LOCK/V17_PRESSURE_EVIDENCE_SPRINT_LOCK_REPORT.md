# V17 Pressure Evidence Sprint Lock

One batch: evidence manifest + evidence preflight + zero-MD evidence-readiness readout + one V18 target decision.

Status: `V17_PRESSURE_EVIDENCE_SPRINT_LOCKED`
Selected V18 target: `p53_TAD_MDM2`
Selected V18 test: `p53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST`
Claim allowed: `False`
New MD executed: `False`
Positive folding evidence targets: `[]`

## Interpretation
V17 is a one-shot pressure-evidence sprint lock. It combines evidence manifest, file preflight, zero-MD evidence-readiness readout, and one selected V18 target decision. It does not run MD, tune thresholds, use native metrics for selection, or upgrade role classification into folding evidence.

## Target rows

### p53_TAD_MDM2
Role class: `disorder_partner_induced_object`
Zero-MD status: `partner_bound_interface_evidence_context_available_zero_md`
Ready for V18: `True`
Next evidence test: `p53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST`

Available evidence:
- `interface_or_contact_evidence`
- `leakage_guard_autonomous_fold_vs_partner_induced_fold`
- `partner_bound_complex_context`
- `state_label_bound_context`

Missing evidence:
- `isolated_TAD_disorder_or_clean_abstain_context`
- `external_couplings_if_available`

### KcsA
Role class: `membrane_pore_oligomer_object`
Zero-MD status: `membrane_pore_context_only_annotations_missing_zero_md`
Ready for V18: `False`
Next evidence test: `KcsA_MEMBRANE_PORE_ROLE_EVIDENCE_READOUT`

Available evidence:
- `KcsA_complex_coordinate_context`
- `chain_and_interface_identity`
- `leakage_guard_against_soluble_core_misread`

Missing evidence:
- `membrane_topology_annotation`
- `pore_selectivity_filter_annotation`
- `oligomer_or_tetramer_assembly_context`
- `transmembrane_helix_roles`
- `external_couplings_if_available`

### XCL1_lymphotactin
Role class: `metamorphic_switch_object`
Zero-MD status: `two_state_context_only_state_specific_evidence_guard_missing_zero_md`
Ready for V18: `False`
Next evidence test: `XCL1_STATE_SPECIFIC_ROLE_SEPARATION_READOUT`

Available evidence:
- `state_A_and_state_B_context_structures`

Missing evidence:
- `state_A_and_state_B_labels`
- `monomer_dimer_context`
- `condition_labels_if_available`
- `state_specific_structural_or_contact_evidence`
- `state_specific_external_couplings_or_constraints_if_available`
- `leakage_guard_preventing_mixed_state_fake_core`
