# V37 Mechanism Question Probes

Status: `V37_MECHANISM_QUESTION_PROBES_PASSED_CLAIM_DISABLED`
Input V36 status: `V36_REAL_EVIDENCE_DOSSIERS_READY_CLAIM_DISABLED`
Assigned targets: `3` / `3`
Partial targets: `0`
Blocked targets: `0`
Target-name masking passed: `True`
Swapped dossier detection passed: `True`
Coordinate-derived source count: `0`
Internal runtime source count: `0`
Controls: `12` / `12`
Claim allowed: `False`
New MD allowed: `False`
Folding solved: `False`
Next action: `use_mechanism_question_maps_to_design_class_specific_readouts_no_claim_no_MD`

## Targets
### KcsA
Mechanism class: `membrane_pore_filter_oligomeric_ion_selectivity`
Required operators found: `['filter_signature_operator', 'ion_selectivity_operator', 'membrane_topology_operator', 'oligomer_or_interface_context_operator']`
Forbidden operators triggered: `[]`

Scientific questions:
- `True` Does the non-coordinate dossier support potassium-channel identity?
- `True` Does it support membrane/pore/filter context?
- `True` Does it separate filter/ion-selectivity grammar from whole-channel folding claims?
- `True` Does it avoid using coordinate contacts?

### XCL1_lymphotactin
Mechanism class: `metamorphic_two_state_fold_switch`
Required operators found: `['no_mixed_state_pooling_operator', 'state_A_chemokine_monomer_operator', 'state_B_beta_sandwich_dimer_operator', 'state_specific_function_operator']`
Forbidden operators triggered: `[]`

Scientific questions:
- `True` Does the dossier support two biologically relevant states?
- `True` Does it prevent mixed-state pooling?
- `True` Does it keep state A and state B function evidence separate?
- `True` Does it refuse a single-fold claim?

### alpha_synuclein_SNCA
Mechanism class: `intrinsic_disorder_contextual_ensemble`
Required operators found: `['context_bound_structure_operator', 'disorder_to_order_context_operator', 'ensemble_not_single_fold_operator', 'intrinsic_disorder_operator']`
Forbidden operators triggered: `[]`

Scientific questions:
- `True` Does the dossier support free-state disorder?
- `True` Does it represent membrane-bound helix as context-dependent, not solved folding?
- `True` Does it block single-native-fold grammar?
- `True` Does it keep aggregation/amyloid context from becoming a native-fold claim?

## Controls
- `baseline_v36_dossiers_assign_all_three_mechanism_classes`: `PASS`
- `target_name_masked_dossiers_still_assign_from_evidence_content`: `PASS`
- `swapped_dossiers_are_detected_not_renamed_into_target_grammar`: `PASS`
- `kcsa_without_ion_filter_evidence_is_partial`: `PASS`
- `xcl1_without_state_b_is_partial`: `PASS`
- `xcl1_forced_mixed_state_pooling_is_blocked`: `PASS`
- `snca_without_disorder_evidence_is_partial`: `PASS`
- `snca_forced_compact_single_fold_is_blocked`: `PASS`
- `generic_only_annotations_do_not_assign_high_confidence_grammar`: `PASS`
- `coordinate_derived_source_supplied_to_v37_is_blocked`: `PASS`
- `internal_runtime_source_supplied_to_v37_is_blocked`: `PASS`
- `placeholder_source_row_supplied_to_v37_is_blocked`: `PASS`

## Failed Checks
- None

## Locked Interpretation
V37 passing means the V36 non-coordinate dossiers can classify the folding-problem type: KcsA as membrane pore/filter ion-selectivity grammar, XCL1 as metamorphic state-switch grammar, and SNCA as intrinsic-disorder ensemble grammar. This is not structure prediction, not MD, and not a folding-solved claim.
