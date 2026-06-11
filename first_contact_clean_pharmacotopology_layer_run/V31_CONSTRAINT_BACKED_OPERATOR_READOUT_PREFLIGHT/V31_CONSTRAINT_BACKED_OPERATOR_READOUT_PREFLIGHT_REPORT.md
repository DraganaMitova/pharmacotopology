# V31 Constraint-Backed Operator Readout Preflight

Status: `V31_CLEAN_ABSTAIN_NO_REAL_EXTERNAL_CONSTRAINTS_FOR_SELECTED_TARGETS`
Selected V31 targets: `['XCL1_lymphotactin', 'KcsA']`
Selected V32 target: `None`
Real external constraint targets in scope: `[]`
Claim allowed: `False`
New MD allowed: `False`

## Allowed-use policy
- `real_external_constraint_or_coupling` -> `allowed_for_constraint_backed_operator_readout`
- `real_external_alignment_source` -> `allowed_for_constraint_derivation_preflight_only`
- `annotation_only_external_context` -> `allowed_for_role_context_only`
- `external_structure_source` -> `allowed_for_structure_context_or_validation_only`
- `generated_internal_report` -> `allowed_for_audit_only`
- `model_prediction_or_msa_free_artifact` -> `allowed_for_audit_only`
- `unverified_constraint_like_file` -> `excluded`
- `unverified_annotation_like_file` -> `excluded`
- `unusable_or_unclassified` -> `excluded`

## Locked interpretation
V31 prevents self-confirmation leakage by classifying every V30 candidate file by provenance and allowed use. Generated runtime reports are audit-only; annotation files are role-context-only; only real external constraint/coupling files may support a V32 constraint-backed operator readout. No MD or folding claim is allowed.

## Target rows
### XCL1_lymphotactin
In selected V31 scope: `True`
Preflight status: `no_real_external_constraint_for_constraint_backed_readout`
Real external constraints: `[]`
Generated internal reports: `57`
Allowed-use counts: `{'allowed_for_audit_only': 57}`

### KcsA
In selected V31 scope: `True`
Preflight status: `no_real_external_constraint_for_constraint_backed_readout`
Real external constraints: `[]`
Generated internal reports: `51`
Allowed-use counts: `{'allowed_for_audit_only': 51}`

### p53_TAD_MDM2
In selected V31 scope: `False`
Preflight status: `no_real_external_constraint_for_constraint_backed_readout`
Real external constraints: `[]`
Generated internal reports: `43`
Allowed-use counts: `{'allowed_for_audit_only': 43}`

### 4AKE
In selected V31 scope: `False`
Preflight status: `real_external_constraint_available_for_v32_preflight`
Real external constraints: `['data/all_locked_real_external_coupling_holdout_manifest_v0.locked.json', 'data/external_coupling_target_manifest_v0.locked.json', 'data/folding_real_coordinate_visual_8_couplings.locked.json', 'data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json']`
Generated internal reports: `16`
Allowed-use counts: `{'allowed_for_audit_only': 22, 'allowed_for_constraint_backed_operator_readout': 4, 'allowed_for_constraint_derivation_preflight_only': 2}`

### 1UBQ
In selected V31 scope: `False`
Preflight status: `real_external_constraint_available_for_v32_preflight`
Real external constraints: `['data/all_locked_real_external_coupling_holdout_manifest_v0.locked.json', 'data/blind_holdout_manifest_10af_plmc_v0.locked.json', 'data/blind_holdout_manifest_10dc_plmc_v0.locked.json', 'data/blind_holdout_manifest_two_target_plmc_v0.locked.json', 'data/folding_real_coordinate_holdout_1ubq_external_couplings.v0.locked.json', 'data/folding_real_coordinate_holdout_1ubq_query_centered_apc_external_couplings.v0.locked.json']`
Generated internal reports: `14`
Allowed-use counts: `{'allowed_for_audit_only': 14, 'allowed_for_constraint_backed_operator_readout': 6}`

### 1CLL
In selected V31 scope: `False`
Preflight status: `real_external_constraint_available_for_v32_preflight`
Real external constraints: `['data/all_locked_real_external_coupling_holdout_manifest_v0.locked.json', 'data/external_coupling_target_manifest_v0.locked.json', 'data/folding_real_coordinate_visual_8_couplings.locked.json', 'data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json']`
Generated internal reports: `15`
Allowed-use counts: `{'allowed_for_audit_only': 15, 'allowed_for_constraint_backed_operator_readout': 4}`

