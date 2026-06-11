# V33 Negative Controls

Status: `V33_NEGATIVE_CONTROLS_PASSED_CLAIM_DISABLED`
Passed controls: `9` / `9`
Claim allowed: `False`
New MD allowed: `False`

## Controls
### baseline_real_kcsa_sources_pass_before_negative_controls
Observed: `V33_CONSTRAINT_BACKED_OPERATOR_READOUT_PASSED_CLAIM_DISABLED`
Expected: `V33_CONSTRAINT_BACKED_OPERATOR_READOUT_PASSED_CLAIM_DISABLED`
Passed: `True`
Reason: Real KcsA pore/filter + assembly/interface imported source rows should pass only as claim-disabled operator readout.

### missing_pore_filter_bucket_abstains
Observed: `V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED`
Expected: `V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED`
Passed: `True`
Reason: KcsA cannot produce a pore/filter operator readout without the pore/filter source bucket.

### missing_assembly_interface_bucket_abstains
Observed: `V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED`
Expected: `V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED`
Passed: `True`
Reason: KcsA cannot produce a tetramer/interface operator readout without the assembly/interface source bucket.

### empty_operator_buckets_abstain
Observed: `V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED`
Expected: `V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED`
Passed: `True`
Reason: No source buckets means no operator readout.

### wrong_target_xcl1_blocked_by_v33_v0_scope
Observed: `V33_BLOCKED_PRECONDITION_OR_PROVENANCE_FAILURE_CLAIM_DISABLED`
Expected: `V33_BLOCKED_PRECONDITION_OR_PROVENANCE_FAILURE_CLAIM_DISABLED`
Passed: `True`
Reason: V33 v0 is deliberately scoped to KcsA after V32 selection; XCL1 needs its own state-specific import/readout.

### internal_runtime_source_poison_blocked_by_v32
Observed: `V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_BLOCKED`
Expected: `V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_BLOCKED`
Passed: `True`
Reason: Files under first_contact_clean_pharmacotopology_layer_run must never become external evidence.

### annotation_only_context_does_not_select_v33
Observed: `V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED`
Expected: `V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED`
Passed: `True`
Reason: Annotation-only rows can give role context, but must not open a constraint-backed operator readout.

### v32_only_pore_bucket_not_ready_for_kcsa_v33
Observed: `V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED`
Expected: `V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED`
Passed: `True`
Reason: KcsA needs both pore/filter and assembly/interface source buckets before V33 selection.

### v32_only_interface_bucket_not_ready_for_kcsa_v33
Observed: `V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED`
Expected: `V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED`
Passed: `True`
Reason: Assembly/interface context alone must not become a KcsA pore/filter readout.

## Locked interpretation
These controls validate V32/V33 evidence gating only. Passing them means the system blocks internal/runtime evidence, annotation-only rows, missing buckets, and wrong-target misuse while keeping claim and MD gates closed. It does not prove de novo folding or universal protein-folding prediction.
