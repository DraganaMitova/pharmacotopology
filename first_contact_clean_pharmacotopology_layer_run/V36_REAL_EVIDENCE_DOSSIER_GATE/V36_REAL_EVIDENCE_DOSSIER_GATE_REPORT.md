# V36 Real Evidence Dossier Gate

Status: `V36_REAL_EVIDENCE_DOSSIERS_READY_CLAIM_DISABLED`
Targets ready: `['KcsA', 'XCL1_lymphotactin', 'alpha_synuclein_SNCA']`
Targets partial: `[]`
Targets blocked: `[]`
External source count: `12`
Coordinate-derived source count: `0`
Internal runtime source count: `0`
Placeholder source count: `0`
Controls: `9` / `9`
Claim allowed: `False`
New MD allowed: `False`
Folding solved: `False`
Next action: `use_V36_dossiers_to_design_target_specific_operator_grammar_no_claim_no_MD`

## Source Counts By Target
- `KcsA`: `4`
- `XCL1_lymphotactin`: `4`
- `alpha_synuclein_SNCA`: `4`

## Source Types By Target
- `KcsA`: `['UniProt sequence/features/function annotations', 'external sequence conservation signatures']`
- `XCL1_lymphotactin`: `['UniProt sequence/features/function annotations', 'literature-derived state/function annotations']`
- `alpha_synuclein_SNCA`: `['DisProt disorder annotations', 'UniProt sequence/features/function annotations', 'experimentally described motifs/state labels', 'literature-derived state/function annotations']`

## Failed Checks
- None

## Controls
### kcsa_v33_v34_coordinate_csvs_blocked
Passed: `True`
Observed status: `V36_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE`
Expected status: `V36_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE`
Reason: KcsA V33/V34 coordinate-derived contact CSVs must be blocked as V36 evidence.

### kcsa_generic_channel_annotation_partial
Passed: `True`
Observed status: `V36_PARTIAL_EVIDENCE_DOSSIERS_CLEAN_ABSTAIN`
Expected status: `V36_PARTIAL_EVIDENCE_DOSSIERS_CLEAN_ABSTAIN`
Reason: A generic channel label is not enough for KcsA readiness.

### xcl1_one_state_only_partial
Passed: `True`
Observed status: `V36_PARTIAL_EVIDENCE_DOSSIERS_CLEAN_ABSTAIN`
Expected status: `V36_PARTIAL_EVIDENCE_DOSSIERS_CLEAN_ABSTAIN`
Reason: XCL1 needs both native-state function contexts.

### xcl1_mixed_state_pooling_blocked
Passed: `True`
Observed status: `V36_BLOCKED_PLACEHOLDER_OR_UNTRUSTED_SOURCES`
Expected status: `V36_BLOCKED_PLACEHOLDER_OR_UNTRUSTED_SOURCES`
Reason: XCL1 state-A/state-B pooling is a blocked grammar error.

### snca_single_fold_grammar_blocked
Passed: `True`
Observed status: `V36_BLOCKED_PLACEHOLDER_OR_UNTRUSTED_SOURCES`
Expected status: `V36_BLOCKED_PLACEHOLDER_OR_UNTRUSTED_SOURCES`
Reason: SNCA must remain an IDP/ensemble problem, not a forced single-fold problem.

### snca_without_disorder_evidence_partial
Passed: `True`
Observed status: `V36_PARTIAL_EVIDENCE_DOSSIERS_CLEAN_ABSTAIN`
Expected status: `V36_PARTIAL_EVIDENCE_DOSSIERS_CLEAN_ABSTAIN`
Reason: SNCA cannot be ready if disorder evidence is removed.

### placeholder_citation_blocked
Passed: `True`
Observed status: `V36_BLOCKED_PLACEHOLDER_OR_UNTRUSTED_SOURCES`
Expected status: `V36_BLOCKED_PLACEHOLDER_OR_UNTRUSTED_SOURCES`
Reason: Placeholder source names or citations must be blocked.

### internal_runtime_source_blocked
Passed: `True`
Observed status: `V36_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE`
Expected status: `V36_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE`
Reason: Runtime reports under first_contact_clean_pharmacotopology_layer_run are not evidence.

### native_coordinate_metrics_before_selection_blocked
Passed: `True`
Observed status: `V36_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE`
Expected status: `V36_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE`
Reason: Native coordinate metrics before selection must be blocked.

## Locked Interpretation
V36 readiness means the repo now has real external, non-coordinate dossiers for KcsA, XCL1, and SNCA operator grammar. It is not a folding prediction win, does not run MD, does not use native coordinates for selection, and does not solve protein folding.
