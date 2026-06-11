# V35 Non-Coordinate Evolutionary Holdout

Status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Source status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Selected target: `KcsA`
Evidence rows: `0`
Non-coordinate external sources: `0`
Coordinate-derived sources: `0`
Internal runtime sources: `0`
Operator candidate found: `False`
Controls: `9` / `9`
Claim allowed: `False`
New MD allowed: `False`
New MD executed: `False`
Positive folding evidence found: `False`
Folding solved: `False`
Next action: `acquire_real_noncoordinate_external_evolutionary_source_before_any_claim_or_MD`

## Selected Sources
- None

## Failed Checks
- None

## Controls
### missing_noncoordinate_source_clean_abstains
Passed: `True`
Observed status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Expected status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Observed failed checks: `[]`
Reason: No local non-coordinate external evidence must produce a clean abstain, not a fake pass.

### coordinate_derived_kcsa_csv_is_blocked
Passed: `True`
Observed status: `V35_BLOCKED_COORDINATE_DERIVED_SOURCE_SUPPLIED`
Expected status: `V35_BLOCKED_COORDINATE_DERIVED_SOURCE_SUPPLIED`
Observed failed checks: `['row_0:coordinate_derived_field_not_false', 'row_0:coordinate_derived_source_supplied', 'row_0:excluded_evidence_type:coordinate_derived_contact', 'row_0:source_path_not_in_kcsa_noncoordinate_evolutionary_boundary']`
Reason: Coordinate-derived KcsA CSVs from V33/V34 must not open V35.

### internal_runtime_source_is_blocked
Passed: `True`
Observed status: `V35_BLOCKED_INTERNAL_RUNTIME_SOURCE_SUPPLIED`
Expected status: `V35_BLOCKED_INTERNAL_RUNTIME_SOURCE_SUPPLIED`
Observed failed checks: `['row_0:excluded_evidence_type:internal_runtime_report', 'row_0:internal_runtime_source_supplied', 'row_0:source_path_not_in_kcsa_noncoordinate_evolutionary_boundary']`
Reason: Generated runtime artifacts are audit outputs, not external evidence.

### placeholder_source_name_or_citation_is_invalid
Passed: `True`
Observed status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Expected status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Observed failed checks: `['row_0:placeholder_or_missing_source_name', 'row_0:placeholder_or_missing_source_url_or_citation']`
Reason: A row without stable provenance must not open V35 even if its content looks plausible.

### annotation_only_source_does_not_open_v35
Passed: `True`
Observed status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Expected status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Observed failed checks: `['row_0:excluded_evidence_type:annotation_only_claim']`
Reason: Annotation-only rows remain context-only and cannot become non-coordinate evidence.

### wrong_target_in_kcsa_v35_path_does_not_open
Passed: `True`
Observed status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Expected status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Observed failed checks: `['row_0:target_not_KcsA']`
Reason: V35 v0 is scoped to KcsA; a wrong target in the KcsA V35 boundary must not pass.

### renamed_potassium_filter_signature_fails_content_validation
Passed: `True`
Observed status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Expected status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Observed failed checks: `['row_0:missing_potassium_channel_filter_signature']`
Reason: A shuffled or renamed filter signature must fail sequence-content validation.

### ion_specificity_relabel_away_from_potassium_fails
Passed: `True`
Observed status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Expected status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Observed failed checks: `['row_0:missing_potassium_ion_specificity']`
Reason: KcsA-like non-coordinate context must not survive relabeling away from K+.

### generic_channel_annotation_is_not_enough
Passed: `True`
Observed status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Expected status: `V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE`
Observed failed checks: `['row_0:missing_external_evolutionary_or_family_basis', 'row_0:missing_potassium_channel_family_specificity', 'row_0:missing_potassium_channel_filter_signature', 'row_0:missing_potassium_ion_specificity']`
Reason: A generic channel label is not enough without filter, ion, and evolutionary/family specificity.

## Locked Interpretation
V35 only opens on real local non-coordinate external evolutionary evidence with stable provenance and KcsA-like potassium-channel sequence/family content. Coordinate-derived contacts, internal runtime reports, annotation-only rows, placeholders, wrong targets, and generic channel labels remain blocked or clean-abstained. Passing V35 is not a folding claim, does not run MD, and does not solve protein folding.
