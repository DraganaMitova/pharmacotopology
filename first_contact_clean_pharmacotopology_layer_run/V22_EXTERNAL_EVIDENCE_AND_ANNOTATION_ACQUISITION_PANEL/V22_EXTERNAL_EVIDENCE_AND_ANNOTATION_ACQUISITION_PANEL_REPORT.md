# V22 External Evidence and Annotation Acquisition Panel

Panel status: `V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL_LOCKED`
Selected V23 target: `p53_TAD_MDM2`
Selected V23 test: `p53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST`
Ready for V23: `['p53_TAD_MDM2']`
Positive folding evidence targets: `[]`
Claim allowed: `False`
New MD executed: `False`

## Locked interpretation
V22 acquires and locks external annotation references, checks local annotation/coupling availability, and selects one V23 target. p53/MDM2 is selected because it has the strongest ready external-evidence contrast: isolated disorder reference context plus partner-bound interface evidence. No MD, threshold tuning, native-metric selection, or folding claim is allowed.

## Target preflight rows
### p53_TAD_MDM2
- Ready for V23: `True`
- Annotation source refs locked: `True`
- Raw local annotation files: `[]`
- Missing evidence: `['external_couplings_if_available', 'external_couplings_or_msa_files_if_available', 'isolated_TAD_disorder_or_clean_abstain_context', 'isolated_TAD_disorder_prior_or_external_annotation', 'raw_local_external_annotation_files_if_needed']`

### KcsA
- Ready for V23: `False`
- Annotation source refs locked: `True`
- Raw local annotation files: `[]`
- Missing evidence: `['biological_tetramer_assembly_annotation_for_tetramer_claim', 'external_couplings_if_available', 'external_couplings_or_msa_files_if_available', 'local_biological_assembly_or_OPM_membrane_annotation_file', 'raw_local_external_annotation_files_if_needed']`

### XCL1_lymphotactin
- Ready for V23: `False`
- Annotation source refs locked: `True`
- Raw local annotation files: `[]`
- Missing evidence: `['condition_labels_if_available', 'external_couplings_or_msa_files_if_available', 'local_state_condition_or_monomer_dimer_annotation_file', 'raw_local_external_annotation_files_if_needed', 'state_specific_external_couplings_or_constraints_if_available']`

