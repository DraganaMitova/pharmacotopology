# V21 Pressure Evidence Panel Summary

Summary status: `V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_LOCKED`
Positive pressure evidence targets: `['p53_TAD_MDM2', 'KcsA', 'XCL1_lymphotactin']`
Positive folding evidence targets: `[]`
Claim allowed: `False`
New MD executed: `False`
Fixed residue cutoff used: `False`

## Locked interpretation
V21 summarizes three zero-MD pressure-context evidence wins: p53/MDM2 disorder-partner pressure, KcsA membrane/pore pressure, and XCL1 metamorphic state-specific pressure. This is not a universal folding claim; positive_folding_evidence_targets remains empty and claim_allowed remains false.

## Target rows
### p53_TAD_MDM2
- Pressure class: `disorder_partner_induced_binding`
- Evidence status: `passed`
- Summary: isolated TAD clean-abstains while MDM2-bound TAD preserves partner-induced interface/helix evidence
- Missing evidence: `['external_couplings_if_available', 'isolated_TAD_disorder_or_clean_abstain_context', 'isolated_TAD_disorder_prior_or_external_annotation']`
- MD readiness: `no_new_MD_recommended_before_missing_external_context`

### KcsA
- Pressure class: `membrane_pore_oligomer_environment`
- Evidence status: `passed`
- Summary: pore/filter diagnostic, potassium-ion shell, TM helix scaffold, and chain/interface context detected without soluble-core misclassification
- Missing evidence: `['biological_tetramer_assembly_annotation_for_tetramer_claim', 'external_couplings_if_available']`
- MD readiness: `no_membrane_MD_recommended_before_assembly_and_external_coupling_context`

### XCL1_lymphotactin
- Pressure class: `metamorphic_fold_switching`
- Evidence status: `passed`
- Summary: state A and state B role evidence kept in separate buckets with no mixed-state fake core or forced single-fold claim
- Missing evidence: `['condition_labels_if_available', 'state_specific_external_couplings_or_constraints_if_available']`
- MD readiness: `no_MD_recommended_before_state_conditions_and_external_constraints`

## Next panel
`V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL`

## Failed checks
`[]`
