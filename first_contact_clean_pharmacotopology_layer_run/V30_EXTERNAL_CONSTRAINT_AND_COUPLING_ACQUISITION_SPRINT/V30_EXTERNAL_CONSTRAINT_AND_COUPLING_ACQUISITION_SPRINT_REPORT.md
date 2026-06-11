# V30 External Constraint and Coupling Acquisition Sprint

Status: `V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT_LOCKED`
Selected next panel: `V31_CONSTRAINT_BACKED_OPERATOR_READOUT_PREFLIGHT`
Selected V31 targets: `['XCL1_lymphotactin', 'KcsA']`
Constraint/coupling-ready targets: `['XCL1_lymphotactin', 'KcsA', 'p53_TAD_MDM2', '4AKE', '1UBQ', '1CLL']`
MD-ready targets: `[]`
Claim allowed: `False`
New MD allowed: `False`

## Locked interpretation
V30 converts the V29 MD-readiness block into a practical acquisition sprint: local target-specific couplings/constraints are scanned, missing constraint layers are made explicit, and XCL1/KcsA are selected for the next external-constraint import/preflight step. No MD or folding claim is allowed.

## Target rows
### XCL1_lymphotactin
Acquisition status: `local_constraint_or_coupling_files_present_ready_for_next_readout_preflight`
Coupling/constraint files: `['first_contact_clean_pharmacotopology_layer_run/V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL/v22_coupling_availability_scan.json', 'first_contact_clean_pharmacotopology_layer_run/V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION/V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION_REPORT.md', 'first_contact_clean_pharmacotopology_layer_run/V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION/v27_xcl1_condition_and_coupling_evidence_acquisition_certificate.json', 'first_contact_clean_pharmacotopology_layer_run/V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION/v27_xcl1_condition_label_manifest.json', 'first_contact_clean_pharmacotopology_layer_run/V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION/v27_xcl1_condition_preflight.json', 'first_contact_clean_pharmacotopology_layer_run/V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION/v27_xcl1_next_decision.json', 'first_contact_clean_pharmacotopology_layer_run/V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION/v27_xcl1_state_specific_coupling_availability.json']`
Missing: `['state_specific_external_couplings_or_constraints_if_available', 'independent_condition_or_oligomerization_context_if_available', 'mixed_state_leakage_guard_preserved']`
Recommended next action: `acquire_or_import_state_specific_external_constraints_then_preflight_no_MD`

### KcsA
Acquisition status: `local_constraint_or_coupling_files_present_ready_for_next_readout_preflight`
Coupling/constraint files: `['first_contact_clean_pharmacotopology_layer_run/V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL/v22_coupling_availability_scan.json', 'first_contact_clean_pharmacotopology_layer_run/V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION/v24_kcsa_external_coupling_availability.json', 'first_contact_clean_pharmacotopology_layer_run/V25_FAST_MECHANISM_EVIDENCE_SPRINT/v25_kcsa_coupling_interface_readout.json']`
Missing: `['external_couplings_if_available', 'pore_filter_coupling_support', 'expanded_biological_assembly_coordinate_or_author_defined_assembly_model_for_tetramer_claim']`
Recommended next action: `acquire_or_import_KcsA_external_couplings_and_assembly_interface_constraints_no_MD`

### p53_TAD_MDM2
Acquisition status: `local_constraint_or_coupling_files_present_ready_for_next_readout_preflight`
Coupling/constraint files: `['first_contact_clean_pharmacotopology_layer_run/V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL/v22_coupling_availability_scan.json']`
Missing: `['external_couplings_if_available', 'raw_local_external_annotation_files_if_needed']`
Recommended next action: `optional_external_coupling_or_raw_annotation_import_no_MD`

### 4AKE
Acquisition status: `local_constraint_or_coupling_files_present_ready_for_next_readout_preflight`
Coupling/constraint files: `['data/all_locked_real_external_coupling_holdout_manifest_v0.locked.json', 'data/external_coupling_target_manifest_v0.locked.json', 'data/folding_real_coordinate_visual_8_couplings.locked.json', 'data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json', 'external_msa/4ake_pfam00406/4ake_pfam00406_focus.fasta', 'external_msa/4ake_pfam00406/PF00406_full.sto', 'external_msa_free_predictors/4ake.fasta', 'external_msa_free_predictors/4ake_esmfold_api.esmfold_api_report.json', 'external_msa_free_predictors/ESMFOLD2_AND_SPIRED_NEXT.md', 'external_msa_free_predictors/README.md', 'external_msa_free_predictors/run_all_4ake_msa_free_tryhard.sh', 'external_msa_free_predictors/tryhard_runs_live/4ake_esmfold_api.pdb']`
Missing: `['external_validation_if_available', 'interdomain_hinge_or_closure_evidence']`
Recommended next action: `keep_locked_until_target_specific_external_constraints_are_available`

### 1UBQ
Acquisition status: `local_constraint_or_coupling_files_present_ready_for_next_readout_preflight`
Coupling/constraint files: `['data/all_locked_real_external_coupling_holdout_manifest_v0.locked.json', 'data/blind_holdout_manifest_10af_plmc_v0.locked.json', 'data/blind_holdout_manifest_10dc_plmc_v0.locked.json', 'data/blind_holdout_manifest_two_target_plmc_v0.locked.json', 'data/folding_real_coordinate_holdout_1ubq_external_couplings.v0.locked.json', 'data/folding_real_coordinate_holdout_1ubq_query_centered_apc_external_couplings.v0.locked.json']`
Missing: `['DCA_background_enrichment', 'cross_target_validation']`
Recommended next action: `keep_locked_until_target_specific_external_constraints_are_available`

### 1CLL
Acquisition status: `local_constraint_or_coupling_files_present_ready_for_next_readout_preflight`
Coupling/constraint files: `['data/all_locked_real_external_coupling_holdout_manifest_v0.locked.json', 'data/external_coupling_target_manifest_v0.locked.json', 'data/folding_real_coordinate_visual_8_couplings.locked.json', 'data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json']`
Missing: `['interdomain_hinge_evidence_if_claim_needed']`
Recommended next action: `keep_locked_until_target_specific_external_constraints_are_available`

