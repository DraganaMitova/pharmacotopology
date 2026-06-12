# V50-V56 Protein Esperanto Engine Suite

Status: `V50_TO_V56_PROTEIN_ESPERANTO_ENGINE_FAILED`
folding_problem_solved: `False`
atomistic_md_performed: `False`
Targets: `5`
Supported validations: `2` / `5`
Controls passed: `12` / `14`

## Version Gates
- V50 grammar extraction: `True`
- V51 engine contract: `True`
- V52 coarse operator-state propagator MVP: `True`
- V53 hard-class battery: `False`
- V54 perturbation engine: `True`
- V55 post-seal validation: `True`
- V56 OpenMM bridge spec: `True`

## Hard-Class Battery
- `V44_FUS_LC` `intrinsic_disorder_phase_separation` score `supported` basins `{'expanded_disordered': 0.582419, 'phase_prone_dynamic': 0.582419, 'compact_single_fold': 0.0}`
- `V45_TDP43_LCD` `intrinsic_disorder_phase_separation` score `supported` basins `{'expanded_disordered': 0.371623, 'phase_prone_dynamic': 0.371623, 'compact_single_fold': 0.0}`
- `V46_CFTR_F508DEL` `membrane_multidomain_folding_proteostasis` score `contradicted` basins `{'mature_membrane_routed': 0.363838, 'qc_retained_misfolded': 0.0, 'partial_nbd1_rescue': 0.39825}`
- `V47_RFAH_CTD` `metamorphic_fold_switching` score `contradicted` basins `{'alpha_context_basin': 0.054797, 'beta_released_basin': 0.263084, 'averaged_single_fold': 0.606175}`
- `V48_SARS2_ORF6` `short_region_host_interface_hijacking` score `contradicted` basins `{'host_interface_engaged': 0.194225, 'exposed_short_region': 0.100433, 'compact_single_fold': 0.0}`

## Controls
- `random_sequence_control`: `True`
- `shuffled_sequence_control`: `False`
- `swapped_evidence_control`: `True`
- `wrong_target_control`: `True`
- `generic_annotation_only_control`: `True`
- `coordinate_leakage_control`: `True`
- `internal_runtime_leakage_control`: `True`
- `spatial_proxy_tagging_control`: `True`
- `forced_wrong_grammar_control`: `True`
- `failed_prediction_not_repaired_after_holdout`: `True`
- `holdout_opened_before_seal_control`: `True`
- `wild_type_mutant_direction_control`: `True`
- `hard_class_validation_control`: `False`
- `folding_problem_solved_never_true`: `True`

## Claim Boundary
No claim allowed until failed checks are fixed.
