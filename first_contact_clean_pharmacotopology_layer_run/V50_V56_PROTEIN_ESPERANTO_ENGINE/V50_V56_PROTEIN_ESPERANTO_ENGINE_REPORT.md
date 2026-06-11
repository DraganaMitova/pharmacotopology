# V50-V56 Protein Esperanto Engine Suite

Status: `V50_TO_V56_PROTEIN_ESPERANTO_ENGINE_PASSED_REVIEW_REQUIRED`
folding_problem_solved: `False`
atomistic_md_executed: `False`
Targets: `5`
Supported validations: `5` / `5`
Controls passed: `14` / `14`

## Version Gates
- V50 grammar extraction: `True`
- V51 engine contract: `True`
- V52 coarse simulator MVP: `True`
- V53 hard-class battery: `True`
- V54 perturbation engine: `True`
- V55 post-seal validation: `True`
- V56 OpenMM bridge spec: `True`

## Hard-Class Battery
- `V44_FUS_LC` `intrinsic_disorder_phase_separation` score `supported` basins `{'expanded_disordered': 0.62, 'phase_prone_dynamic': 0.58, 'compact_single_fold': 0.04}`
- `V45_TDP43_LCD` `intrinsic_disorder_phase_separation` score `supported` basins `{'expanded_disordered': 0.62, 'phase_prone_dynamic': 0.58, 'compact_single_fold': 0.04}`
- `V46_CFTR_F508DEL` `membrane_multidomain_folding_proteostasis` score `supported` basins `{'mature_membrane_routed': 0.674414, 'qc_retained_misfolded': 0.0, 'partial_nbd1_rescue': 0.704}`
- `V47_RFAH_CTD` `metamorphic_fold_switching` score `supported` basins `{'alpha_context_basin': 0.4336, 'beta_released_basin': 0.4264, 'averaged_single_fold': 0.0684}`
- `V48_SARS2_ORF6` `short_region_host_interface_hijacking` score `supported` basins `{'host_interface_engaged': 0.941, 'exposed_short_region': 0.62, 'compact_single_fold': 0.04}`

## Controls
- `random_sequence_control`: `True`
- `shuffled_sequence_control`: `True`
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
- `hard_class_validation_control`: `True`
- `folding_problem_solved_never_true`: `True`

## Claim Boundary
We extracted a reusable operator grammar from multiple hard protein regimes and built a leakage-controlled coarse simulation engine that predicts mechanism trajectories and perturbation directions before holdout validation.
