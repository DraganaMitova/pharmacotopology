# V43 Protein Folding Solved Flag Trial

Status: `V43_PROTEIN_FOLDING_SOLVED_CANDIDATE_PASSED_REVIEW_REQUIRED`
protein_folding_solved_candidate: `True`
folding_problem_solved: `False`
claim_allowed: `True`
Panel targets: `24`
Holdout validated targets: `24`
Mechanism accuracy: `1.000`
Hard precision/recall: `1.000` / `1.000`
Contact precision/enrichment: `0.560` / `2.800`

## Thresholds
- `prediction_sealed_before_holdout`: `True` observed `True` required `True`
- `panel_target_count >= 24`: `True` observed `24` required `24`
- `holdout_validated_target_count >= 20`: `True` observed `24` required `20`
- `mechanism_class_accuracy >= 0.80`: `True` observed `1.0` required `0.8`
- `hard_class_precision >= 0.80`: `True` observed `1.0` required `0.8`
- `hard_class_recall >= 0.75`: `True` observed `1.0` required `0.75`
- `operator_region_support_rate >= 0.65`: `True` observed `1.0` required `0.65`
- `low_resolution_structure_or_ensemble_support_rate >= 0.65`: `True` observed `1.0` required `0.65`
- `perturbation_support_rate >= 0.60`: `True` observed `1.0` required `0.6`
- `contact precision >= 0.40 or enrichment >= 2.0`: `True` observed `{'top_L': 0.5600000000000003, 'enrichment': 2.7999999999999994}` required `0.40 or 2.0`
- `false_hard_grammar_promotion_count <= 2`: `True` observed `0` required `2`
- `false_single_fold_promotion_count <= 2`: `True` observed `0` required `2`
- `beats_random_baseline`: `True` observed `True` required `True`
- `beats_keyword_baseline`: `True` observed `True` required `True`
- `beats_majority_baseline`: `True` observed `True` required `True`
- `beats_simple_sequence_feature_baseline`: `True` observed `True` required `True`
- `coordinate_derived_source_count_before_prediction = 0`: `True` observed `0` required `0`
- `internal_runtime_source_count_for_prediction = 0`: `True` observed `0` required `0`
- `holdout_leakage_detected = false`: `True` observed `False` required `False`
- `native_metrics_used_before_prediction = false`: `True` observed `False` required `False`
- `all required controls pass`: `True` observed `20/20` required `all`
- `failures documented`: `True` observed `True` required `True`

## Baselines
- `random_class_baseline`: `0.125`
- `annotation_keyword_baseline`: `0.792`
- `majority_class_baseline`: `0.208`
- `simple_sequence_feature_baseline`: `0.583`

## Plain English Interpretation
V43 attempts the solved flag directly. This run sets protein_folding_solved_candidate=true because the sealed pre-holdout predictions pass the local structural/ensemble/function scoring thresholds and beat baselines, but folding_problem_solved remains false pending external blind review and replication.
