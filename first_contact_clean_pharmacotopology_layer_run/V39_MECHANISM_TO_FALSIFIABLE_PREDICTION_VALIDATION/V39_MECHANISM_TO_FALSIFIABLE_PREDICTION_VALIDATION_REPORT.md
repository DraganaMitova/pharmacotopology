# V39 Mechanism To Falsifiable Prediction Validation

Status: `V39_MECHANISM_PREDICTIONS_VALIDATED_CLAIM_DISABLED`
Validated targets: `3` / `3`
Controls: `12` / `12`
Coordinate-derived source count: `0`
Internal runtime source count: `0`
Answer key used for prediction: `False`

## Target Summary
### KcsA
Validation level: `supported_by_independent_holdout`
Prediction count: `4`
Holdout sources: `4`
Supported buckets: `['filter_or_signature_holdout', 'membrane_topology_holdout', 'oligomer_or_interface_holdout', 'perturbation_holdout']`
Missing buckets: `[]`

### XCL1_lymphotactin
Validation level: `supported_by_independent_holdout`
Prediction count: `5`
Holdout sources: `5`
Supported buckets: `['metamorphic_two_state_holdout', 'perturbation_holdout', 'pooling_rule_holdout', 'state_A_state_B_holdout', 'state_function_holdout']`
Missing buckets: `[]`

### alpha_synuclein_SNCA
Validation level: `supported_by_independent_holdout`
Prediction count: `5`
Holdout sources: `5`
Supported buckets: `['context_bound_holdout', 'disorder_holdout', 'disorder_to_order_holdout', 'perturbation_holdout', 'single_fold_block_holdout']`
Missing buckets: `[]`

## Plain English Interpretation
V39 passing means non-coordinate mechanism grammar produced falsifiable operator-level predictions for KcsA, XCL1, and SNCA, and independent holdout evidence supported them without coordinates, runtime evidence, answer-key leakage, MD, or any folding-solved claim.
