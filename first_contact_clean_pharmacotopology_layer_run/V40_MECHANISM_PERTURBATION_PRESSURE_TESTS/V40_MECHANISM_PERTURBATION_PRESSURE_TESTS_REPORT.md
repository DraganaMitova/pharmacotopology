# V40 Mechanism Perturbation Pressure Tests

Status: `V40_MECHANISM_PERTURBATION_PRESSURE_PASSED_CLAIM_DISABLED`
Supported targets: `3` / `3`
Controls: `14` / `14`
Coordinate-derived source count: `0`
Internal runtime source count: `0`
Answer key used for prediction: `False`

## Target Summary
### KcsA
Validation level: `perturbation_supported_by_holdout`
Perturbations: `4`
Supported perturbation buckets: `['filter_signature_perturbation', 'ion_selectivity_perturbation', 'membrane_topology_context_perturbation', 'oligomer_interface_context_perturbation']`
Scientist answer: Filter and ion-selectivity perturbations are causal pressure points; generic channel annotation alone is insufficient. Supported perturbation buckets: filter_signature_perturbation, ion_selectivity_perturbation, membrane_topology_context_perturbation, oligomer_interface_context_perturbation.

### XCL1_lymphotactin
Validation level: `perturbation_supported_by_holdout`
Perturbations: `5`
Supported perturbation buckets: `['mixed_state_pooling_error', 'state_A_loss_or_weakening', 'state_B_loss_or_weakening', 'state_function_decoupling', 'two_state_balance_shift']`
Scientist answer: State separation is causal pressure; deleting either state or pooling states weakens or invalidates the metamorphic mechanism. Supported perturbation buckets: mixed_state_pooling_error, state_A_loss_or_weakening, state_B_loss_or_weakening, state_function_decoupling, two_state_balance_shift.

### alpha_synuclein_SNCA
Validation level: `perturbation_supported_by_holdout`
Perturbations: `5`
Supported perturbation buckets: `['aggregation_context_overpromotion', 'context_bound_helix_overpromotion', 'disorder_evidence_loss', 'disorder_to_order_context_shift', 'ensemble_to_single_fold_forcing']`
Scientist answer: Disorder evidence and context boundaries are causal pressure; compact single-fold and native-fold overpromotion are blocked. Supported perturbation buckets: aggregation_context_overpromotion, context_bound_helix_overpromotion, disorder_evidence_loss, disorder_to_order_context_shift, ensemble_to_single_fold_forcing.

## Plain English Interpretation
V40 passing means the non-coordinate mechanism grammar identifies causal perturbation pressure for KcsA, XCL1, and SNCA: which operators should break, weaken, shift, or block the mechanism, while keeping coordinates, MD, answer-key leakage, and folding-solved claims disabled.
