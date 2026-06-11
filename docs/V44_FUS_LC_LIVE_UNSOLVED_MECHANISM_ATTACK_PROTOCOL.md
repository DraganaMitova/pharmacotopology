# V44 FUS-LC Live Unsolved Mechanism Attack Protocol

## Purpose

V44 stops the retrospective benchmark loop and targets a live hard protein-mechanism problem: human FUS low-complexity domain, especially FUS-LC residues 1-214.

This is not a single-structure prediction task. The required output is a sealed de novo mechanism-language solution packet covering ensemble grammar, operator regions, perturbation predictions, condition shifts, experiments, and falsification criteria.

## Target

- Target: human FUS low-complexity domain / FUS-LC
- Accession: UniProt `P35637`
- Region: residues `1-214`
- Expected mechanism class: `intrinsic_disorder_phase_separation_contextual_ensemble`

## Allowed Prediction Inputs

- UniProt sequence and non-coordinate feature annotations
- DisProt disorder-only annotations
- low-complexity and compositional sequence features
- sequence-derived sticker/spacer patterning
- phosphorylation site annotations
- non-coordinate state/function annotations kept separate from validation holdouts

## Blocked Prediction Inputs

- PDB/mmCIF coordinates before sealing
- coordinate-derived contacts or native contact maps
- AlphaFold/ESMFold/RoseTTAFold coordinates before sealing
- internal runtime reports as biological evidence
- validation holdouts before prediction sealing
- answer key/class labels
- target-name-only assignment

## Sealing

The runner writes:

- `sealed_prediction_packet.json`
- `prediction_hash`
- `prediction_timestamp`
- `prediction_inputs_manifest.json`
- `blocked_inputs_manifest.json`
- `no_holdout_access_before_hash: true`

Only after the prediction packet hash exists does V44 write post-seal holdouts and validation scores.

## Required Solution Packet Fields

- `target_id`
- `sequence_region_scope`
- `mechanism_class`
- `operator_regions`
- `operator_region_rationale`
- `predicted_state_grammar`
- `predicted_context_switches`
- `perturbation_predictions`
- `low_resolution_ensemble_prediction`
- `proposed_experimental_tests`
- `falsification_criteria`
- `claim_boundary`

## Operator Buckets

Required buckets include:

- `low_complexity_disorder_operator`
- `aromatic_sticker_pattern_operator`
- `polar_spacer_context_operator`
- `phosphorylation_shift_operator`
- `LLPS_self_association_operator`
- `fibril_or_gel_state_shift_operator`
- `context_bound_RNA_or_protein_interaction_operator`

Forbidden buckets:

- `compact_single_native_fold_operator`
- `solved_atomic_structure_operator`
- `de_novo_global_coordinate_solution_operator`
- `AlphaFold_confidence_proxy_operator`
- `coordinate_contact_operator`

## Pass Conditions

V44 passes only if:

- prediction is sealed before holdout validation
- no coordinate leakage occurs before sealing
- no internal runtime evidence is promoted as biology
- mechanism class is `intrinsic_disorder_phase_separation_contextual_ensemble`
- at least six operator buckets are produced
- at least eight perturbation predictions are produced
- at least five perturbation predictions are supported or partially supported by holdouts
- compact single-fold grammar is rejected
- contradicted prediction count is at most two
- all leakage controls pass
- claim boundary remains honest

V44 may set:

- `live_unsolved_target_solution_packet = true`
- `protein_folding_solved_candidate_strengthened = true`

only if pass conditions are met.

V44 must never set:

- `folding_problem_solved = true`

## Safe Claim Boundary

Allowed if passed:

`We have a sealed live solution packet for the FUS-LC mechanism-language problem: a source-separated ensemble/phase/perturbation model for residues 1-214, supported by post-seal holdout evidence. This is not a universal protein-folding solved claim.`

Forbidden:

- we solved the universal protein-folding problem
- FUS-LC has one solved atomic native structure
- coordinates were predicted de novo for all FUS-LC states
- external review is unnecessary
- all IDP phase-separation mechanisms are solved
