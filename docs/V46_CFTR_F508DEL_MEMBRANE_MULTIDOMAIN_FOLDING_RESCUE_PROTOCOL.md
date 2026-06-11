# V46 CFTR F508del Membrane Multidomain Folding Rescue Protocol

## Purpose

V46 leaves the IDP/phase-separation family and attacks a membrane multidomain disease-folding target: human CFTR F508del.

The goal is not atomic structure prediction. V46 asks whether the mechanism language can produce a sealed packet for a coupled membrane, NBD1, interdomain-interface, trafficking/proteostasis, and corrector-rescue defect.

## Target

- Target: human CFTR / cystic fibrosis transmembrane conductance regulator
- Accession: UniProt `P13569`
- Region: full-length CFTR with focus on F508del in NBD1
- Expected mechanism class: `membrane_multidomain_folding_assembly_proteostasis_defect`

## Allowed Prediction Inputs

- sequence
- UniProt annotations
- non-coordinate topology/domain architecture
- mutation/function annotations
- non-coordinate biochemical, trafficking, proteostasis, and family/domain signatures

## Blocked Prediction Inputs

- PDB/mmCIF coordinates before sealing
- PDB-derived contacts or native contact maps
- AlphaFold/ESMFold/RoseTTAFold coordinates
- CFTR coordinate models as prediction evidence
- internal runtime reports as biological evidence
- validation holdouts before prediction sealing
- answer key/class labels

## Sealing

The runner writes:

- `sealed_prediction_packet.json`
- `prediction_hash`
- `prediction_timestamp`
- `prediction_inputs_manifest.json`
- `blocked_inputs_manifest.json`
- `no_holdout_access_before_hash: true`

Only after the hash exists does V46 write post-seal holdouts and validation scores.

## Required Solution Packet Fields

- `target_id`
- `sequence_region_scope`
- `mechanism_class`
- `operator_regions`
- `operator_region_rationale`
- `predicted_defect_grammar`
- `predicted_rescue_or_corrector_grammar`
- `perturbation_predictions`
- `low_resolution_structure_or_ensemble_prediction`
- `proposed_experimental_tests`
- `falsification_criteria`
- `claim_boundary`

## Operator Buckets

Required buckets include:

- `membrane_domain_operator`
- `NBD1_stability_operator`
- `F508del_local_destabilization_operator`
- `interdomain_interface_coupling_operator`
- `trafficking_quality_control_operator`
- `corrector_or_rescue_context_operator`
- `multidomain_assembly_operator`
- `channel_function_context_operator`

Forbidden buckets:

- `generic_channel_annotation_only_operator`
- `single_local_mutation_only_operator`
- `compact_single_domain_fold_operator`
- `solved_atomic_structure_operator`
- `coordinate_contact_operator`
- `AlphaFold_confidence_proxy_operator`

## Pass Conditions

V46 passes only if:

- prediction is sealed before holdout validation
- no coordinate leakage occurs before sealing
- no internal runtime evidence is promoted as biology
- mechanism class is `membrane_multidomain_folding_assembly_proteostasis_defect`
- at least seven operator buckets are produced
- at least ten perturbation predictions are produced
- at least seven perturbation predictions are supported or partially supported by holdouts
- interdomain coupling support rate is at least `0.6`
- rescue logic support rate is at least `0.6`
- forbidden generic-channel-only grammar is rejected
- forbidden single-local-mutation-only grammar is rejected
- contradicted prediction count is at most two
- all leakage controls pass
- claim boundary remains honest

V46 may set:

- `live_membrane_solution_packet = true`
- `protein_folding_solved_candidate_strengthened = true`

only if pass conditions are met.

V46 must never set:

- `folding_problem_solved = true`

## Safe Claim Boundary

Allowed if passed:

`We have a sealed live solution packet for CFTR F508del as a membrane multidomain folding, assembly, trafficking, and rescue problem. This is not a universal protein-folding solved claim or an atomic-coordinate claim.`

Forbidden:

- we solved the universal protein-folding problem
- F508del-CFTR is explained by one local deletion only
- CFTR coordinates were predicted de novo
- generic channel annotation solves CFTR F508del
- external review is unnecessary
