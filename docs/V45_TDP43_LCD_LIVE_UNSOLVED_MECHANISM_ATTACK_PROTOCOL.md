# V45 TDP-43 LCD Live Unsolved Mechanism Attack Protocol

## Purpose

V45 is the adversarial sibling replication after V44 FUS-LC. The target is human TDP-43 / TARDBP low-complexity C-terminal domain, residues 274-414.

The goal is not to replay the FUS-LC Q/G/S/Y tyrosine-ladder grammar. V45 must produce a sealed TDP-43-specific mechanism-language packet for sparse aromatic, glycine-rich, prion-like, helix-prone, nucleic-acid-context-sensitive ensemble behavior.

## Target

- Target: human TDP-43 low-complexity C-terminal domain
- Accession: UniProt `Q13148`
- DisProt: `DP01108`
- Region: residues `274-414`
- Expected mechanism class: `intrinsic_disorder_prion_like_phase_separation_contextual_ensemble`

## Allowed Prediction Inputs

- UniProt sequence and non-coordinate feature annotations
- DisProt disorder-only annotations
- low-complexity and compositional sequence features
- sequence-derived sparse aromatic, glycine-rich, Q/N, methionine, and serine patterning
- PTM annotations that are non-coordinate

## Blocked Prediction Inputs

- PDB/mmCIF coordinates before sealing
- coordinate-derived contacts or native contact maps
- AlphaFold/ESMFold/RoseTTAFold coordinates before sealing
- internal runtime reports as biological evidence
- validation holdouts before prediction sealing
- answer key/class labels
- target-name-only assignment
- FUS-LC solution transfer as biological evidence

## Sealing

The runner writes:

- `sealed_prediction_packet.json`
- `prediction_hash`
- `prediction_timestamp`
- `prediction_inputs_manifest.json`
- `blocked_inputs_manifest.json`
- `no_holdout_access_before_hash: true`

Only after the prediction packet hash exists does V45 write post-seal holdouts and validation scores.

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
- `glycine_rich_prion_like_lcd_operator`
- `sparse_aromatic_sticker_operator`
- `methionine_hydrophobic_cluster_operator`
- `alpha_helical_llps_segment_operator`
- `phosphorylation_pathology_shift_operator`
- `LLPS_self_association_operator`
- `amyloid_or_solid_state_maturation_operator`
- `nucleic_acid_context_switch_operator`
- `disease_mutation_modulation_operator`

Forbidden buckets:

- `compact_single_native_fold_operator`
- `solved_atomic_structure_operator`
- `de_novo_global_coordinate_solution_operator`
- `AlphaFold_confidence_proxy_operator`
- `coordinate_contact_operator`
- `fus_lc_tyrosine_ladder_transfer_operator`

## Pass Conditions

V45 passes only if:

- prediction is sealed before holdout validation
- no coordinate leakage occurs before sealing
- no internal runtime evidence is promoted as biology
- mechanism class is `intrinsic_disorder_prion_like_phase_separation_contextual_ensemble`
- at least eight operator buckets are produced
- at least nine perturbation predictions are produced
- at least six perturbation predictions are supported or partially supported by holdouts
- compact single-fold grammar is rejected
- FUS-LC transfer grammar is rejected
- contradicted prediction count is at most two
- all leakage controls pass
- claim boundary remains honest

V45 may set:

- `live_unsolved_target_solution_packet = true`
- `protein_folding_solved_candidate_strengthened = true`

only if pass conditions are met.

V45 must never set:

- `folding_problem_solved = true`

## Safe Claim Boundary

Allowed if passed:

`We have a sealed live solution packet for the TDP-43 LCD mechanism-language problem: a source-separated sparse-aromatic/prion-like ensemble and perturbation model for residues 274-414, supported by post-seal holdout evidence. This is not a universal protein-folding solved claim.`

Forbidden:

- we solved the universal protein-folding problem
- TDP-43 LCD has one solved atomic native structure
- FUS-LC grammar automatically solves all IDP LCDs
- coordinates were predicted de novo for all TDP-43 states
- external review is unnecessary
