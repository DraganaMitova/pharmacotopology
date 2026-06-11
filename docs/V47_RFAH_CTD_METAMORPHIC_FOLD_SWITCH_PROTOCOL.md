# V47 RfaH-CTD Metamorphic Fold-Switch Protocol

V47 tests whether the mechanism-language workflow can handle a true metamorphic fold-switching target.

The target is the Escherichia coli RfaH C-terminal domain, residues 101-162. The scientific question is not one static structure. It is whether a sealed packet can separate an NTD-bound alpha-helical hairpin/autoinhibited grammar from a released beta-barrel or beta-roll translation-coupling grammar.

## Scope

- Target: Escherichia coli RfaH C-terminal domain
- UniProt: P0AFW0 / RFAH_ECOLI
- Region: residues 101-162
- Mechanism class: `metamorphic_context_dependent_alpha_beta_fold_switch`
- Runtime claim: no MD, no coordinate prediction, no universal folding-solved claim

## Allowed Prediction Inputs

- RfaH sequence
- UniProt non-coordinate function and interaction annotations
- domain architecture annotations
- literature-derived state and partner-context labels
- non-coordinate fold-switch annotations
- derived sequence-composition features

## Blocked Prediction Inputs

- PDB/mmCIF coordinates before sealing
- coordinate-derived contacts, distances, or native contact maps
- AlphaFold/ESMFold/RoseTTAFold coordinates before sealing
- UniProt secondary-structure records when their evidence is PDB-derived
- internal runtime reports as biological evidence
- validation holdouts before sealing
- target-name-only or generic transcription-factor assignment
- swapped XCL1/KaiB/Mad2 evidence as RfaH validation

## Sealing Rule

V47 writes `sealed_prediction_packet.json` before any holdout validation.

The sealed packet includes:

- `prediction_hash`
- `prediction_timestamp`
- `prediction_inputs_manifest.json`
- `blocked_inputs_manifest.json`
- `no_holdout_access_before_hash: true`

Only after the hash exists does V47 write post-seal holdouts and validation scores.

## Required Operators

- `alpha_helical_hairpin_state_operator`
- `beta_barrel_or_beta_roll_state_operator`
- `NTD_bound_autoinhibition_operator`
- `CTD_release_switch_operator`
- `transcription_activation_context_operator`
- `translation_coupling_context_operator`
- `partner_context_refolding_operator`
- `no_single_consensus_fold_operator`

## Forbidden Operators

- `compact_single_native_fold_operator`
- `generic_two_state_annotation_only_operator`
- `intrinsic_disorder_phase_separation_operator`
- `membrane_channel_operator`
- `solved_atomic_structure_operator`
- `coordinate_contact_operator`
- `AlphaFold_confidence_proxy_operator`

## Required Perturbation Logic

V47 must emit at least 10 explicit perturbation predictions. The required classes are:

- weakening NTD-bound context destabilizes alpha/autoinhibited grammar
- CTD release favors beta/translation-coupling grammar
- forcing one consensus fold is blocked
- removing beta-state evidence weakens or invalidates the packet
- removing alpha-state evidence weakens or invalidates the packet
- removing partner/context evidence weakens the full switch packet
- generic transcription-factor annotation alone fails
- IDP phase-separation grammar fails
- membrane-pore grammar fails
- intermediate evidence can support switch-path grammar but not an atomic pathway-solved claim

## Pass Conditions

V47 passes only if:

- prediction is sealed before holdout validation
- coordinate and internal-runtime leakage counts are zero
- mechanism class is exactly `metamorphic_context_dependent_alpha_beta_fold_switch`
- at least 7 operator buckets are produced
- at least 10 perturbation predictions are produced
- at least 7 perturbations are supported or partially supported by post-seal holdouts
- `state_separation_support_rate >= 0.7`
- `partner_context_support_rate >= 0.6`
- single-consensus-fold grammar is rejected
- wrong-class grammars are rejected
- contradictions are at most 2
- all leakage and shortcut controls pass
- claim boundary remains honest

V47 may set:

```text
live_fold_switch_solution_packet = true
protein_folding_solved_candidate_strengthened = true
```

only when the pass conditions are met.

V47 must never set:

```text
folding_problem_solved = true
```

## Honest Claim Boundary

Allowed wording:

`We have a sealed live fold-switch solution packet for RfaH-CTD as a context-dependent alpha/beta metamorphic mechanism. This is not a universal protein-folding solved claim or an atomic transition-pathway claim.`

Forbidden wording:

- universal protein folding is solved
- RfaH-CTD has one solved consensus fold in all contexts
- RfaH-CTD transition coordinates were predicted de novo
- generic two-state annotation solves RfaH
- external review is unnecessary
