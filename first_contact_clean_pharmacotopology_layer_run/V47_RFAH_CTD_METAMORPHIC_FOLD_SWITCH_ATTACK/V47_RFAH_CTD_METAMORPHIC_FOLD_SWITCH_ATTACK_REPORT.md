# V47 RfaH-CTD Metamorphic Fold-Switch Attack

Status: `V47_RFAH_CTD_METAMORPHIC_FOLD_SWITCH_PASSED_REVIEW_REQUIRED`
Target and region: `Escherichia coli RfaH C-terminal domain / RfaH P0AFW0 CTD residues 101-162`
live_fold_switch_solution_packet: `True`
protein_folding_solved_candidate_strengthened: `True`
folding_problem_solved: `False`
Mechanism class: `metamorphic_context_dependent_alpha_beta_fold_switch`
V46 committed: `True`
V46 commit hash: `6a9ae65be5ed300158d6a9cc119541a7260fc218`
V47 committed: `False`

## Operator Buckets
- `alpha_helical_hairpin_state_operator`
- `beta_barrel_or_beta_roll_state_operator`
- `NTD_bound_autoinhibition_operator`
- `CTD_release_switch_operator`
- `transcription_activation_context_operator`
- `translation_coupling_context_operator`
- `partner_context_refolding_operator`
- `no_single_consensus_fold_operator`

## Perturbation Validation
- `V47_PERT_001` `supported`: remove or weaken NTD-bound context
- `V47_PERT_002` `supported`: release CTD from NTD context
- `V47_PERT_003` `supported`: force one consensus CTD fold
- `V47_PERT_004` `supported`: remove beta-state evidence
- `V47_PERT_005` `supported`: remove alpha-state evidence
- `V47_PERT_006` `supported`: remove partner/context evidence
- `V47_PERT_007` `supported`: use generic transcription-factor annotation alone
- `V47_PERT_008` `supported`: treat RfaH-CTD as an IDP phase-separation target
- `V47_PERT_009` `supported`: treat RfaH-CTD as a membrane-pore target
- `V47_PERT_010` `partially_supported`: add transition or intermediate evidence
- `V47_PERT_011` `partially_supported`: disrupt released-CTD S10/ribosome partner context
- `V47_PERT_012` `supported`: stabilize the NTD-CTD interface

## Contradicted Predictions
- none

## Proposed Experiments
- compare CTD secondary-structure signatures in isolated, NTD-bound, and released/partner contexts using non-coordinate spectroscopy
- mutate or weaken NTD-CTD interface residues and measure alpha-state loss versus beta-state gain
- test released CTD binding to S10/ribosome-context partners and measure translation-coupling readouts
- force alpha- or beta-stabilizing substitutions and verify context-dependent functional tradeoffs
- measure kinetic intermediates after CTD release without using them as atomic-pathway proof
- run negative controls with generic transcription-factor and membrane-channel annotations only

## Falsification Criteria
- RfaH-CTD reproducibly has one context-independent consensus fold across NTD-bound and released contexts
- removing NTD-bound context does not change alpha-hairpin/autoinhibition grammar
- released CTD context does not support beta-state or translation-coupling grammar
- generic transcription-factor annotation alone predicts the full RfaH packet
- IDP phase-separation or membrane-channel grammar explains the packet better than fold switching
- coordinate evidence before sealing is required for the packet
- independent holdouts contradict more than two explicit perturbation predictions

## Leakage Status
- coordinate-derived before prediction: `0`
- internal runtime before prediction: `0`
- holdout leakage detected: `False`
- native metrics before prediction: `False`
- coordinate truth before prediction: `False`

## Controls
Passed `17` / `17`.

## Plain English Interpretation
V47 produced a sealed, source-separated RfaH-CTD packet in a third hard folding class: metamorphic alpha/beta fold switching. The packet separates NTD-bound alpha-hairpin/autoinhibition, CTD release, beta-state translation-coupling, partner-context refolding, and single-consensus-fold rejection without using coordinates before sealing.
