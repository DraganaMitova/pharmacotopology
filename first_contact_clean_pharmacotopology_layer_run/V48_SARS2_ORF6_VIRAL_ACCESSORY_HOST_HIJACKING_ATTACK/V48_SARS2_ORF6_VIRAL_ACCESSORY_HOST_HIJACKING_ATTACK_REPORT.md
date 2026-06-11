# V48 SARS-CoV-2 ORF6 Viral Accessory Host-Hijacking Attack

Status: `V48_SARS2_ORF6_VIRAL_ACCESSORY_HOST_HIJACKING_PASSED_REVIEW_REQUIRED`
Target and region: `SARS-CoV-2 ORF6 / accessory protein 6 / ns6 / ORF6 P0DTC6 full length 1-61 with C-terminal host-interface focus 38-61`
live_viral_host_hijacking_solution_packet: `True`
protein_folding_solved_candidate_strengthened: `True`
folding_problem_solved: `False`
Mechanism class: `viral_accessory_short_region_host_hijacking_disorder_interface`
V47 committed: `True`
V47 commit hash: `9ab5d600cd7e5f5abeebeae6f74aed9623be2b5d`
V48 committed: `False`

## Operator Buckets
- `C_terminal_host_interaction_operator`
- `RAE1_NUP98_binding_context_operator`
- `nuclear_transport_disruption_operator`
- `interferon_antagonism_context_operator`
- `short_linear_motif_or_MoRF_operator`
- `disorder_to_interface_operator`
- `localization_context_operator`
- `no_globular_single_fold_operator`

## Perturbation Validation
- `V48_PERT_001` `supported`: delete or disrupt ORF6 C-terminal residues 38-61
- `V48_PERT_002` `supported`: map RAE1/NUP98 interaction evidence
- `V48_PERT_003` `supported`: use generic viral accessory annotation only
- `V48_PERT_004` `supported`: measure nuclear import/export disruption after ORF6 expression
- `V48_PERT_005` `supported`: measure IFN/STAT/IRF antagonism
- `V48_PERT_006` `partially_supported`: add short-motif or MoRF-like evidence
- `V48_PERT_007` `supported`: force compact single-fold grammar
- `V48_PERT_008` `supported`: swap ORF8, ORF3a, or NSP evidence into validation
- `V48_PERT_009` `supported`: remove RAE1/NUP98 evidence
- `V48_PERT_010` `supported`: remove C-terminal/motif evidence
- `V48_PERT_011` `supported`: mutate C-terminal M58 or nearby acidic motif residues
- `V48_PERT_012` `partially_supported`: alter ER/Golgi localization region 18-24

## Contradicted Predictions
- none

## Proposed Experiments
- delete ORF6 residues 38-61 and test RAE1/NUP98 binding, nuclear transport, and IFN readouts
- mutate M58 and nearby acidic C-terminal residues and quantify NUP98-RAE1 binding
- compare N-terminal, central, and C-terminal deletions for host-transport blockade
- separate ER/Golgi localization mutations from C-terminal RAE1/NUP98 binding mutations
- measure STAT1/IRF3 nuclear translocation and ISG expression after C-terminal perturbation
- run negative controls with generic viral accessory, ORF8, ORF3a, and membrane-channel grammars

## Falsification Criteria
- C-terminal disruption does not weaken RAE1/NUP98 binding or nuclear transport disruption
- RAE1/NUP98 evidence does not localize to the ORF6 C-terminal region
- generic viral accessory annotation alone predicts the full packet
- ORF6 is best explained as a compact stable globular fold independent of host-interface context
- IFN antagonism proves an atomic fold rather than a functional consequence
- ORF8/ORF3a/NSP evidence validates ORF6-specific C-terminal predictions
- coordinate evidence before sealing is required for the packet
- independent holdouts contradict more than two explicit perturbation predictions

## Leakage Status
- coordinate-derived before prediction: `0`
- internal runtime before prediction: `0`
- holdout leakage detected: `False`
- native metrics before prediction: `False`
- coordinate truth before prediction: `False`

## Controls
Passed `18` / `18`.

## Plain English Interpretation
V48 produced a sealed, source-separated SARS-CoV-2 ORF6 packet in a fourth hard mechanism class: viral short-region host hijacking. The packet localizes the mechanism to a C-terminal RAE1/NUP98 interface, treats nuclear transport disruption and interferon antagonism as functional consequences, keeps localization as context, and rejects stable globular-fold and generic viral-accessory shortcuts without using coordinates before sealing.
