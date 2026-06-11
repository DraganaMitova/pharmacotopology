# V46 CFTR F508del Membrane Multidomain Folding Rescue Attack

Status: `V46_CFTR_F508DEL_MEMBRANE_MULTIDOMAIN_SOLUTION_PACKET_PASSED_REVIEW_REQUIRED`
Target and region: `human CFTR F508del membrane multidomain folding defect / full-length CFTR P13569 with F508del focus in NBD1`
live_membrane_solution_packet: `True`
protein_folding_solved_candidate_strengthened: `True`
folding_problem_solved: `False`
Mechanism class: `membrane_multidomain_folding_assembly_proteostasis_defect`

## Operator Buckets
- `membrane_domain_operator`
- `NBD1_stability_operator`
- `F508del_local_destabilization_operator`
- `interdomain_interface_coupling_operator`
- `trafficking_quality_control_operator`
- `corrector_or_rescue_context_operator`
- `multidomain_assembly_operator`
- `channel_function_context_operator`

## Perturbation Validation
- `V46_PERT_001` `supported`: delete F508 in NBD1
- `V46_PERT_002` `supported`: track F508del maturation and glycan processing
- `V46_PERT_003` `partially_supported`: stabilize NBD1 alone
- `V46_PERT_004` `supported`: restore NBD1-MSD/ICL interface coupling
- `V46_PERT_005` `supported`: apply corrector or rescue context
- `V46_PERT_006` `supported`: use generic channel annotation only
- `V46_PERT_007` `supported`: force soluble compact single-domain grammar
- `V46_PERT_008` `supported`: remove interdomain/interface evidence
- `V46_PERT_009` `supported`: remove NBD1 stability evidence
- `V46_PERT_010` `partially_supported`: remove trafficking/proteostasis evidence
- `V46_PERT_011` `partially_supported`: measure rescued channel opening at the membrane
- `V46_PERT_012` `supported`: combine NBD1, interface, and proteostasis correctors

## Contradicted Predictions
- none

## Proposed Experiments
- compare WT and F508del NBD1 stability/folding using non-coordinate biochemical stability assays
- measure full-length F508del glycan maturation and ER-to-surface trafficking with and without NBD1 stabilizers
- test interface/domain-coupling suppressors or correctors for synergy with NBD1 stabilization
- compare corrector classes for NBD1-only, interface, and proteostasis/maturation rescue signatures
- measure rescued surface channel opening separately from maturation rescue
- ablate proteostasis/quality-control components and quantify maturation versus degradation changes

## Falsification Criteria
- F508del is fully explained by one local residue deletion with no NBD1 stability effect
- F508del correction is complete with NBD1 stabilization alone in all relevant maturation/function assays
- interdomain/interface evidence does not affect correction or maturation readouts
- trafficking/proteostasis readouts do not change despite folding/maturation rescue
- generic channel annotation alone predicts the F508del rescue logic
- coordinate evidence before sealing is required for the mechanism packet
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
V46 produced a sealed, source-separated CFTR F508del packet outside the IDP family. The mechanism grammar is membrane multidomain folding: NBD1 destabilization, NBD1-MSD/ICL interface coupling, maturation/trafficking quality control, corrector/rescue logic, and downstream channel function context.
