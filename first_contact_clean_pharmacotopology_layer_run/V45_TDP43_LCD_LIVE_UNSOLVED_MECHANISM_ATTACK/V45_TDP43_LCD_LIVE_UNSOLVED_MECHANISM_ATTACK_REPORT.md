# V45 TDP-43 LCD Live Unsolved Mechanism Attack

Status: `V45_TDP43_LCD_LIVE_UNSOLVED_SOLUTION_PACKET_PASSED_REVIEW_REQUIRED`
Target and region: `human TDP-43 low-complexity C-terminal domain / TDP-43 Q13148 residues 274-414`
live_unsolved_target_solution_packet: `True`
protein_folding_solved_candidate_strengthened: `True`
folding_problem_solved: `False`
Mechanism class: `intrinsic_disorder_prion_like_phase_separation_contextual_ensemble`

## Operator Buckets
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

## Perturbation Validation
- `V45_PERT_001` `supported`: mutate sparse W/F/Y aromatic stickers
- `V45_PERT_002` `supported`: disrupt W334-centered and F-rich sticker neighborhoods
- `V45_PERT_003` `supported`: mutate or break the 321-343 helix-prone CTD segment
- `V45_PERT_004` `supported`: remove glycine-rich/prion-like low-complexity identity
- `V45_PERT_005` `supported`: introduce ALS-linked CTD mutations near 287-393
- `V45_PERT_006` `partially_supported`: increase phosphorylation or phosphomimetic charge in the CTD/LCD
- `V45_PERT_007` `supported`: add RNA/DNA or restore multidomain nucleic-acid-binding context
- `V45_PERT_008` `supported`: isolate C-terminal fragments or extend stress exposure
- `V45_PERT_009` `supported`: force compact single-native-fold grammar
- `V45_PERT_010` `supported`: transfer FUS-LC dense tyrosine-ladder grammar onto TDP-43 LCD
- `V45_PERT_011` `partially_supported`: alter methionine-rich mid-CTD and C-terminal tail clusters

## Contradicted Predictions
- none

## Proposed Experiments
- mutate W/F/Y stickers, especially W334 and nearby F-rich neighborhoods, and measure turbidity, microscopy, FRAP, and NMR shifts
- break the 321-343 helix-prone CTD segment and compare phase behavior to wild type
- compare ALS-linked CTD mutants across LLPS, aggregation, and nucleic-acid-binding contexts
- test phosphomimetic and phospho-null variants in CTD condensate and aggregate assays
- add RNA/DNA or restore multidomain context and quantify viscosity, elasticity, recruitment, and self-association shifts
- compare full CTD 274-414 to C-terminal fragments under aging/stress for amyloid or solid-like maturation

## Falsification Criteria
- TDP-43 LCD residues 274-414 reproducibly adopt one compact single native fold under ordinary conditions
- W/F/Y aromatic sticker mutations do not alter phase behavior or self-association
- 321-343 helix-prone segment perturbation has no effect on LLPS or ensemble signatures
- nucleic-acid or multidomain context does not shift TDP-43 condensate material properties when binding/context is present
- ALS-linked CTD mutations do not perturb LLPS or aggregation readouts in any reproducible direction
- FUS-LC dense tyrosine-ladder grammar validates TDP-43 LCD predictions without TDP-specific evidence
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
V45 produced a sealed, source-separated TDP-43 LCD live solution packet that is not just FUS-LC replay: it uses sparse aromatics, helix-prone CTD grammar, glycine-rich prion-like disorder, mutation modulation, nucleic-acid context, and amyloid/solid maturation as distinct operator language.
