# V44 FUS-LC Live Unsolved Mechanism Attack

Status: `V44_FUS_LC_LIVE_UNSOLVED_SOLUTION_PACKET_PASSED_REVIEW_REQUIRED`
Target and region: `human FUS low-complexity domain / FUS P35637 residues 1-214`
live_unsolved_target_solution_packet: `True`
protein_folding_solved_candidate_strengthened: `True`
folding_problem_solved: `False`
Mechanism class: `intrinsic_disorder_phase_separation_contextual_ensemble`

## Operator Buckets
- `low_complexity_disorder_operator`
- `aromatic_sticker_pattern_operator`
- `polar_spacer_context_operator`
- `phosphorylation_shift_operator`
- `LLPS_self_association_operator`
- `fibril_or_gel_state_shift_operator`
- `context_bound_RNA_or_protein_interaction_operator`

## Perturbation Validation
- `V44_PERT_001` `partially_supported`: phosphomimetic or phosphorylation-like increase at N-terminal serine cluster
- `V44_PERT_002` `supported`: tyrosine/aromatic sticker disruption across FUS-LC
- `V44_PERT_003` `supported`: tyrosine pattern redistribution without changing total length
- `V44_PERT_004` `supported`: replace or delete the Q/G/S/Y low-complexity identity
- `V44_PERT_005` `supported`: force compact single-fold grammar
- `V44_PERT_006` `supported`: add RNA or phase-modulating protein partner
- `V44_PERT_007` `supported`: isolate or mutate fibril-core-prone windows
- `V44_PERT_008` `supported`: change pH, salt, concentration, temperature, or phosphorylation state
- `V44_PERT_009` `supported`: remove the N-terminal core-prone sticker region
- `V44_PERT_010` `partially_supported`: remove phosphorylation capacity with serine-to-alanine substitutions

## Contradicted Predictions
- none

## Proposed Experiments
- compare wild-type FUS-LC against S-to-E phosphomimetic and S-to-A phospho-null panels in turbidity, microscopy, FRAP, and fibril assays
- mutate or scramble tyrosine stickers while preserving Q/G/S composition and measure LLPS boundary shifts
- test 1-95, 39-95, 110-150, and 155-190 constructs for gel/fibril maturation and fuzzy-flank behavior
- measure salt, pH, concentration, temperature, and phosphorylation-state phase diagrams
- add RNA and Kapbeta2-like protein context and quantify condensate suppression, recruitment, or material-state redirection
- use NMR/HDX/crosslinking/single-molecule assays to test ensemble shifts without demanding a single coordinate model

## Falsification Criteria
- FUS-LC residues 1-214 reproducibly adopt one compact single native fold under ordinary conditions
- aromatic sticker disruption has no effect on self-association or phase behavior
- phosphorylation/phosphomimetic changes have no reproducible effect on LLPS, aggregation, or fibril competence
- RNA/protein partners do not shift ensemble or condensate behavior when binding is observed
- fibril/gel evidence requires a whole-domain native fold rather than local state-dependent ordering
- independent holdouts contradict more than two explicit perturbation predictions

## Leakage Status
- coordinate-derived before prediction: `0`
- internal runtime before prediction: `0`
- holdout leakage detected: `False`
- native metrics before prediction: `False`
- coordinate truth before prediction: `False`

## Controls
Passed `15` / `15`.

## Plain English Interpretation
V44 produced a sealed, source-separated FUS-LC live solution packet for an intrinsically disordered phase-separating ensemble problem. It supports a mechanism-language claim for FUS-LC residues 1-214, not a universal solved-protein-folding or atomic-coordinate claim.
