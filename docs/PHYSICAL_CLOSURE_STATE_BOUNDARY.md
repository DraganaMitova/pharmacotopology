# Physical Closure State Boundary

This layer instantiates graph-selected closure candidates as coarse physical
states.

It does not solve protein folding, simulate atomistic dynamics, generate
proteins, design molecules, infer clinical meaning, or unlock global fold-class
claims.

## What It Adds

Previous layers scored:

```text
pair contacts
closure events
competitive closure selections
closure graph cores
matched decoys
```

This layer adds a coarse physical-state packet for each graph-selected closure:

```text
loop_strain
steric_clash_score
burial_gain
unsatisfied_polar_penalty
future_frustration_score
physical_state_score
```

The state is still sequence-only and coarse. It is not an atomistic model, a
rotamer model, a solvent simulation, or a folding trajectory.

## Non-Leakage Rule

Physical scoring may use sequence-derived properties:

```text
segment residue classes
loop composition
hydrophobic/aromatic/polar/charged fractions
closure span
existing sequence-only closure features
```

Physical scoring must not use:

```text
native contact labels
PDB coordinate contacts
source/fold labels as truth
raw sequences in exported artifacts
```

Native contacts are attached only after physical scoring for audit metrics.

## Current Result

The physical evaluator successfully builds all candidate states:

```text
candidate_state_count = 320
state_build_success_count = 320
state_build_failure_count = 0
```

Physical scoring does beat matched decoys by score:

```text
mean_physical_state_score = 0.41558
mean_decoy_physical_state_score = 0.354731
real_vs_decoy_physical_enrichment_ratio = 1.171535
real_beats_decoy_physical_score_rate = 0.65625
physical_state_rank_enrichment_at_25 = 1.1264
physical_enrichment_target_met = true
```

But the native/contact gates still reject the law:

```text
post_physical_false_nucleus_rate = 0.609375
post_physical_contact_cluster_precision = 0.043311
post_physical_long_range_contact_recall = 0.372496
post_physical_false_nucleus_rate_target_met = false
post_physical_contact_cluster_precision_target_met = false
post_physical_long_range_contact_recall_target_met = true
physical_state_law_survives = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

## Interpretation

This is progress, but not a solution:

```text
physical score sees something decoys do not
physical score does not yet separate native-positive nuclei well enough
false nuclei remain too common
contact-cluster precision remains too low
```

The next missing variable is likely more physical than this coarse evaluator:
orientation, side-chain packing, solvent exposure, explicit polar satisfaction,
or kinetic order.

