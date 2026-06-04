# Active Physical Selection Boundary

This layer makes the coarse physical closure-state evaluator active.

Earlier physical scoring was an audit layer over graph-selected closure events.
This layer asks whether the same physical terms can control selection:

```text
graph selector only
graph + active physical rerank
graph + physical viability gate
graph + physical viability gate + future-frustration gate
```

It does not solve protein folding, discover a folding mechanism, generate
proteins, design molecules, infer clinical meaning, or unlock global fold-class
claims.

## What It Adds

The active selector uses these coarse terms:

```text
burial_gain
loop_strain
steric_clash_score
unsatisfied_polar_penalty
future_frustration_score
matched-decoy physical margin
```

It also runs term ablations:

```text
without_loop_strain
without_steric_clash
without_burial_gain
without_unsatisfied_polar_penalty
without_future_frustration
without_decoy_margin
```

Native contacts are attached only after selection for scoring.

## Non-Leakage Rule

Active physical selection may use sequence-derived closure features and coarse
physical-state terms.

It must not use:

```text
native contact labels
PDB coordinate contacts
source/fold labels as truth
raw sequences in exported artifacts
```

## Current Result

The active selector comparison is mixed and rejected:

```text
graph_only_false_nucleus_rate = 0.609375
physical_rerank_false_nucleus_rate = 0.559375
physical_gate_false_nucleus_rate = 0.316338
future_frustration_false_nucleus_rate = 0.309455

graph_only_cluster_precision = 0.043311
physical_rerank_cluster_precision = 0.050488
physical_gate_cluster_precision = 0.069038
future_frustration_cluster_precision = 0.071567

graph_only_long_range_recall = 0.372496
physical_rerank_long_range_recall = 0.405008
physical_gate_long_range_recall = 0.088436
future_frustration_long_range_recall = 0.070742

physical_rerank_real_vs_decoy_enrichment_ratio = 1.585882
active_physical_selection_survives = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

The interpretation is:

```text
physical rerank improves decoy enrichment and long-range recall
physical gates reduce false nuclei and improve precision
physical gates destroy too much long-range recall
no active selector passes all survival gates
```

## Term Ablation

Current ablation says only one term is clearly useful under this proxy:

```text
best_physical_term = burial_gain
worst_physical_term = future_frustration
physical_terms_with_positive_ablation_effect = burial_gain
physical_terms_rejected_as_noise =
  loop_strain
  steric_clash
  unsatisfied_polar_penalty
  future_frustration
  decoy_margin
```

This does not prove that the rejected terms are biologically irrelevant. It
only says the current coarse proxy does not use them reliably enough.

## Interpretation

This is a useful negative result:

```text
physical terms are informative enough to rerank decoys
physical terms are not yet sufficient to select native-positive closure nuclei
the current future-frustration gate is too blunt
the current active selector does not survive
```

The next missing variable is probably not another simple threshold. It is more
likely orientation, side-chain packing, solvent exposure, explicit polar
satisfaction, or a better kinetic closure order.
