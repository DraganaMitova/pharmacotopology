# Contact Law Threshold Boundary

This layer tests whether a sequence-only contact law can survive a stricter
real-coordinate benchmark.

It does not claim that a folding law has been discovered.

## Question

The tested question is narrow:

```text
Can one threshold band over sequence-only contact features generalize across
held-out coordinate-backed proteins?
```

The benchmark starts from the real-coordinate visual 8-row surface. Native
contacts are derived from locked C-alpha coordinate traces, but only after
sequence-only feature generation.

## Feature Boundary

For every residue pair, the feature matrix computes:

```text
sequence_separation
normalized_separation
local_i_to_i4_support
helix_window_support
beta_window_support
hydrophobic_pair_support
aromatic_anchor_support
opposite_charge_support
same_charge_penalty
breaker_penalty
loop_entropy_cost
cluster_neighbor_support
parallel_contact_support
isolation_penalty
```

The native-contact label is attached only after those sequence-only features
exist:

```text
native_truth_used_before_feature_generation = false
native_label_attached_after_feature_generation = true
```

## Tested Scores

The threshold sweep tests:

```text
current_scalar_score
pair_only_score
pair_plus_cluster_score
pair_plus_entropy_score
pair_plus_cluster_plus_entropy_score
```

The leave-one-protein-out validation forbids row-specific threshold tuning:

```text
train on 7 proteins
choose threshold on train rows only
test on 1 held-out protein
```

## Current Result

The current scalar score is rejected as a law:

```text
current_scalar_score_best_global_threshold = 0.52
current_scalar_score_best_global_f1 = 0.150906
current_scalar_score_best_global_micro_f1 = 0.147772
current_scalar_score_threshold_stable = false
current_scalar_score_law_rejected = true
```

The best held-out candidate in this batch is:

```text
best_law_candidate_model = pair_plus_entropy_score
best_law_candidate_loo_mean_test_f1 = 0.233569
best_law_candidate_loo_threshold_std = 0.007071
best_law_candidate_survives = false
law_generalizes = false
```

It improves held-out F1 and reduces false contacts, but it fails the survival
rule because long-range recall collapses. That is useful falsification, not a
failure of the project.

## Required Claim Boundary

Required flags:

```text
law_search_completed = true
artifact_reproducible = true
current_scalar_score_law_rejected = true
native_truth_used_before_feature_generation = false
row_specific_thresholds_forbidden = true
mechanism_discovery_claim_allowed = false
global_folding_claim_allowed = false
folding_problem_solved = false
```

Allowed statement:

```text
The current scalar contact heuristic fails as a stable threshold law on the
real-coordinate benchmark. A pair-plus-entropy candidate improves held-out F1
but does not survive as a general contact law because it loses long-range
contact recovery.
```

Forbidden statement:

```text
The project discovered a folding law.
The project solved contact prediction.
The project can infer full protein folds.
```
