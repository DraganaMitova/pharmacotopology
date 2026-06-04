# Competitive Nucleus Selection Boundary

This layer audits the folding-nucleus closure search with a competitive selection pass.

It does not solve protein folding, discover a folding mechanism, infer clinical value,
design molecules, or generate protein sequences.

## What It Adds

The previous nucleus closure search produced many candidate closures:

```text
candidate_closure_event_count = 5041
accepted_event_count = 3846
false_nucleus_rate = 0.692721
contact_cluster_precision = 0.032261
long_range_contact_recall_after_nucleus = 0.913242
```

That result recovered long-range native regions, but it did so by accepting too many
closures. The competitive layer adds:

```text
closure compatibility graph
sequence-only frustration filter
geometry rejection counts
overlap rejection counts
trap rejection counts
competitive budgeted selection
native scoring only after selection
```

## Non-Leakage Rule

Selection may use sequence-derived closure features:

```text
nucleus_score
closure_event_stability
contact_cluster_gain
registry_support
loop_entropy_cost
geometry_violation_cost
frustration_cost
isolation_penalty
```

Selection must not use:

```text
native_contact_count_after_scoring
native_long_range_contact_count_after_scoring
PDB-derived native contacts
fold labels
source labels as truth
raw sequences in exported artifacts
```

Native coordinates are attached only after event generation and after competitive
selection so the layer can be audited.

## Current Result

The current competitive selector reduces the event flood and keeps useful long-range
coverage:

```text
pre_competition_event_count = 3846
post_competition_selected_event_count = 686
event_reduction_ratio = 0.821633
pre_long_range_contact_recall = 0.913242
post_long_range_contact_recall = 0.588907
```

But the law does not survive because false closures and contact-cluster precision are
still weak:

```text
pre_false_nucleus_rate = 0.692721
post_false_nucleus_rate = 0.594592
pre_contact_cluster_precision = 0.032261
post_contact_cluster_precision = 0.044608
nucleus_competition_law_survives = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

## Interpretation

This is a useful failure. It says:

```text
competitive selection can reduce closure flooding
competitive selection can preserve some long-range recovery
competitive selection is not yet precise enough to claim a folding mechanism
```

The next biology-facing work should improve false-contact rejection or cluster
precision without using native labels before selection.

