# Nucleus Graph Selectivity Boundary

This layer tests whether selected closure events behave like a small
graph-level folding nucleus rather than a ranked list of plausible closures.

It does not solve protein folding, discover a folding law, generate proteins,
design molecules, infer clinical meaning, or unlock global fold-class claims.

## Model Object

The graph view is:

```text
segments = nodes
closure events = edges
candidate nucleus = selected closure graph core
```

Selection may use sequence-derived closure features and graph proxies:

```text
competition score
nucleus score
cluster precision proxy
false-contact risk proxy
mutual support count
overlap abuse count
topology conflict count
trap graph pressure
```

Selection must not use:

```text
native contact labels
PDB coordinate contacts
source/fold labels as truth
raw sequences in exported artifacts
```

Native labels are attached after graph selection for scoring and falsification.

## Decoy Falsification

Each selected graph event is compared with a matched decoy from the same
protein. Decoys are matched by sequence-only properties such as:

```text
normalized span
hydrophobic burial gain
contact cluster gain
registry support
sequence span
```

A graph law would need selected nuclei to beat these matched decoys without
using native truth during selection.

## Current Result

The graph layer reduces the selected surface:

```text
pre_graph_selected_event_count = 686
post_graph_selected_event_count = 320
post_graph_selected_event_target_met = true
```

It keeps enough long-range recovery to pass the weak long-range target:

```text
pre_long_range_contact_recall = 0.588907
post_long_range_contact_recall = 0.372496
post_long_range_contact_recall_target_met = true
```

But it fails the important selectivity gates:

```text
pre_false_nucleus_rate = 0.594592
post_false_nucleus_rate = 0.609375
pre_contact_cluster_precision = 0.044608
post_contact_cluster_precision = 0.043311
real_vs_decoy_enrichment_ratio = 0.806452
decoy_enrichment_target_met = false
nucleus_graph_law_survives = false
```

## Interpretation

This is a useful falsification:

```text
closure events are the right object to inspect
competition reduces flooding
graph-core scoring reduces count further
matched decoys beat the selected graph
false nuclei remain too common
```

The current sequence-only graph selector is not mechanistic enough. The next
allowed work should look for a missing variable that separates real cooperative
closures from matched trap closures without using native labels before
selection.

