# Visual Contact Repair Boundary

This document defines the boundary for contact-topology repair and native-gap
analysis.

The repair layer is allowed to ask:

```text
Which native contacts were missed after prediction?
Which predicted contacts were false after prediction?
Which missed contacts form clusters?
Which false contacts form clusters?
Did the model fail by local-only collapse, beta-pairing weakness,
premature compaction, membrane mis-topology, disorder over-collapse,
or architecture fragmentation?
```

It is not allowed to claim that protein folding is solved.

## Prediction Boundary

The baseline predictor may use only:

```text
raw sequence inside the locked local benchmark file
sequence length
sequence-only residue classes
sequence-only contact heuristics
```

The repair predictor may use only:

```text
baseline sequence-only contact candidates
sequence-only compact-anchor pressure
sequence-only beta-registry pressure
sequence-only disorder-like and membrane-like blockers
```

The repair predictor may not use:

```text
native contact pairs
native contact-map hash
truth axes
source label
PDB identity as an answer key
clinical data
drug, molecule, or protein design objectives
```

## Native-Gap Analysis Boundary

Native contacts enter only after baseline and repaired predictions exist.

Allowed post-prediction analysis:

```text
missed native contact clusters
false predicted contact clusters
long-range vs short-range recall split
beta-pairing contact recall
closure timing analysis
premature compaction flag
failure mechanism cohort
```

Forbidden post-prediction interpretation:

```text
global fold solved
accurate atomistic trajectory
clinical relevance
drug design
molecule generation
protein sequence design
```

## Required Metrics

The repair report must include:

```text
long_range_contact_recall
beta_pairing_contact_recall
premature_compaction_count
false_contact_cluster_count
native_cluster_miss_count
visual_failure_cohort_count
visible_partial_success_count
visible_failure_count
```

The safety certificate must include:

```text
native_truth_used_before_prediction = false
native_truth_used_before_repair = false
raw_sequence_exposed = false
global_folding_claim_allowed = false
folding_problem_solved = false
```

## Current Claim

The allowed claim is:

```text
The repository now explains visible contact failures and applies sequence-only
contact-topology repairs that improve long-range and beta-pairing recovery on
the locked visual benchmark.
```

The locked claim remains:

```text
This is not a solved protein folding engine.
```
