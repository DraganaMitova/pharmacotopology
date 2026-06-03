# Visual Folding Mechanism Boundary

This document defines the boundary for the visual folding mechanism workbench.

The workbench is allowed to show coarse mechanism evidence:

```text
sequence-only contact candidates
locked coarse native-contact targets
contact-map overlap
coarse closure trajectory
relative energy curve
visible failure cohorts
```

It is not allowed to claim that protein folding is solved.

## Purpose

Earlier folding batches made the benchmark safer by refusing collapsed class
claims and unsafe axis claims. That was the right safety move, but it also made
the project mostly tabular.

This layer adds a first visual mechanism target:

```text
Can a sequence-only contact-candidate model reconstruct any visible part of a
locked coarse native-contact map?
```

The answer is allowed to be partial. Failure is part of the artifact.

## What The Workbench Generates

The locked 12-row benchmark creates one visual packet per row:

```text
native_contact_map.svg
predicted_contact_map.svg
contact_map_overlay.svg
folding_trajectory.html
energy_curve.svg
contact_closure_curve.svg
coarse_grain_final.svg
```

The root artifacts are:

```text
visual_mechanism_12_report.json
visual_mechanism_12_rows.csv
visual_mechanism_12_contact_metrics.csv
visual_mechanism_12_failure_cohorts.csv
visual_mechanism_12_dashboard.html
visual_mechanism_12_certificate.json
```

## Prediction Boundary

Prediction may use:

```text
raw sequence from the locked data file
sequence length
sequence-only residue classes
sequence-only local and long-range contact heuristics
```

Prediction may not use:

```text
native contact pairs
native contact-map hash
truth axes
source label
PDB identity as an answer key
```

The report must keep:

```text
native_truth_used_before_prediction = false
raw_sequence_exposed = false
```

## Scoring Boundary

Native contacts enter only after prediction.

Allowed scoring:

```text
native_contact_recall
native_contact_precision
contact_map_f1
long_range_contact_recall
short_range_contact_recall
false_contact_rate
```

Allowed visual interpretation:

```text
visible_partial_success
low_recall_partial_failure
beta_long_range_pairing_failure
disorder_control_overclosure_failure
membrane_contact_topology_failure
architecture_boundary_contact_failure
```

Forbidden interpretation:

```text
global fold solved
accurate atomistic trajectory
clinical relevance
drug design
molecule generation
protein sequence design
```

## Dashboard Required Language

The dashboard must explicitly show:

```text
This Is A Mechanism Visualization, Not A Solved Folding Engine
Contact Maps Are The First Proof Target
Trajectory Is Coarse-Grained
Native Contacts Are Used Only After Prediction
Failures Are Visual Evidence
Global Folding Claim Remains Locked
```

## Success Criteria

This batch succeeds only if:

```text
visual artifacts generated for 12/12 rows
contact_map_f1 computed for 12/12 rows
at least 3 visible partial successes
failures visualized
native_truth_used_before_prediction = false
raw_sequence_exposed = false
global_folding_claim_allowed = false
folding_problem_solved = false
```

## Current Claim

The allowed claim is narrow:

```text
The repository now has a visual contact-map workbench for inspecting coarse
folding mechanism hypotheses and their failures.
```

The locked claim remains:

```text
This is not a protein folding engine.
```
