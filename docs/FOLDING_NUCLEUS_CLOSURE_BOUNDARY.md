# Folding Nucleus Closure Boundary

This layer pivots from residue-pair thresholding to cooperative segment
closure events.

It does not claim a folding law has been discovered.

## Why This Layer Exists

The contact-law threshold search rejected the current scalar pair score:

```text
current_scalar_score_law_rejected = true
law_generalizes = false
```

The strongest failure was conceptual. Individual residue contacts are not
independent decisions. Long-range native contacts become plausible when
neighborhoods of residues can close cooperatively.

So this layer asks a different question:

```text
Can sequence-only segment closure events recover long-range native contact
regions better than pair-level thresholds?
```

## Event Definition

Each candidate closure event connects two sequence segments:

```text
segment A
segment B
candidate contact region between A and B
```

The event is scored from sequence-only features:

```text
contact_cluster_gain
secondary_structure_compatibility
hydrophobic_burial_gain
registry_support
- loop_entropy_cost
- geometry_violation_cost
- frustration_cost
- isolation_penalty
```

Coordinate-native labels are attached only after the closure event exists.

## Current Result

The nucleus closure search found a real signal, but not a law:

```text
candidate_closure_event_count = 5041
selected_threshold = 0.3
native_nucleus_recall = 0.50826
long_range_contact_recall_after_nucleus = 0.913242
pair_level_mean_long_range_contact_recall = 0.0
long_range_recall_delta_vs_pair_level = 0.913242
false_nucleus_rate = 0.692721
contact_cluster_precision = 0.032261
nucleus_level_long_range_beats_pair_level = true
nucleus_law_survives = false
```

Interpretation:

```text
Cooperative closure is a better object than isolated pair scoring for
recovering long-range native regions.

But the false closure rate is far too high, so this is not a discovered
folding nucleus law.
```

## Required Boundary Flags

```text
nucleus_level_scoring_completed = true
native_truth_used_before_event_generation = false
native_label_attached_after_event_generation = true
row_specific_nucleus_thresholds_forbidden = true
mechanism_discovery_claim_allowed = false
global_folding_claim_allowed = false
folding_problem_solved = false
```

## Allowed Claim

```text
Segment-level closure events recover long-range native contact regions better
than the tested pair-level threshold baseline, but they overgenerate false
nuclei.
```

## Forbidden Claim

```text
The project discovered the folding nucleus law.
The project predicts folding pathways.
The project solved protein folding.
```

## Next Target

The next useful work is not more contact-pair scoring. It is reducing false
closure events without destroying long-range recall:

```text
steric feasibility
motif orientation
domain emergence order
solvent burial consistency
non-native trap rejection
```
