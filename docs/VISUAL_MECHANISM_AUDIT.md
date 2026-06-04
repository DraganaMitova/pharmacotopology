# Visual Mechanism Audit

This audit freezes the current visual-contact benchmark boundary.

## Current Truth

The repository has a useful visual mechanism workbench:

```text
native contact maps
predicted contact maps
native vs predicted overlays
coarse folding trajectories
energy and closure curves
visible failure cohorts
```

It also has a contact-repair pass that improves the tiny locked benchmark:

```text
visible partial successes: 5 -> 8
visible failures: 7 -> 4
mean contact-map F1: 0.104031 -> 0.263729
long-range contact recall improved
beta-pairing contact recall improved
```

This is useful, but it is not mechanism discovery.

## Frozen Claim Boundary

Allowed claim:

```text
The project can visualize and audit coarse contact-map hypotheses and failures.
```

Forbidden claim:

```text
The project discovered how proteins fold.
```

The visual 12-row benchmark is:

```text
toy
coarse
internal
contact-map-only
not a physical folding benchmark
not a global fold benchmark
not a clinical, drug, molecule, or protein-design benchmark
```

## Overfit Risk

The contact repair uses sequence-only heuristics, but some are hardcoded
template families:

```text
compact beta registry centers
fixed cross-sheet anchors
compact long-range anchors
```

Those templates can fit the locked 12-row contact targets too well. The repair
improvement must therefore be read as:

```text
A coarse repair heuristic improved a tiny locked contact-map toy benchmark.
```

Not:

```text
The model discovered protein-folding mechanism.
```

## Required Audit Flags

The generated audit must keep:

```text
visual_12_is_toy_benchmark = true
contact_repair_overfit_risk_reported = true
artifact_reproducible = true
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
global_folding_claim_allowed = false
raw_sequence_exposed = false
native_truth_used_before_prediction = false
native_truth_used_before_repair = false
```

## Archive Rule

Share the project with a tracked-file archive:

```bash
git archive --format=zip --output /private/tmp/pharmacotopology-clean.zip HEAD
```

Do not use a Finder zip as reproducibility evidence. Finder zips can include
local cache files or omit the exact tracked state.

## Next Allowed Work

Before adding more folding logic, the next work must preserve this audit:

```text
keep visual 12 marked toy/internal
report template overfit risk
keep native truth after prediction
keep raw sequences out of generated reports
keep folding_problem_solved = false
```
