# V59 Real Folding-Process Replication Panel

V59 tests folding process, not only final structural compatibility.

It uses a panel of real proteins with public sequence plus process-level evidence from folding kinetics, phi-value analysis, mutation kinetics, intermediate/HDX/NMR-style literature, and fast-folder ensemble sources. V50-V58 grammar targets are excluded.

## Engine Revision

The first V59 run blocked and required engine hardening:

- generic `alpha_beta` metadata was too easily routed into metamorphic switching
- small helix / mini-protein folding-kinetics targets could abstain despite enough non-coordinate process metadata
- process route prediction needed to live in the engine rather than in a one-off runner

The revision is explicit and bounded:

- operator set unchanged
- mechanism-class set unchanged
- coordinate and AlphaFold leakage still blocked before sealing
- process holdouts still open only after sealed hashes exist

## Process Levels

V59 scores:

- P1 process class
- P2 early/late route order
- P3 mutation-sensitive region family
- P4 calibrated reliability

The passed run has:

```text
10 real process targets
10 accepted
3 accepted_with_caution
0 abstain_recommended
accepted_process_accuracy = 1.0
raw_process_accuracy = 1.0
controls = 14/14
engine_biology_modified = true
```

## Boundary

V59 does not set `folding_problem_solved = true`. It supports bounded process replication on this real process-evidence panel after the V59-required engine hardening.
