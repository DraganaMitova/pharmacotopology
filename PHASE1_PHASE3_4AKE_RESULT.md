# Phase 1 / Phase 3 4AKE Diagnostic Result

Target checked: `4AKE:A`, row index `6`, row id `coord_007_pdb_4AKE_A_adenylate_kinase`.

Important correction: in this archive, `rows[7]` is `1CLL:A`, not `4AKE:A`.

## Phase 1 result

Raw candidate-pool region ceiling for 4AKE is `1.000`.

That means the broad candidate windows can cover all native long-range contacts. The immediate bottleneck is not absence of windows; it is the selector/frontier. However, the external DCA/coupling signal is weak: only `19/193` native long-range contacts are direct coupling hits after audit (`0.098446`).

## Frontier bottleneck

- Full candidate pool: `1225` events, ceiling `1.000`
- Competitive pool: `100` events, ceiling `0.217617`
- Seed trace-loop frontier: `12` events, ceiling `0.155440`
- Self-deciding expanded frontier: `14` events, ceiling `0.196891`
- Self-deciding generated frontier: `13` events, ceiling `0.165803`

## Phase 3 aggressive frontier probe

Implemented selector name:

```text
coupling_trace_loop_aggressive_relative_frontier
```

The selector uses the row-local seed median generation score as a baseline, then promotes broad candidate events above `0.5 * baseline` if the frontier grows by more than `1.2x`.

For 4AKE:

- Seed median baseline: `0.35414`
- Effective threshold: `0.17707`
- Aggressive frontier: `1078` events
- Raw region ceiling: `1.000`

## Collapse result

Aggressive frontier is rejected as a default solution:

- Collapsed long-range recall rises to `0.284974`
- Collapsed precision collapses to `0.011380`
- Collapsed long-range precision is `0.011687`

This is not a solved selector. It opens the map, but it does not select a clean folding nucleus.

## Decision

`rejected_precision_collapse`

Phase 3 proved the candidate windows exist, but also proved that the relative-threshold aggressive frontier is too broad. The honest next move is not another cheap frontier threshold. The next useful batch should either add an independent source of candidates/confidence, or switch to structure-level scoring/ensemble validation.

## Validation run

Focused validation passed:

```text
PYTHONPATH=src pytest -q tests/test_aggressive_relative_frontier_v0.py tests/test_self_deciding_frontier_generation_v0.py tests/test_contact_collapse_pipeline_integration_v0.py
7 passed in 16.24s
```

Full pytest was attempted but did not finish within the tool timeout in this environment.
