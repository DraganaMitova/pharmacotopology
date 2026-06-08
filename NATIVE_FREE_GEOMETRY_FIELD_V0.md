# Native-free Geometry Field V0

This batch tests the next honest target after the global factor-graph result:

> Can the project generate a global geometry field without AlphaFold/ESMFold,
> templates, or any trained structure predictor?

## What changed

Added a bounded native-free geometry-field layer:

- `src/pharmacotopology/folding_native_free_geometry_field.py`
- `scripts/run_native_free_geometry_field_v0.py`
- `tests/test_native_free_geometry_field_v0.py`

The layer has two claim-relevant modes:

1. `pure_sequence_symbolic_geometry`
   - sequence chemistry
   - lightweight symbolic secondary/topological priors
   - global all-or-none factor graph
   - no MSA/DCA, no learned geometry, no template

2. `external_dca_geometry_field`
   - safe external evolutionary coupling anchors
   - sequence symbolic topology
   - contact-map-grid anchor diffusion
   - global all-or-none factor graph
   - no AlphaFold/ESMFold/template/trained structure model
   - not a universal physical-law claim because DCA comes from MSA/evolutionary history

Native coordinates/native contacts are attached only after selection for audit.

## Focused result on locked 8-row benchmark

Runner:

```bash
PYTHONPATH=src python3 scripts/run_native_free_geometry_field_v0.py --mode both
```

No full suite, no GIF, no predictor subprocesses, no internet.

### Pure sequence symbolic geometry

```text
mean_native_contact_precision_after_audit = 0.145181
mean_native_contact_recall_after_audit = 0.143068
mean_long_range_contact_recall_after_audit = 0.038526
mean_contact_map_f1_after_audit = 0.143925
native_free_geometry_field_claim_allowed = false
universal_physical_law_claim_allowed = false
folding_problem_solved = false
```

### External DCA geometry field

```text
mean_native_contact_precision_after_audit = 0.176080
mean_native_contact_recall_after_audit = 0.171665
mean_long_range_contact_recall_after_audit = 0.215862
mean_contact_map_f1_after_audit = 0.173540
native_free_geometry_field_claim_allowed = false
universal_physical_law_claim_allowed = false
folding_problem_solved = false
```

## Interpretation

This is not the crack.

The external DCA geometry field improves the long-range signal versus the earlier
sequence-only factor-graph result, but it is still far below the strict folding
claim gate. The missing piece is not merely a factor graph and not merely a
hand-written geometry field. The evidence now points to a stronger requirement:

a dense, accurate global geometry field that is not recoverable from the current
sequence-only symbolic priors or the available weak external DCA anchors.

## Validation

```text
CPU cleanup before focused tests: matched_count = 0
focused geometry-field tests: 3 passed in 1.31s
bounded runner completed
compileall: OK
final CPU cleanup: matched_count = 0
```
