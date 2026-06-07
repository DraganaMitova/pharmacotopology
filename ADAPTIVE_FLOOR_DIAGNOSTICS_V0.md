# ADAPTIVE_FLOOR_DIAGNOSTICS_V0

## What changed

This patch replaces the fixed inter-lobe future-preservation floor with a row-level adaptive profile.

The adaptive profile is computed from non-oracle row features:

- phase mode from top external coupling pair geometry (`diagonal`, `square`, `point`, etc.)
- sequence complexity as normalized Shannon entropy
- external coupling depth via `effective_sequence_count_over_length`
- target coverage
- row future-preservation ceiling as a noise guard

The old safety default remains `0.42`, but the effective floor can move down to `0.35` only when the row has enough phase/complexity/MSA support. A separate low-signal rescue is allowed only when the row profile is strong and the event still has enough physical/coupling support.

## Important result: 1CLL is not a perfect contact map

The old `7/7` result means seven frontier segment events pass the adaptive gate. It does **not** mean exact contact-map precision is perfect.

Generated diagnostic:

- file: `first_contact_clean_pharmacotopology_layer_run/1cll_contact_map_precision_v0.json`
- frontier events: `7`
- predicted region pairs: `432`
- true-positive region pairs: `21`
- false-positive region pairs: `411`
- region contact precision: `0.048611`
- region contact recall: `0.06422`
- perfect contact map: `false`

So: 1CLL has a useful frontier hit, not a perfect contact map.

## 4AKE / 1MBN / 1TIM compromise

Generated diagnostic:

- file: `first_contact_clean_pharmacotopology_layer_run/adaptive_floor_diagnostics_v0.json`
- table: `first_contact_clean_pharmacotopology_layer_run/adaptive_floor_diagnostics_rows_v0.csv`

Key row-level results from the locked 8-protein benchmark:

| Protein | phase | MSA depth/L | future ceiling | adaptive floor | adaptive pass | low-signal rescue | note |
|---|---:|---:|---:|---:|---:|---:|---|
| 1MBN:A | square | 13.078431 | 0.357241 | 0.355 | 10 | 2 | now attacked by adaptive floor/rescue |
| 4AKE:A | diagonal | 9.350467 | 0.256588 | 0.35 | 1 | 1 | weak but now has a guarded rescue path |
| 1TIM:A | diagonal | 8.101215 | 0.196174 | 0.42 | 0 | 0 | protected by low future-ceiling noise guard |
| 1CLL:A | diagonal | 13.52027 | 0.566687 | 0.35 | 97 | 1 | gap threshold still dominates at 0.481924 |

## Dataset expansion boundary

The archive contains 11 complete locked external-coupling rows across the all-locked manifest:

- original 8 benchmark rows
- 1UBQ holdout
- blind 10DC
- blind 10AF

Only two rows are labeled `multidomain_or_segmented` in the locked data currently available:

- `4AKE:A`
- `1CLL:A`

So this patch scans all complete locked rows in the archive, but it does not fabricate new multi-domain data and does not fetch new proteins from the internet.

## Verification run

Passed targeted checks:

```bash
PYTHONPATH=src python3 -m compileall -q src/pharmacotopology/folding_evolutionary_constraints.py src/pharmacotopology/folding_external_coupling_importer.py src/pharmacotopology/folding_coupling_nucleus_selector.py
PYTHONPATH=src python3 -m pytest -q tests/test_adaptive_coupling_floor_v0.py tests/test_perfect_7_7_1cll.py tests/test_1cll_a_gate_relaxation_verification.py
PYTHONPATH=src python3 -m pytest -q tests/test_real_external_sequence_to_dca_build_v0.py tests/test_blind_external_holdout_quality_switch_v0.py
PYTHONPATH=src python3 -m pytest -q tests/test_external_evolutionary_coupling_trace_loop_v0.py::test_external_importer_preregisters_rows_and_quality_gates tests/test_external_evolutionary_coupling_trace_loop_v0.py::test_coupling_dataset_reuses_constraint_row_grouping tests/test_external_evolutionary_coupling_trace_loop_v0.py::test_external_importer_rejects_disallowed_source_kind tests/test_external_evolutionary_coupling_trace_loop_v0.py::test_external_importer_requires_provenance_fields tests/test_external_evolutionary_coupling_trace_loop_v0.py::test_external_importer_rejects_duplicate_pairs tests/test_external_evolutionary_coupling_trace_loop_v0.py::test_external_importer_marks_coordinate_taint_without_laundering
```

Full `pytest -q` was also attempted with `PYTHONPATH=src`, but it did not finish inside the execution timeout. That appears consistent with the pre-existing long-running external trace-loop test behavior in this repo.
