# Evolutionary Coupling Nucleus Boundary

This layer tests the missing native-selective constraint channel.

Earlier active physical selection showed that burial helps, but fake
hydrophobic burial remains too easy. This benchmark adds a locked coupling
constraint layer and asks whether closure selection improves when candidate
nuclei must:

```text
carry direct coupling support
preserve future coupling-supported closure paths
beat matched decoys by coupling score
build a compatible trace loop instead of only reranking events
```

It does not solve protein folding, discover a folding mechanism, generate
proteins, design molecules, infer clinical meaning, or unlock global fold-class
claims.

## What It Adds

New files:

```text
data/folding_real_coordinate_visual_8_couplings.locked.json
src/pharmacotopology/folding_evolutionary_constraints.py
src/pharmacotopology/folding_coupling_nucleus_selector.py
src/pharmacotopology/folding_coupling_decoy_falsification.py
scripts/run_evolutionary_coupling_nucleus_benchmark.py
tests/test_evolutionary_coupling_nucleus_benchmark.py
```

Generated artifacts:

```text
first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_report.json
first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_selectors.csv
first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_selected_events.csv
first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_assessments.csv
first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_decoys.csv
first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_dashboard.html
first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_certificate.json
```

## Critical Safety Boundary

The checked-in coupling file is an oracle-control layer:

```text
coupling_source_kind = coordinate_oracle_surrogate_for_missing_evolutionary_channel_v1
external_evolutionary_couplings_used = false
coordinate_truth_used_to_build_constraints = true
native_truth_used_before_coupling_selection = true
oracle_constraint_control = true
```

That is intentional. It answers:

```text
If a native-selective coupling channel existed, would the selector improve?
```

It does not answer:

```text
Can the current sequence-only model infer that channel?
```

Before any mechanism-discovery claim, replace the oracle-control file with
external MSA/DCA couplings and rerun the exact same benchmark.

## Current Result

The coupling rerank improves the broad search:

```text
physical_rerank_false_nucleus_rate = 0.559375
coupling_rerank_false_nucleus_rate = 0.34375

physical_rerank_cluster_precision = 0.050488
coupling_rerank_cluster_precision = 0.072461

physical_rerank_long_range_recall = 0.405008
coupling_rerank_long_range_recall = 0.565164

coupling_rerank_real_vs_decoy_enrichment_ratio = 1.931956
```

The iterative trace loop is the important new shape:

```text
coupling_trace_loop_selected_event_count = 54
coupling_trace_loop_false_nucleus_rate = 0.0
coupling_trace_loop_cluster_precision = 0.164128
coupling_trace_loop_long_range_recall = 0.397896
coupling_trace_loop_constraint_recall = 0.466536
coupling_trace_loop_real_vs_decoy_enrichment_ratio = 1.693887
coupling_selector_targets_met = true
```

Interpretation:

```text
native-selective constraints can kill fake burial
an iterative closure trace is better than another static score
the current result is oracle-bounded, not sequence-only discovery
the next real test is external evolutionary coupling input
```

## Non-Leakage Rule

Exports still must not include raw sequences. The coupling artifact uses only
row IDs, source accessions, residue-pair indexes, confidence values, and safety
metadata.

For the checked-in oracle-control benchmark, native truth is allowed only
because the benchmark explicitly labels that truth dependency and keeps all
mechanism claims disabled.

For a real external-coupling run, these fields must change:

```text
external_evolutionary_couplings_used = true
coordinate_truth_used_to_build_constraints = false
native_truth_used_before_coupling_selection = false
oracle_constraint_control = false
```

Only then can the same selector be evaluated as a non-oracle folding signal.
