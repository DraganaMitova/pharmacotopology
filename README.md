# Protein Folding External Coupling Trace Loop

This orphan branch is a clean protein-folding workspace. It transfers only the
files needed for the current real-coordinate benchmark, the external
evolutionary coupling trace-loop selector, its provenance checks, matched
negative controls, selector sweep, and focused tests.

It intentionally does not include the earlier non-folding pharmacotopology
layer, profile dashboards, or legacy benchmark families. The preserved source
branch with the broader history is `protein-folding-hypothesis-lab`.

## Current Focus

```text
EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_V0
```

The active question is whether an external MSA/DCA-style coupling channel can
recover the native-selective constraint signal that the oracle positive control
showed would remove fake nuclei.

Current checked-in external run, strict rank-consistent selector:

```text
selected_event_count = 23
false_nucleus_rate = 0.0
cluster_precision = 0.162044
long_range_recall = 0.208635
vs_control_enrichment = 1.291737
vs_adversarial_calibrated_enrichment = 1.028713
hard_adversarial_calibrated_probe_passed = true
folding_problem_solved = false
claim_allowed = false
```

KNOT/particle-inspired persistence recovery now adds one calibrated trace event
that the hard cluster gate excluded:

```text
persistent_selected_event_count = 24
persistent_false_nucleus_rate = 0.0
persistent_cluster_precision = 0.167578
persistent_long_range_recall = 0.231865
persistent_vs_control_enrichment = 1.245257
persistent_vs_adversarial_calibrated_enrichment = 0.991698
persistent_selector_score_vs_control_enrichment = 1.262676
persistent_selector_score_vs_adversarial_calibrated_enrichment = 1.049653
persistent_recovered_event_count = 1
persistent_recovered_native_long_range_contact_count = 21
persistent_probe_passed = false
persistent_selector_score_probe_passed = true
```

This is a recovery diagnostic, not a solved-folding claim: persistence improves
false-nucleus control, precision, and recall. Raw coupling-only enrichment
remains just below the gates, while the selector's full non-oracle
coupling-nucleus score clears matched and adversarial enrichment controls.

Current recall frontier:

```text
persistent_recall_frontier_count = 23
score_margin_expansion_candidate_count = 2
score_margin_expansion_row_count = 2
score_margin_expansion_native_long_range_row_count = 2
score_margin_expansion_candidate_native_long_range_contacts = 38
score_margin_expansion_false_candidate_count = 0
score_margin_expansion_max_matched_control_candidate_count = 1
score_margin_expansion_max_matched_control_row_count = 1
score_margin_expansion_max_matched_control_native_long_range_contacts = 22
score_margin_expansion_margin_vs_matched_controls = +1 candidate, +1 row, +16 native long-range contacts
score_margin_expansion_max_adversarial_candidate_count = 1
score_margin_expansion_max_adversarial_row_count = 1
score_margin_expansion_max_adversarial_native_long_range_contacts = 22
score_margin_expansion_margin_vs_adversarial_controls = +1 candidate, +1 row, +16 native long-range contacts
score_margin_expansion_repeated_independent_row_signal_seen = true
score_margin_expansion_claim_allowed = false
```

Those two candidates showed where the next recall gain sat. The guarded
selector version now admits them using the same non-native score-margin gate:

```text
score_margin_expanded_selected_event_count = 26
score_margin_expanded_added_event_count = 2
score_margin_expanded_added_native_long_range_contacts = 38
score_margin_expanded_added_false_event_count = 0
score_margin_expanded_false_nucleus_rate = 0.0
score_margin_expanded_cluster_precision = 0.179688
score_margin_expanded_long_range_recall = 0.252423
score_margin_expanded_long_range_recall_delta_vs_persistent = 0.020558
score_margin_expanded_long_range_recall_margin_vs_matched_controls = 0.117799
score_margin_expanded_long_range_recall_margin_vs_adversarial_controls = 0.174647
score_margin_expanded_beats_matched_controls = true
score_margin_expanded_beats_adversarial_calibrated_controls = true
score_margin_expanded_claim_allowed = false
```

The KNOT Boundary Field Lab side project suggested one more bounded continuity
check: treat a low-cluster trace event as admissible only when it has score
margin, future preservation, low blocked-future pressure, and enough local
secondary-structure compatibility to avoid an overconfident weak-shape leak.
That selector adds a small recall gain while preserving the false-nucleus lock:

```text
boundary_continuity_expanded_selected_event_count = 28
boundary_continuity_expanded_added_event_count = 2
boundary_continuity_expanded_added_native_long_range_contacts = 2
boundary_continuity_expanded_added_false_event_count = 0
boundary_continuity_expanded_false_nucleus_rate = 0.0
boundary_continuity_expanded_cluster_precision = 0.174316
boundary_continuity_expanded_long_range_recall = 0.259572
boundary_continuity_expanded_long_range_recall_delta_vs_score_margin = 0.007149
boundary_continuity_expanded_long_range_recall_margin_vs_matched_controls = 0.124948
boundary_continuity_expanded_long_range_recall_margin_vs_adversarial_controls = 0.119626
boundary_continuity_expanded_beats_matched_controls = true
boundary_continuity_expanded_beats_adversarial_calibrated_controls = true
boundary_continuity_expanded_claim_allowed = false
```

KNOT-Core's transfer/falsification and stability-accumulation doctrine then
motivated a stricter edge-continuity probe: admit only modest-score,
high-cluster edge candidates with direct support, preserved future signal, low
blocked-future pressure, and enough local secondary-structure compatibility.

```text
edge_continuity_expanded_selected_event_count = 30
edge_continuity_expanded_added_event_count = 2
edge_continuity_expanded_added_native_long_range_contacts = 3
edge_continuity_expanded_added_false_event_count = 0
edge_continuity_expanded_false_nucleus_rate = 0.0
edge_continuity_expanded_cluster_precision = 0.172363
edge_continuity_expanded_long_range_recall = 0.269988
edge_continuity_expanded_long_range_recall_delta_vs_boundary_continuity = 0.010416
edge_continuity_expanded_long_range_recall_margin_vs_matched_controls = 0.135364
edge_continuity_expanded_long_range_recall_margin_vs_adversarial_controls = 0.130042
edge_continuity_expanded_beats_matched_controls = true
edge_continuity_expanded_beats_adversarial_calibrated_controls = true
edge_continuity_expanded_claim_allowed = false
```

These are controlled recall gains, not universal improvements: cluster
precision drops from `0.179688` to `0.174316`, then to `0.172363`, because the
added events are small native-positive contacts. This remains an 8-row
external-coupling benchmark result, not a solved-folding claim.

Claim mode remains locked. A folding-solved claim is refused unless the data and
per-constraint provenance stay external:

```text
external_evolutionary_couplings_used = true
coordinate_truth_used_to_build_constraints = false
native_truth_used_before_coupling_selection = false
oracle_constraint_control = false
```

## Useful Commands

Run the focused tests:

```bash
python3 -m pytest tests/test_external_evolutionary_coupling_trace_loop_v0.py tests/test_real_external_sequence_to_dca_build_v0.py
```

Regenerate the current external trace-loop artifacts:

```bash
python3 scripts/run_external_evolutionary_coupling_trace_loop_benchmark.py \
  --external-coupling-file data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json
```

Sweep selector variants without running the full historical suite:

```bash
python3 scripts/sweep_external_trace_loop_selector_v0.py \
  --max-configs 24 \
  --output /private/tmp/external_trace_loop_selector_sweep_v0.csv
```

The oracle coupling file is retained only as a positive control:

```text
data/folding_real_coordinate_visual_8_couplings.locked.json
```
