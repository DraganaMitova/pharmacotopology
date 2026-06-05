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

Current checked-in external run:

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
