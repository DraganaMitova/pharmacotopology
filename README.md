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

The particle-simulation side branch then suggested a narrow pressure-release
probe: high blocked-future pressure can be rescued only when the external
couplings show dense direct support, at least one top-rank direct constraint,
bounded score, decoy margin, cluster support, and compatible local shape.

```text
pressure_release_expanded_selected_event_count = 32
pressure_release_expanded_added_event_count = 2
pressure_release_expanded_added_native_long_range_contacts = 13
pressure_release_expanded_added_false_event_count = 0
pressure_release_expanded_false_nucleus_rate = 0.0
pressure_release_expanded_cluster_precision = 0.169434
pressure_release_expanded_long_range_recall = 0.287035
pressure_release_expanded_long_range_recall_delta_vs_edge_continuity = 0.017047
pressure_release_expanded_long_range_recall_margin_vs_matched_controls = 0.152411
pressure_release_expanded_long_range_recall_margin_vs_adversarial_controls = 0.147089
pressure_release_expanded_beats_matched_controls = true
pressure_release_expanded_beats_adversarial_calibrated_controls = true
pressure_release_expanded_claim_allowed = false
```

The next residual probe is a registry-extension selector. It adds trace
candidates only when they extend an already selected strong local register, or
when a low-future tail has positive coupling margin plus enough direct external
constraint density.

```text
registry_extension_expanded_selected_event_count = 35
registry_extension_expanded_added_event_count = 3
registry_extension_expanded_added_native_long_range_contacts = 13
registry_extension_expanded_added_false_event_count = 0
registry_extension_expanded_false_nucleus_rate = 0.0
registry_extension_expanded_cluster_precision = 0.160862
registry_extension_expanded_long_range_recall = 0.292161
registry_extension_expanded_long_range_recall_delta_vs_pressure_release = 0.005126
registry_extension_expanded_long_range_recall_margin_vs_matched_controls = 0.157537
registry_extension_expanded_long_range_recall_margin_vs_adversarial_controls = 0.152215
registry_extension_expanded_beats_matched_controls = true
registry_extension_expanded_beats_adversarial_calibrated_controls = true
registry_extension_expanded_claim_allowed = false
```

The terminal-bridge selector adds the final compatible long-range leftovers
that survive the current controls: one low-score bridge with direct external
support and one high-pressure tail with enough future-compatible support.

```text
terminal_bridge_expanded_selected_event_count = 37
terminal_bridge_expanded_added_event_count = 2
terminal_bridge_expanded_added_native_long_range_contacts = 5
terminal_bridge_expanded_added_false_event_count = 0
terminal_bridge_expanded_false_nucleus_rate = 0.0
terminal_bridge_expanded_cluster_precision = 0.155979
terminal_bridge_expanded_long_range_recall = 0.296334
terminal_bridge_expanded_long_range_recall_delta_vs_registry_extension = 0.004173
terminal_bridge_expanded_long_range_recall_margin_vs_matched_controls = 0.16171
terminal_bridge_expanded_long_range_recall_margin_vs_adversarial_controls = 0.156388
terminal_bridge_expanded_beats_matched_controls = true
terminal_bridge_expanded_beats_adversarial_calibrated_controls = true
terminal_bridge_expanded_claim_allowed = false
```

These are controlled recall gains, not universal improvements: cluster
precision drops from `0.179688` to `0.174316`, then to `0.172363`, because the
added boundary and edge events are small native-positive contacts. The
pressure-release rescue drops precision again to `0.169434` while raising
long-range recall to `0.287035`. Registry extension raises recall to
`0.292161`, but precision drops to `0.160862`. Terminal bridge raises recall
to `0.296334`, while precision drops to `0.155979`. This remains an 8-row
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

Run the default non-hanging focused suite. Use `python3 -m pytest` plus disabled
plugin autoload; bare `pytest` can inherit globally installed plugins that keep
processes alive after tests finish.

```bash
PYTHONPATH=src PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest -q
```

Visual/GIF tests are skipped by default. They run only when explicitly requested:

```bash
PYTHONPATH=src PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest -q --run-visual
```

Every test gets a hard timeout; override only when deliberately doing a long run:

```bash
PYTHONPATH=src PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
PHARMACOTOPOLOGY_TEST_TIMEOUT_SECONDS=120 \
python3 -m pytest -q --run-full-suite
```

Run the fast 4AKE AlphaFold source probe. Historical probe variants are off by
default; use `--suite-mode full` only when intentionally running the old broader
probe set. Each subprocess has its own timeout.

```bash
PYTHONPATH=src python3 scripts/run_4ake_alphafold_real_source_probe_v0.py \
  --suite-mode fast \
  --timeout-seconds 30
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

## 4AKE no-AlphaFold sequence-physics self-consistency V0

The no-AlphaFold loop now adds sequence-only physical priors to the contact
ensemble: pair contact energy, lightweight secondary-structure compatibility,
and current-graph degree consistency. These are independent of AlphaFold and
PDB coordinates, but they are still weak priors rather than a solved-folding
oracle.

Current checked regression result:

```text
sequence_physical_prior_kind = sequence_only_energy_secondary_structure_degree_prior_v0
final_pair_count = 195
final_long_range_pair_count = 95
final_contact_precision = 0.179487
final_contact_recall = 0.050505
final_precision_delta_vs_seed = +0.001822
final_recall_delta_vs_seed = -0.006734
long_range_recall = 0.098446
long_range_recall_delta_vs_seed = 0.000000
mean_final_contact_energy_score = 0.528571
mean_final_secondary_structure_score = 0.440308
mean_final_degree_consistency_score = 0.886718
mean_final_physical_prior_score = 0.591629
folding_problem_solved = false
external_independent_claim_allowed = false
```

That is an honest small precision lift, not a 4AKE crack. The loop now has the
missing physics vote, but the no-AlphaFold data channel remains too weak for a
solved-folding claim.

## DONE-mode 4AKE AlphaFold contract V0

Run the honest contact predictor contract without visual rendering:

```bash
PYTHONPATH=src python3 scripts/run_done_mode_contact_prediction_v0.py \
  --source-accession 4AKE:A \
  --predicted-pdb data/independent_contact_sources/AF-P69441-F1-model_v4.pdb \
  --predicted-pdb-chain A \
  --predicted-source-id alphafold_db_AF-P69441-F1-model_v4
```

GIF generation is off by default. A GIF is produced only with the explicit
`--render-visual-proof` opt-in:

```bash
PYTHONPATH=src python3 scripts/run_done_mode_contact_prediction_v0.py \
  --source-accession 4AKE:A \
  --predicted-pdb data/independent_contact_sources/AF-P69441-F1-model_v4.pdb \
  --predicted-pdb-chain A \
  --predicted-source-id alphafold_db_AF-P69441-F1-model_v4 \
  --render-visual-proof
```

A benchmark claim is allowed only when independent evidence exists and leakage
gates remain closed. Without independent evidence the runner abstains instead of
guessing.

Current locked 4AKE AlphaFold-source result:

```text
benchmark_claim_allowed = true
claim_rejection_reason = none
long_range_precision = 0.744681
long_range_recall = 0.906736
contact_precision = 0.774834
contact_recall = 0.418605
```

See `DONE_MODE_VISUAL_PROOF_V0.md` for the exact boundary.

## MSA-free single-sequence structure-source probe

4AKE without AlphaFold needs an independent source beyond DCA and hand-coded sequence physics. The repo now includes a bounded adapter for non-AlphaFold single-sequence predictors such as ESMFold/OmegaFold-style local PDB outputs:

```bash
PYTHONPATH=src python3 scripts/run_msa_free_single_sequence_structure_probe_v0.py \
  --source-accession 4AKE:A \
  --predicted-pdb /path/to/single_sequence_prediction.pdb \
  --predicted-source-id omegafold_single_sequence_4ake \
  --predicted-pdb-chain A
```

No predictor is run unless `--prediction-command` is explicitly provided, and that subprocess is bounded by `--timeout-seconds`. AlphaFold-like source IDs are rejected by default in this MSA-free probe.

