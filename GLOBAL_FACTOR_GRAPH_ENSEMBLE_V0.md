# Global Factor-Graph + Ensemble V0

## Question

The five-axis physics challenge showed that adding energy/free-energy,
entropy, local cooperativity, dynamics, and context is still not enough when
contacts are selected locally/pairwise.

This batch adds the next real missing layer: a global contact factor graph with
all-or-none cooperative contact groups, residue-degree mutex constraints, and a
deterministic top-K ensemble of global solutions.

## What landed

- `src/pharmacotopology/folding_global_factor_graph_ensemble.py`
  - sequence-only five-axis factor graph mode
  - learned-geometry factor graph mode for external non-AlphaFold predicted PDBs
  - all-or-none cooperative factors
  - residue-degree mutex constraints
  - top-K ensemble global solver
  - consensus contact probability selection
  - strict claim gates that keep universal physics locked unless sequence-only metrics justify it
- `scripts/run_global_factor_graph_ensemble_v0.py`
  - bounded runner
  - no internet calls
  - no predictor subprocesses
  - no GIF generation
  - no full-suite run
- `tests/test_global_factor_graph_ensemble_v0.py`
  - sequence-only factor graph leakage checks
  - sequence-only locked-subset rejection check
  - real 4AKE ESMFold learned-geometry factor graph success check

## Result A: sequence-only global factor graph on locked 8-row benchmark

The factor graph does **not** magically crack folding when it only receives the
five-axis sequence scores.

```text
source_mode = sequence_only_five_axis_factor_graph
row_count = 8
mean_native_contact_precision_after_audit = 0.107253
mean_native_contact_recall_after_audit = 0.108590
mean_long_range_contact_recall_after_audit = 0.041641
mean_contact_map_f1_after_audit = 0.106055
factor_graph_ensemble_claim_allowed = false
universal_physical_law_claim_allowed = false
folding_problem_solved = false
claim_rejection_reason = sequence_only_global_factor_graph_claim_rejected_for_rows:1MBN:A,2LZM:A,1TEN:A,1CSP:A,1TIM:A,1PGA:A,4AKE:A,1CLL:A
```

Interpretation: global grouping alone is not enough if the input field is still
sequence-only and weak. This rejects the false claim that forcing cooperative
patches can recover long-range structure from local physics alone.

## Result B: 4AKE with real non-AlphaFold ESMFold global geometry prior

When the factor graph receives the already-included real ESMFold PDB as an
external MSA-free learned global geometry prior, 4AKE passes strongly.

```text
source_mode = learned_geometry_prior_factor_graph
row_count = 1
mean_native_contact_precision_after_audit = 0.880068
mean_native_contact_recall_after_audit = 0.932021
mean_long_range_contact_recall_after_audit = 0.901554
mean_contact_map_f1_after_audit = 0.905300
factor_graph_ensemble_claim_allowed = true
universal_physical_law_claim_allowed = false
folding_problem_solved = true
claim_rejection_reason = global_factor_graph_ensemble_target_claim_survived_learned_geometry_gate_not_universal_physics
```

Interpretation: this is a valid target solve over a non-AlphaFold, no-MSA learned
single-sequence geometry prior. It is **not** a universal physical-law solve.

## Leakage boundary

For both modes:

```text
coordinate_truth_used_before_selection = false
native_truth_used_before_selection = false
raw_sequence_exposed = false
```

For the learned-geometry 4AKE mode:

```text
alphafold_used_before_selection = false
msa_used_before_selection = false
template_used_before_selection = false
learned_geometry_prior_used_before_selection = true
```

## Validation commands actually run

No full suite was run. CPU cleanup was run before focused test commands.

```bash
PYTHONPATH=src python scripts/kill_project_cpu_workers_v0.py --min-cpu 1.0
PYTHONPATH=src timeout 90s python -m pytest -q tests/test_global_factor_graph_ensemble_v0.py
PYTHONPATH=src timeout 120s python scripts/run_global_factor_graph_ensemble_v0.py --ensemble-size 6
PYTHONPATH=src timeout 60s python -m compileall -q src/pharmacotopology/folding_global_factor_graph_ensemble.py scripts/run_global_factor_graph_ensemble_v0.py
PYTHONPATH=src python scripts/kill_project_cpu_workers_v0.py --min-cpu 1.0
```

Validation output:

```text
tests/test_global_factor_graph_ensemble_v0.py: 3 passed in 7.95s
compileall: OK
final CPU cleanup: matched_count = 0
```

## Final honest decision

- Factor graph + ensemble is now implemented.
- It is the right next structural layer after five-axis local physics.
- It does **not** crack sequence-only universal folding.
- It does cleanly solve 4AKE when given the non-AlphaFold ESMFold global geometry prior.
- The next true unknown is not another local physics term; it is how to generate
a strong global geometry field without relying on learned priors or MSAs.
