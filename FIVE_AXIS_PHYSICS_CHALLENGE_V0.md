# Five-Axis Physics Challenge V0

This is the honest “not another patch” layer for the current protein-folding diagnosis.

It does **not** claim that folding is solved. It asks one falsifiable question:

> If we add the five missing physical axes — energy/free energy, entropy, cooperativity, dynamics, and environment/context — can a sequence-only physical rule recover native contacts on the locked real-coordinate benchmark without learned priors, MSA, templates, native labels, or coordinate leakage?

## What landed

Added:

- `src/pharmacotopology/folding_five_axis_physics.py`
- `scripts/run_five_axis_physics_challenge_v0.py`
- `scripts/kill_project_cpu_workers_v0.py`
- `tests/test_five_axis_physics_challenge_v0.py`

Hardened:

- `scripts/run_fast_suite_v0.py` now uses a subprocess timeout and exits `124` if pytest hangs.

## Five axes now represented

For every candidate contact, before native labels are attached, the layer computes:

1. `energy_support_score`
2. `entropy_retention_score`
3. `cooperative_neighbour_score`
4. `environmental_context_score`
5. `dynamic_ensemble_score`
6. `free_energy_support_score` from energy + entropy

The dynamic proxy has compact/open/hinge state scores. This is not molecular dynamics. It is a deterministic sequence-only ensemble proxy designed to expose whether the missing axes are enough to create a stronger contact law.

## Leakage boundary

The challenge packet keeps these locked false before selection:

```text
coordinate_truth_used_before_selection = false
native_truth_used_before_selection = false
learned_prior_used_before_selection = false
msa_used_before_selection = false
template_used_before_selection = false
raw_sequence_exposed = false
```

Native contacts are attached only after selection for audit.

## Locked benchmark result

Command:

```bash
python3 scripts/kill_project_cpu_workers_v0.py --min-cpu 5
PYTHONPATH=src timeout 120s python3 scripts/run_five_axis_physics_challenge_v0.py
```

Result on `data/folding_real_coordinate_visual_8.locked.json`:

```text
row_count = 8
mean_native_contact_precision_after_audit = 0.387362
mean_native_contact_recall_after_audit = 0.039198
mean_long_range_contact_recall_after_audit = 0.007464
mean_contact_map_f1_after_audit = 0.069241
mean_f1_margin_vs_best_control = 0.013417
mean_long_range_recall_margin_vs_best_control = 0.002868
universal_physical_law_claim_allowed = false
folding_problem_solved = false
```

Per-row audit:

```text
1MBN:A  selected=23   precision=0.782609  recall=0.045570  long_range_recall=0.000000  f1=0.086124
2LZM:A  selected=82   precision=0.402439  recall=0.077103  long_range_recall=0.000000  f1=0.129412
1TEN:A  selected=18   precision=0.055556  recall=0.004167  long_range_recall=0.000000  f1=0.007752
1CSP:A  selected=14   precision=0.214286  recall=0.016949  long_range_recall=0.000000  f1=0.031414
1TIM:A  selected=120  precision=0.166667  recall=0.027586  long_range_recall=0.044164  f1=0.047337
1PGA:A  selected=20   precision=0.350000  recall=0.051095  long_range_recall=0.000000  f1=0.089172
4AKE:A  selected=107  precision=0.252336  recall=0.048301  long_range_recall=0.015544  f1=0.081081
1CLL:A  selected=16   precision=0.875000  recall=0.042813  long_range_recall=0.000000  f1=0.081633
```

Claim rejection:

```text
global_folding_claim_rejected_five_axis_gate_failed_for_rows:1MBN:A,2LZM:A,1TEN:A,1CSP:A,1TIM:A,1PGA:A,4AKE:A,1CLL:A
```

## Interpretation

This layer proves the engineering gap is now explicitly represented, but it also proves the current sequence-only physics is not the missing universal law.

The strongest signal is local/chemistry precision in some rows. The failure is long-range recall and global contact-map completeness. That matches the diagnosis: sequence-only local physics can identify plausible contacts, but it does not supply the learned/evolutionary/global geometry needed to choose the right long-range cooperative contact network.

## Honest conclusion

Not cracked.

What is now better:

- the missing axes are wired;
- the test is native-free before selection;
- the failure is quantified instead of guessed;
- the claim gate refuses a solved-folding statement;
- the run path is bounded and has a CPU-worker cleanup script.

The next real research move is not another threshold tweak. It is to make cooperativity global: a factor-graph/ensemble constraint system where contacts are selected as a jointly consistent network, then falsified against matched long-range controls.
