# 4AKE no-AlphaFold iterative distance-geometry diffusion final attempt

This package adds a bounded Path B implementation: native-free iterative distance-geometry diffusion. It does not use AlphaFold, templates, native coordinates, or raw FASTA persistence before prediction. Native 4AKE contacts are attached only after selection for audit metrics.

## New implementation

- `src/pharmacotopology/folding_iterative_distance_geometry_diffusion.py`
- `scripts/run_4ake_iterative_distance_geometry_diffusion_v0.py`
- `tests/test_iterative_distance_geometry_diffusion_v0.py`

The method combines safe external DCA anchors, sequence-only energy/secondary-structure/degree priors, candidate closure events, an anchor-field diffusion layer, and repeated graph shortest-path/classical-MDS distance geometry. Default execution is bounded to avoid hanging; heavier mode sweep is opt-in with `--mode-sweep`.

## 4AKE result after native audit

Selected mode by native-free internal quality: `balanced_dg_diffusion`

```text
folding_problem_solved = false
native_contact_precision = 0.181818
native_contact_recall = 0.143113
contact_map_f1 = 0.160160
long_range_precision = 0.147826
long_range_recall = 0.088083
predicted_contact_count = 440
true_positive_contacts = 80
native_contact_count = 559
seed_anchor_pair_count = 214
claim_rejection_reason = native_free_iterative_diffusion_did_not_reach_0_70_precision_and_recall
```

This is an improvement over the earlier sparse no-AF path in recall, but it still does not solve 4AKE at the required 0.70 precision and 0.70 recall threshold.

## Safety / no-hang behavior

```text
gif_generation_used = false
alphafold_used_as_evidence = false
structure_template_used_as_evidence = false
internet_required = false
native_truth_used_before_selection = false
coordinate_truth_used_before_selection = false
raw_sequence_persisted = false
default mode sweep = off
```

Validation logs included:

- `ITERATIVE_DG_DIFFUSION_RUN_VALIDATION.log`
- `TEST_SUITE_ITERATIVE_DG_NOHANG_VALIDATION.log`
- `first_contact_clean_pharmacotopology_layer_run/4ake_iterative_distance_geometry_diffusion_v0/4ake_iterative_distance_geometry_diffusion_report.json`

## Commands

Run the bounded 4AKE attempt:

```bash
PYTHONPATH=src PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
python3 -u scripts/run_4ake_iterative_distance_geometry_diffusion_v0.py
```

Run the no-hang default test suite:

```bash
PYTHONPATH=src PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
PHARMACOTOPOLOGY_TEST_TIMEOUT_SECONDS=70 \
timeout 180s python3 -m pytest -q -vv
```

Observed validation:

```text
5 passed, 2 skipped in 13.86s
```

## Scientific conclusion

The full iterative distance-geometry diffusion path was implemented and tested. It did not solve 4AKE without AlphaFold or an MSA-free learned structure predictor. The system correctly abstains instead of claiming a solved folding benchmark.
