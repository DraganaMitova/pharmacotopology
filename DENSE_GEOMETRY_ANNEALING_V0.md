# Dense Geometry Annealing V0

This batch tests the next honest hypothesis after the native-free geometry-field
result: whether sparse DCA anchors can be amplified into a dense, accurate global
geometry field without ESMFold/AlphaFold, templates, a trained structure model,
or native/coordinate truth before selection.

## Boundary

- No GIF generation.
- No internet.
- No predictor subprocesses.
- No full pytest suite.
- Native/contact truth is attached only after selection for audit metrics.
- CPU cleanup is run before focused tests/runners and after final validation.

## Implemented

- `src/pharmacotopology/folding_dense_geometry_annealing.py`
- `scripts/run_dense_geometry_annealing_v0.py`
- `tests/test_dense_geometry_annealing_v0.py`
- Artifacts under `first_contact_clean_pharmacotopology_layer_run/dense_geometry_annealing_v0/`

## Modes tested

### 1. Pure sequence coarse annealing

Sequence priors + bounded coarse chain relaxation, without DCA/MSA anchors.

Result:

```text
mean_precision = 0.304095
mean_recall = 0.238340
mean_long_range_recall = 0.000000
mean_f1 = 0.266799
folding_problem_solved = false
universal_physical_law_claim_allowed = false
```

Interpretation: coarse chain relaxation can produce a denser local/medium contact
field, but it does not recover the long-range fold topology.

### 2. External DCA compact annealing

Safe external DCA anchors pull a coarse C-alpha trace through bounded spring
relaxation and steric/compactness constraints.

Result:

```text
mean_precision = 0.262841
mean_recall = 0.205385
mean_long_range_recall = 0.000000
mean_f1 = 0.230234
folding_problem_solved = false
universal_physical_law_claim_allowed = false
```

Interpretation: compact annealing increases density and some contact recovery,
but it collapses toward local/medium compactness and loses long-range topology.

### 3. External DCA multifield annealing

Union of compact annealing plus a long-range anchor-temperature channel.

Result:

```text
mean_precision = 0.173463
mean_recall = 0.269943
mean_long_range_recall = 0.246702
mean_f1 = 0.210899
folding_problem_solved = false
universal_physical_law_claim_allowed = false
```

4AKE row:

```text
precision = 0.213785
recall = 0.327370
long_range_recall = 0.088083
f1 = 0.258657
folding_problem_solved = false
```

Interpretation: long-range exploration recovers more long-range signal than the
compact channel, but it is not accurate enough.  The field becomes denser, not
accurate.

## Final diagnosis

Dense geometry generation has split into two failure modes:

1. **Compact channel:** better precision/F1, but long-range recall collapses.
2. **Long-range channel:** better long-range recall, but precision and F1 fall.

This means the missing object is not merely "more density".  The missing object
is a dense field that is simultaneously accurate and globally ordered.  ESMFold
provides both density and ordering; native-free DCA/annealing does not.

## Claim decision

```text
dense_geometry_annealing_claim_allowed = false
universal_physical_law_claim_allowed = false
folding_problem_solved = false
```

This is not a crack.  It is a narrower falsification: sparse DCA anchors plus
bounded coarse annealing can create signal, but not a dense accurate global
geometry prior comparable to the learned ESMFold prior.
