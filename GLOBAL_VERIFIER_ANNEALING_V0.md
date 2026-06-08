# GLOBAL_VERIFIER_ANNEALING_V0

## Purpose

This batch tests the next falsifiable hypothesis after geodesic-consensus filtering failed:

> A path-level geodesic filter is not enough; the missing layer may be a native-free global verifier that scores whole candidate C-alpha structures/contact fields before native audit.

The run stays structure-first:

```text
safe DCA/geodesic restraints
→ multiple whole-structure C-alpha candidates
→ native-free global verifier
→ selected structure/contact map frozen
→ native coordinate contacts attached only after selection for audit
```

## Included verifier terms

- Statistical-potential surrogate from sequence chemistry and final geometry.
- Degree-distribution coherence.
- Loop/patch contact coherence.
- Restraint satisfaction.
- Compactness and collision penalties.
- Monte-Carlo-style deterministic multi-candidate global selection.

No AlphaFold, ESMFold, template structure, native contacts, or native coordinates are used before selection.

## Bounded validation

```text
focused tests: 3 passed in 1.18s
bounded runner: completed
full_suite_run = false
gif_generation_used = false
predictor_subprocesses_used = false
atomistic_md_engine_used = false
```

## Main 8-row result

```text
mean_precision = 0.109675
mean_recall = 0.104113
mean_long_range_recall = 0.058541
mean_F1 = 0.106693
mean_selected_global_score = 0.608361
mean_statistical_potential_score = 0.691851
mean_degree_coherence_score = 0.923557
mean_restraint_satisfaction_score = 0.879919

global_verifier_claim_allowed = false
universal_physical_law_claim_allowed = false
folding_problem_solved = false
```

## 4AKE result

```text
safe_restraint_count = 608
raw_extracted_contact_count = 2449
globally_filtered_contact_count = 524
globally_filtered_long_range_contact_count = 279
precision = 0.158397
recall = 0.148479
long_range_recall = 0.041451
F1 = 0.153278
solved = false
```

## Interpretation

The global verifier can produce internally coherent-looking maps:

```text
4AKE degree_coherence_score = 0.936177
4AKE restraint_satisfaction_score = 0.823982
```

But that coherence does not correspond to the native topology. This is the important negative result:

```text
native-free global coherence surrogate ≠ true global geometry prior
```

The verifier avoids the previous contact explosion, improving 4AKE precision over loose geodesic consensus (`0.093 → 0.158`), but it also collapses true long-range recovery (`0.181 → 0.041`).

So the current failure mode is not just missing a global score. It is missing a global score that is calibrated to real fold topology rather than generic compact/coherent protein-like geometry.
