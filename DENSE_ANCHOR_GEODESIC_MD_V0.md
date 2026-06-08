# Dense Anchor Geodesic MD V0

This batch tests the requested structure-first continuation:

```text
safe DCA anchors
+ sequence-closure / MI-surrogate support
+ beta-bridge support
+ anchor chaining
+ geodesic interpolation
-> dense restraints
-> bounded coarse C-alpha global relaxation
-> extract contacts from final geometry
-> native audit only after selection
```

It is not an OpenMM/GROMACS atomistic run. It is a dependency-free, bounded
coarse-grain falsifier designed to avoid hanging subprocesses and to keep the CPU
clean. No ESMFold, AlphaFold, template, native coordinate, or native contact
truth is used before selection.

Important limitation: raw MSA files are not bundled in this benchmark archive, so
this batch does not perform a fresh mutual-information computation. The MI
channel is recorded as a surrogate based on safe DCA raw/APC support plus
sequence-only closure agreement.

Claim gate remains the same: precision, recall, and long-range recall must all
reach 0.70 before any solve claim is allowed.

## Validated result

```text
mean_dense_anchor_count = 551.125
min_dense_anchor_count = 224
mean_precision = 0.081669
mean_recall = 0.627898
mean_long_range_recall = 0.491544
mean_f1 = 0.144079
folding_problem_solved = false
universal_physical_law_claim_allowed = false
claim_rejection_reason = dense_anchor_geodesic_claim_rejected_precision_below_0_70
```

For 4AKE:A specifically:

```text
dense_anchor_count = 822
precision = 0.095254
recall = 0.524150
long_range_recall = 0.191710
F1 = 0.161210
solved = false
```

Interpretation: dense anchors and geodesic interpolation increased global coverage
on the full 8-row benchmark, but precision collapsed. The method creates a dense
restraint field, but the field is not accurate enough to replace ESMFold-like
global geometry.
