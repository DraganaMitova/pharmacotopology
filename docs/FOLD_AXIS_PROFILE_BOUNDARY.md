# Fold Axis Profile Boundary

This layer adds axis-safe profile coverage recovery to the locked 50-row folding
benchmark. It does not create new global fold-class predictions. It only emits
partial axis claims where the existing sequence-only gate path already supports
that axis.

## Why This Exists

The previous safety repair closed forced same-axis conflicts:

```text
forced_same_axis_conflict_count = 0
forced_order_axis_conflict_count = 0
forced_secondary_axis_conflict_count = 0
high_confidence_wrong_count_after_axis_scoring = 0
```

That closure was correct, but it lowered collapsed class coverage. The next
safe move is not to force more classes. It is to say which axes are supported
and which remain unknown:

```text
I know this axis.
I do not know that axis.
I refuse the global fold class.
```

## Axes

```text
secondary_structure_axis:
  alpha_rich
  beta_rich
  alpha_beta_mixed
  weak_or_unknown

architecture_axis:
  compact_single_domain
  multidomain_or_segmented
  repeat_like
  fragment_scope
  unknown

order_axis:
  ordered
  disordered_flexible
  mixed_or_uncertain

environment_axis:
  soluble_like
  membrane_like
  unknown
```

## Profile Rules

The profile layer can recover:

```text
secondary axis:
  only from specific alpha-periodic or beta-pairing secondary gates

architecture axis:
  only fragment_scope from the small_peptide_or_fragment regime

order axis:
  from source forced folded/disordered classes
  from foldable compact gate evidence
  from the folded-domain-mimic guard as ordered
  from the alpha/mixed ambiguity guard as ordered

environment axis:
  only membrane_like from the membrane_like regime
```

The profile layer does not recover:

```text
global fold class
compact vs multidomain architecture
repeat architecture
soluble_like environment
raw protein sequences
protein designs
molecules
drug or clinical claims
```

## Run

```bash
python3 scripts/run_fold_axis_profile_benchmark.py
```

It writes:

```text
first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_profile_report.json
first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_profile_rows.csv
first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_profile_abstentions.csv
first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_profile_recovery_candidates.csv
first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_profile_dashboard.html
first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_profile_certificate.json
```

## Current Result

```text
benchmark_size = 50
collapsed_class_coverage = 0.28
axis_profile_coverage = 0.86
secondary_axis_coverage = 0.16
architecture_axis_coverage = 0.14
order_axis_coverage = 0.84
environment_axis_coverage = 0.04
safe_axis_recovered_count = 39
safe_axis_recovered_row_count = 29
unsafe_class_recovery_count = 0
guard_override_count = 0
forced_same_axis_conflict_count = 0
axis_profile_same_axis_conflict_count = 0
high_confidence_wrong_count_after_axis_scoring = 0
global_fold_class_claim_allowed = false
axis_profile_claim_allowed = true
folding_problem_solved = false
```

The result is intentionally asymmetric. It recovers many order-axis claims,
some fragment-scope claims, two secondary-axis claims, and two membrane-axis
claims. It does not try to repair compact-vs-multidomain architecture yet.

## Guard Cases

```text
pdb_1BFG_A_basic_fgf
profile_global_fold_class = insufficient_topology_evidence
profile_secondary_structure_axis = weak_or_unknown
profile_order_axis = ordered
safe_axis_recovered_axes = order_axis
guard_override = false

pdb_1BTA_A_barstar
profile_global_fold_class = insufficient_topology_evidence
profile_secondary_structure_axis = weak_or_unknown
profile_order_axis = ordered
safe_axis_recovered_axes = order_axis
guard_override = false
```

These are not solved fold classifications. They are narrower, axis-level
statements that preserve the abstention boundary.

## Interpretation

This layer makes the project more useful without becoming less honest:

```text
BAD:
  abstained row -> force alpha_rich

GOOD:
  abstained row -> secondary_axis = alpha_rich
                 -> architecture_axis = unknown
                 -> global_class_claim_allowed = false
```

The next evidence target became `docs/FOLD_ARCHITECTURE_AXIS_BOUNDARY.md`,
which adds a separate sequence-only architecture evidence packet while keeping
global fold-class claims refused.
