# External Fold-Family Holdout Boundary

This layer adds a locked 100-row external fold-family holdout benchmark for the
current safe folding-axis stack. It is a falsification batch, not a repair
batch.

The runner evaluates the existing stack as-is:

```text
regime routing
fold-axis safety guards
axis-profile recovery
architecture-axis adjudication
```

No thresholds, guards, architecture rules, profile recovery rules, or predictor
decision paths are changed in this batch.

## Locked Holdout

The locked data file is:

```text
data/folding_benchmarks_external_fold_family_100.locked.json
```

It contains 100 rows with this distribution:

```text
20 alpha-rich / helical soluble rows
20 beta-rich rows
20 alpha-beta mixed rows
15 membrane-like rows
15 disordered/flexible or mixed-order rows
10 architecture stress rows
```

Each row stores:

```text
row_id
source_id
source_kind
sequence
sequence_sha256
length
external_family_id
external_family_name
external_family_group
holdout_split
truth_axes
truth_scope
evidence_notes
```

Raw sequences are present only in the locked data file. The generated report,
CSV files, dashboard, and certificate expose sequence hashes and metadata, not
raw sequences.

## Non-Overlap Guard

The holdout runner checks the locked holdout against the development 50:

```text
same row_id
same source_id
same sequence_sha256
same external_family_id when present
```

Current result:

```text
holdout_row_count = 100
holdout_unique_sequence_count = 100
holdout_unique_family_count = 23
development_overlap_count = 0
development_sequence_overlap_count = 0
development_family_overlap_count = 0
holdout_lock_valid = true
holdout_non_overlap_valid = true
```

## Run

```bash
python3 scripts/run_external_fold_family_holdout_benchmark.py
```

It writes:

```text
first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_report.json
first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_rows.csv
first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_family_summary.csv
first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_axis_conflicts.csv
first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_abstentions.csv
first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_failure_cohorts.csv
first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_dashboard.html
first_contact_clean_pharmacotopology_layer_run/external_fold_family_100_certificate.json
```

## Current Holdout Result

```text
collapsed_class_coverage = 0.37
axis_profile_coverage = 0.78
secondary_axis_coverage = 0.05
architecture_axis_coverage = 0.17
order_axis_coverage = 0.76
environment_axis_coverage = 0.17

axis_profile_same_axis_conflict_count = 22
architecture_axis_same_axis_conflict_count = 14
forced_same_axis_conflict_count = 22
high_confidence_wrong_count_after_axis_scoring = 17

safe_axis_claim_count = 79
unsafe_axis_claim_count = 36
safe_axis_recovered_count = 78
unsafe_class_recovery_count = 0
guard_override_count = 0

family_level_failure_count = 14
family_level_abstention_count = 22
family_level_conflict_count = 13
family_generalization_status = conflicts_detected

global_fold_class_claim_allowed = false
axis_profile_claim_allowed = true
architecture_axis_claim_allowed = false
folding_problem_solved = false
claim_allowed = false
```

This is not a failure of the batch. It is the point of the batch. The current
stack generalizes enough to make partial axis claims, but the holdout exposes
unsafe axis conflicts and a weaker architecture-axis surface.

## Dominant Failure Signals

The first repair candidates should come from the generated failure cohorts, not
from guessing. The current holdout highlights:

```text
architecture_axis_conflict
order_axis_conflict
beta_rich_low_complexity_confusion
alpha_beta_boundary_ambiguity
repeat_like_under_claimed
large_chain_multidomain_under_claimed
```

Those are cohorts, not fixes. This commit does not repair them.

The next layer, documented in `docs/EXTERNAL_AXIS_REPAIR_BOUNDARY.md`, repairs
the two unsafe overclaim classes by abstention:

```text
disordered_flexible order-axis projection -> quarantine unless disorder evidence is strong
repeat_like architecture projection -> quarantine unless recurrence evidence is strong
```

The holdout artifacts remain useful as the before picture. The repair artifacts
record the after picture.

## Boundary

This layer does not:

```text
tune thresholds
change architecture rules
change regime rules
change safety guards
add row-specific exceptions
recover global fold-class coverage
generate protein sequences as a design task
design molecules, binders, drugs, therapies, or dosing logic
claim protein folding is solved
```

It only answers:

```text
The current stack was tested against unseen fold-family rows without tuning.
Here is where it generalizes.
Here is where it fails.
```
