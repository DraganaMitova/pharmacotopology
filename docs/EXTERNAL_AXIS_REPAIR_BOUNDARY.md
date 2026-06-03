# External Axis Repair Boundary

This layer responds to the locked external fold-family holdout by adding
external-safe axis quarantines. It is a safety calibration layer, not an
accuracy or folding-solved layer.

The holdout exposed two overclaims:

```text
disordered_flexible was projected too easily onto order_axis
repeat_like architecture was claimed from weak/generic periodicity
```

The repair is intentionally conservative:

```text
unsafe axis claim -> unknown / abstention
```

It does not flip rows to better-looking labels.

## Run

```bash
python3 scripts/run_external_axis_repair_benchmark.py
```

It writes:

```text
first_contact_clean_pharmacotopology_layer_run/external_axis_repair_report.json
first_contact_clean_pharmacotopology_layer_run/external_axis_repair_rows.csv
first_contact_clean_pharmacotopology_layer_run/external_axis_repair_conflict_delta.csv
first_contact_clean_pharmacotopology_layer_run/external_axis_repair_abstention_delta.csv
first_contact_clean_pharmacotopology_layer_run/external_axis_repair_quarantine_rows.csv
first_contact_clean_pharmacotopology_layer_run/external_axis_repair_family_summary.csv
first_contact_clean_pharmacotopology_layer_run/external_axis_repair_dashboard.html
first_contact_clean_pharmacotopology_layer_run/external_axis_repair_certificate.json
```

## Order-Axis Quarantine

A forced `disordered_flexible` class is no longer sufficient by itself to make
an order-axis claim.

The repair layer allows:

```text
order_axis = disordered_flexible
```

only when disorder evidence is strong and folded-domain mimic evidence is weak:

```text
high disorder pressure
non-trivial disorder-run evidence
strong breaker / low-complexity evidence
weak compact closure
weak beta pairing support
weak folded beta/mixed mimic pressure
```

Otherwise the row is quarantined:

```text
order_axis = mixed_or_uncertain
order_axis_claim_allowed = false
order_axis_abstention_reason = external_order_axis_folded_mimic_quarantine
```

The repair does not convert these rows to `ordered`. That would be a hidden
truth repair. The safe move is abstention.

## Architecture-Axis Quarantine

The architecture repair blocks `repeat_like` claims when the evidence looks
like compact-domain hydrophobic periodicity rather than repeat architecture.

The repeat-vs-compact guard looks for:

```text
small/medium length
repeat pressure close to compact-domain pressure
non-trivial long-range closure
weak composition contrast
weak recurrence support
hydrophobic-periodicity-only risk
```

When this pattern appears, the repaired architecture packet says:

```text
architecture_axis_prediction = unknown
architecture_axis_claim_allowed = false
architecture_axis_abstention_reason = repeat_compact_single_domain_ambiguity_quarantine
```

Repeat-like requires recurrence, not just periodic hydrophobic clusters.

## Current Repair Result

```text
pre_repair_axis_profile_coverage = 0.78
post_repair_axis_profile_coverage = 0.55

pre_repair_axis_profile_same_axis_conflict_count = 22
post_repair_axis_profile_same_axis_conflict_count = 0

pre_repair_architecture_axis_same_axis_conflict_count = 14
post_repair_architecture_axis_same_axis_conflict_count = 0

pre_repair_unsafe_axis_claim_count = 36
post_repair_unsafe_axis_claim_count = 0

pre_repair_high_confidence_wrong_count_after_axis_scoring = 17
post_repair_high_confidence_wrong_count_after_axis_scoring = 0

order_axis_folded_mimic_quarantined_count = 22
repeat_compact_ambiguity_quarantined_count = 14
coverage_loss_from_external_safety = 0.23
external_safety_repair_successful = true
```

Coverage dropped. That is the intended cost. The repair closes unsafe external
axis claims by refusing to overstate what sequence-only evidence can support.

## Certificate Boundary

The generated certificate asserts:

```json
{
  "external_axis_repair_complete": true,
  "legacy_axis_artifacts_reproducible": true,
  "legacy_axis_profile_artifacts_reproducible": true,
  "order_axis_folded_mimic_quarantine_active": true,
  "repeat_compact_ambiguity_quarantine_active": true,
  "post_repair_axis_profile_same_axis_conflict_count": 0,
  "post_repair_architecture_axis_same_axis_conflict_count": 0,
  "post_repair_unsafe_axis_claim_count": 0,
  "post_repair_high_confidence_wrong_count_after_axis_scoring": 0,
  "unsafe_class_recovery_count": 0,
  "guard_override_count": 0,
  "global_fold_class_claim_allowed": false,
  "folding_problem_solved": false,
  "claim_allowed": false
}
```

## Boundary Rules

This layer does not:

```text
turn beta/mixed rows into ordered claims
tune directly to row IDs
use truth_axes before prediction
improve collapsed_class_coverage
recover global fold class
claim folding solved
export raw sequences in generated artifacts
generate proteins, molecules, binders, therapies, or dosing logic
```

It only says:

```text
The external holdout showed unsafe axis overclaims.
Those overclaims are now quarantined.
The remaining useful output is axis-safe partial truth with abstention.
```
