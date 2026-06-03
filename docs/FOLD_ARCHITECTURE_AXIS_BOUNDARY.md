# Fold Architecture Axis Boundary

This layer adds sequence-only architecture-axis evidence adjudication to the
locked 50-row folding benchmark. It is not a fold-class predictor and it does
not use labels, PDB truth, or structure-derived classes before prediction.

## Why This Exists

The axis-profile layer recovered partial biological truth safely:

```text
collapsed_class_coverage = 0.28
axis_profile_coverage = 0.86
architecture_axis_coverage = 0.14
axis_profile_same_axis_conflict_count = 0
unsafe_class_recovery_count = 0
```

The remaining weak point was architecture. The system could often say an order
axis, but it still needed a stricter way to decide whether a sequence supports:

```text
compact_single_domain
multidomain_or_segmented
repeat_like
fragment_scope
unknown
```

## Evidence Packet

Each row now gets an architecture evidence packet:

```text
row_id
sequence_hash
length_band
segmentation_pressure
repeat_pressure
compact_domain_pressure
fragment_scope_pressure
multidomain_pressure
membrane_segmentation_pressure
disorder_segmentation_pressure
evidence_strength
architecture_axis_prediction
architecture_axis_confidence
architecture_axis_claim_allowed
architecture_axis_abstention_reason
```

The CSV and dashboard do not export raw protein sequences.

## Claim Rules

`fragment_scope` is allowed only for short sequence scope where the row is
better treated as a solved fragment than as global protein architecture.

`compact_single_domain` is allowed only when small-to-medium length, compact
closure, low segmentation, low repeat pressure, low disorder pressure, and low
membrane segmentation all agree.

`multidomain_or_segmented` is allowed only for large modular evidence with
boundary pressure that is not better explained as repeat-like, membrane-like,
or disorder segmentation.

`repeat_like` is allowed only for strong stable repeat evidence. In the current
50-row run, no row clears that threshold.

The external-safe repair layer adds an opt-in repeat-vs-compact quarantine for
the external holdout. It refuses `repeat_like` when small/medium compact-domain
closure, weak composition contrast, weak recurrence support, and
hydrophobic-periodicity-only risk explain the signal better than repeat
architecture.

That guard is documented in `docs/EXTERNAL_AXIS_REPAIR_BOUNDARY.md` and is used
by `scripts/run_external_axis_repair_benchmark.py`.

## Forbidden Moves

```text
Do not use labels or PDB truth before prediction.
Do not treat secondary structure as architecture evidence.
Do not claim multidomain just because the protein is long.
Do not claim multidomain just because the sequence is mixed alpha/beta.
Do not claim multidomain for membrane helix segmentation.
Do not recover global fold-class coverage.
Do not override abstention guards.
Do not claim protein folding solved.
Do not design molecules, drugs, binders, therapies, or dosing logic.
```

## Run

```bash
python3 scripts/run_architecture_axis_benchmark.py
```

It writes:

```text
first_contact_clean_pharmacotopology_layer_run/real_folding_50_architecture_axis_report.json
first_contact_clean_pharmacotopology_layer_run/real_folding_50_architecture_axis_rows.csv
first_contact_clean_pharmacotopology_layer_run/real_folding_50_architecture_axis_conflicts.csv
first_contact_clean_pharmacotopology_layer_run/real_folding_50_architecture_axis_abstentions.csv
first_contact_clean_pharmacotopology_layer_run/real_folding_50_architecture_axis_dashboard.html
first_contact_clean_pharmacotopology_layer_run/real_folding_50_architecture_axis_certificate.json
```

## Current Result

```text
benchmark_size = 50
previous_profile_architecture_axis_coverage = 0.14
architecture_axis_coverage = 0.32
architecture_axis_claim_allowed_count = 16
architecture_axis_abstained_count = 34
architecture_axis_conflict_count = 0
architecture_axis_safe_claim_count = 16
architecture_axis_same_axis_conflict_count = 0
architecture_axis_accuracy = 1.0
fragment_scope_detected_count = 7
compact_single_domain_claim_count = 6
multidomain_claim_count = 3
multidomain_abstained_count = 9
repeat_like_claim_count = 0
architecture_claim_without_secondary_leakage_count = 0
architecture_claim_without_label_leakage_count = 0
architecture_axis_recovered_from_profile_count = 7
global_fold_class_claim_allowed = false
axis_profile_claim_allowed = true
claim_allowed = false
folding_problem_solved = false
```

This improves architecture-axis coverage from `0.14` to `0.32` without adding
same-axis conflicts.

## Important Abstentions

Rhodopsin has membrane-like segmentation pressure, so the architecture layer
does not claim multidomain from that signal:

```text
pdb_1F88_A_rhodopsin
architecture_axis_prediction = unknown
architecture_axis_abstention_reason = membrane segmentation can reflect helices rather than global domains
```

GroEL has strong modular pressure, but repeat/segmentation signals are mixed.
The layer refuses the multidomain claim:

```text
pdb_1AON_A_groel
architecture_axis_prediction = unknown
architecture_axis_abstention_reason = repeat/segmentation signals are mixed
```

Those abstentions are the point of the batch. Architecture evidence is allowed
to improve coverage only where it stays separated from secondary structure,
length-only inference, membrane topology, and label leakage.
