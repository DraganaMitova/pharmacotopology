# Fold Axis Truth Boundary

This layer adds orthogonal fold-axis adjudication to the locked 50-row folding
benchmark. It is not a new folding predictor. It does not retune the regime
router or hierarchy gates. It only asks whether the existing single-class
benchmark labels are collapsing different kinds of truth into one mutually
exclusive class.

## Why This Exists

The previous 50-row regime layer improved safety by reducing high-confidence
wrong predictions from 6 to 0, but it also revealed a benchmark ontology
problem. The single fold class mixes several axes:

```text
secondary-structure composition
architecture / segmentation
ordering / stability
environment / topology regime
```

Those axes are not mutually exclusive. A protein can be membrane-like,
alpha-rich, ordered, and segmental at the same time. So a row where one channel
says `multidomain_boundary` and another says `alpha_rich` may be an orthogonal
axis projection, not a true contradiction.

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

## Boundary

Prediction-side axes are projected only from existing sequence-only outputs:
the regime-routed prediction, the hierarchy output, and sequence-derived regime
signals already emitted before scoring.

The layer does not use these truth channels before prediction:

```text
external labels
PDB coordinates
CATH / SCOP classes
DisProt truth labels
structure-derived fold classes
reference topology signatures
```

Truth axes are adjudicated only after prediction from structure-derived class,
external label class, structure evidence metadata, and environment evidence
where available. Raw protein sequences are not exported in the CSV or
dashboard.

## Run

```bash
python3 scripts/run_fold_axis_adjudication_benchmark.py
```

It writes:

```text
first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_adjudication_report.json
first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_rows.csv
first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_conflicts.csv
first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_manual_review.csv
first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_confusion_matrices.csv
first_contact_clean_pharmacotopology_layer_run/real_folding_50_axis_dashboard.html
```

## Current Result

```text
benchmark_size = 50
single_class_taxonomy_collapse_detected = true
structure_label_disagreement_count = 17
orthogonal_axis_disagreement_count = 10
true_same_axis_conflict_count = 8
structure_label_same_axis_conflict_count = 7
forced_same_axis_conflict_count = 0
forced_order_axis_conflict_count = 0
forced_secondary_axis_conflict_count = 0
abstained_axis_conflict_count = 8
regime_axis_conflict_count = 8
folded_domain_mimic_abstained_count = 1
secondary_axis_ambiguity_abstained_count = 2
coverage_loss_from_safety_guards = 3
axis_unscorable_count = 147
high_confidence_wrong_count_after_axis_scoring = 0
artifact_reproducible = true
claim_allowed = false
folding_problem_solved = false
```

The environment axis is mostly unscorable from the current truth channels, so
`axis_unscorable_count` is expected to be high. That is useful: it shows which
claims the benchmark does not yet have evidence to score.

## Interpretation

The dashboard explicitly preserves these distinctions:

```text
Single-Class Benchmark Is Lossy
Architecture ≠ Secondary Structure
Membrane Regime ≠ Fold Class
Disorder ≠ Beta/Alpha Absence
Fragment Evidence ≠ Global Fold Truth
```

The rhodopsin-style case is now interpreted as taxonomy collapse:

```text
membrane_like + alpha_rich + multidomain/segmental
```

The two previously forced low-confidence model errors are now abstentions, not
new confident classifications:

```text
pdb_1BFG_A_basic_fgf
guard: folded_domain_mimic_disorder_conflict
predicted fold class: insufficient_topology_evidence
truth order axis: ordered

pdb_1BTA_A_barstar
guard: secondary_structure_gate:abstained_alpha_mixed_ambiguity
predicted fold class: insufficient_topology_evidence
truth secondary axis: alpha_beta_mixed
```

The next modeling work should target those areas directly, but only after
evidence improves enough to classify rather than merely abstain: folded beta
proteins that mimic disorder from sequence-only features, and weak
alpha-vs-mixed secondary-axis separation.

The immediate coverage-recovery layer after this boundary is
`docs/FOLD_AXIS_PROFILE_BOUNDARY.md`. It keeps global fold class claims refused
while emitting only safe partial axis profiles.
