# Current Project Status

This is the current truth of the repository after the contact-topology repair
and native-gap analysis milestone.

## Latest Safety Baseline

```text
latest recorded safety commit = 019378b Add visual folding mechanism workbench
current target = Add contact-topology repair and native-gap analysis
```

The current target does not unlock fold-class claims. It explains visible
contact failures after prediction and tests sequence-only contact repairs.

## Internal 50-Row Status

The internal 50-row benchmark is closed as a development safety benchmark:

```text
axis_profile_coverage = 0.86
axis_profile_same_axis_conflict_count = 0
forced_same_axis_conflict_count = 0
high_confidence_wrong_count_after_axis_scoring = 0
global_fold_class_claim_allowed = false
folding_problem_solved = false
claim_allowed = false
```

This does not mean folding is solved. It means the old internal same-axis
safety failures are closed.

## External 100-Row Status

The external holdout before repair exposed unsafe generalization:

```text
axis_profile_coverage = 0.78
axis_profile_same_axis_conflict_count = 22
architecture_axis_same_axis_conflict_count = 14
unsafe_axis_claim_count = 36
high_confidence_wrong_count_after_axis_scoring = 17
```

The external-safe repair closes those unsafe claims by abstention:

```text
post_repair_axis_profile_coverage = 0.55
post_repair_axis_profile_same_axis_conflict_count = 0
post_repair_architecture_axis_same_axis_conflict_count = 0
post_repair_unsafe_axis_claim_count = 0
post_repair_high_confidence_wrong_count_after_axis_scoring = 0
order_axis_folded_mimic_quarantined_count = 22
repeat_compact_ambiguity_quarantined_count = 14
coverage_loss_from_external_safety = 0.23
external_safety_repair_successful = true
```

Coverage loss is the intended cost of the repair.

## Visual 12-Row Mechanism Status

The visual workbench is a locked contact-map mechanism layer:

```text
visual_artifacts_generated_for_rows = 12
contact_map_f1_computed_count = 12
visible_partial_success_count = 5
visible_failure_count = 7
mean_contact_map_f1 = 0.104031
native_truth_used_before_prediction = false
raw_sequence_exposed = false
global_folding_claim_allowed = false
folding_problem_solved = false
```

The visual workbench does not claim atomistic folding. It shows coarse native
contact targets, sequence-only contact candidates, overlays, closure curves,
coarse trajectories, and failure cohorts.

## Contact-Repair 12-Row Status

The repair workbench adds native-gap analysis after prediction and sequence-only
long-range/beta-pairing repairs:

```text
baseline_visible_partial_success_count = 5
visible_partial_success_count = 8
baseline_visible_failure_count = 7
visible_failure_count = 4
baseline_mean_contact_map_f1 = 0.104031
repaired_mean_contact_map_f1 = 0.263729
long_range_contact_recall_delta = 0.472223
beta_pairing_contact_recall_delta = 0.6541
premature_compaction_count = 5
visual_failure_cohort_count = 4
native_truth_used_before_prediction = false
native_truth_used_before_repair = false
raw_sequence_exposed = false
global_folding_claim_allowed = false
folding_problem_solved = false
```

Remaining repaired failure mechanisms:

```text
bad_beta_pairing = 1
disorder_over_collapse = 1
membrane_mis_topology = 1
premature_compaction = 1
```

This is a contact-topology repair milestone, not a solved folding milestone.

## Canonical Active Stack

```text
regime routing
fold-axis adjudication
axis-profile recovery
architecture-axis adjudication
external holdout
external-safe axis repair
visual contact-map mechanism workbench
contact-topology repair and native-gap analysis
```

Canonical runners:

```text
scripts/run_fold_axis_adjudication_benchmark.py
scripts/run_fold_axis_profile_benchmark.py
scripts/run_architecture_axis_benchmark.py
scripts/run_external_fold_family_holdout_benchmark.py
scripts/run_external_axis_repair_benchmark.py
scripts/run_visual_folding_mechanism_benchmark.py
scripts/run_contact_topology_repair_benchmark.py
```

Active artifacts:

```text
real_folding_50_axis_*
real_folding_50_axis_profile_*
real_folding_50_architecture_axis_*
external_fold_family_100_*
external_axis_repair_*
visual_mechanism_12_*
visuals/*/*
contact_topology_repair_12_*
contact_repair_visuals/*/*
```

Archived legacy artifacts:

```text
first_contact_clean_pharmacotopology_layer_run/archived_legacy/
```

## Known Limits

```text
protein folding is not solved
global fold-class coverage is still locked
axis profile coverage is partial
secondary-structure evidence is still thin on the external holdout
architecture evidence is conservative after repair
visual contact prediction is coarse and low-F1
contact repair is benchmark-local and still coarse
disorder, membrane, and mixed beta/compaction failures remain visible
legacy feature modules still contain useful low-level primitives
```

The current state is useful as a safe falsification, abstention, visual
mechanism-inspection, and contact-repair stack, not as a predictive folding
system.

## Next Allowed Research Target

Do not add more broad biology claims before the canonical surface stays clean.

The next research target can be:

```text
visual contact mechanism calibration
```

Allowed shape:

```text
improve contact evidence packets
preserve truth-after-prediction scoring
keep global fold-class claims locked
track coverage cost
record uncertainty, failure cohorts, and visible contact-map drift
reduce premature compaction without hiding disorder/membrane failures
```

Forbidden shape:

```text
force beta/mixed/ordered labels to improve coverage
use holdout truth during prediction
recover collapsed-class coverage
claim folding solved
generate proteins, molecules, binders, therapies, or dosing logic
```

## Clean Zip

For sharing, export from tracked files only:

```bash
git archive --format=zip --output pharmacotopology-clean.zip HEAD
```

This avoids local caches and runtime directories.
