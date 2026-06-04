# Current Project Status

This is the current truth of the repository after the real-coordinate visual
contact benchmark milestone.

## Latest Safety Baseline

```text
latest recorded safety commit = e8cc5eb Audit visual mechanism claims and freeze toy-contact benchmark
current target = Add real-coordinate visual contact benchmark
```

The current target does not add another repair heuristic. It upgrades the
visual proof target from toy locked contacts to native contact maps derived
from locked C-alpha coordinate traces.

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

## Visual Mechanism Audit Status

The audit compares the baseline visual workbench and repaired contact workbench:

```text
visual_12_is_toy_benchmark = true
coarse_native_contacts_only = true
baseline_visible_partial_success_count = 5
repaired_visible_partial_success_count = 8
baseline_visible_failure_count = 7
repaired_visible_failure_count = 4
hardcoded_beta_registry_pair_templates_detected = true
contact_repair_overfit_risk_reported = true
overfit_risk_row_count = 5
beta_template_success_gain_row_count = 3
artifact_reproducible = true
clean_archive_required = true
finder_zip_allowed = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

The audit is intentionally conservative. It says the repair improvement is
real on the tiny locked contact-map benchmark, but not proof of folding
mechanism discovery.

## Real-Coordinate Visual 8-Row Status

The coordinate benchmark uses eight RCSB PDB chains from the locked real-10
reference set. Native contacts are derived from locked C-alpha coordinate
traces at runtime:

```text
real_coordinate_native_contacts_extracted = true
toy_locked_contact_targets_used = false
coarse_ca_only = true
full_atomic_folding_available = false
benchmark_size = 8
contact_map_f1_computed_count = 8
mean_contact_map_f1 = 0.051198
max_contact_map_f1 = 0.089776
mean_native_contact_precision = 0.158448
mean_native_contact_recall = 0.030616
visible_partial_success_count = 3
visible_failure_count = 5
native_truth_used_before_prediction = false
coordinate_truth_used_before_prediction = false
raw_sequence_exposed = false
repair_heuristic_applied = false
mechanism_discovery_claim_allowed = false
global_folding_claim_allowed = false
folding_problem_solved = false
```

This is a stricter contact-map proof target than the 12-row toy visual
benchmark. The lower score is not a regression; it is the point of using a more
realistic target.

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
visual mechanism claim audit
real-coordinate visual contact benchmark
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
scripts/run_visual_mechanism_audit.py
scripts/run_real_coordinate_visual_benchmark.py
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
visual_mechanism_audit_*
real_coordinate_visual_8_*
real_coordinate_visuals/*/*
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
visual 12 is toy/coarse/internal
beta-registry repair has explicit overfit risk
real-coordinate visual 8 is C-alpha only
real-coordinate mean contact-map F1 is low
coordinate-native contacts are not atomistic folding truth
legacy feature modules still contain useful low-level primitives
```

The current state is useful as a safe falsification, abstention, visual
mechanism-inspection, and contact-repair stack, not as a predictive folding
system.

## Next Allowed Research Target

Do not add more broad biology claims before the canonical surface stays clean.

The next research target can be:

```text
real-coordinate contact failure analysis
```

Allowed shape:

```text
improve contact evidence packets
preserve truth-after-prediction scoring
keep global fold-class claims locked
track coverage cost
record uncertainty, failure cohorts, and visible contact-map drift
analyze coordinate-native missed contact cohorts
avoid repair changes until failures are explained on the coordinate benchmark
keep visual_12_is_toy_benchmark = true
keep contact_repair_overfit_risk_reported = true
keep toy_locked_contact_targets_used = false on real-coordinate surfaces
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
