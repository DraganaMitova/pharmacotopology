# Pharmacotopology

A bounded research workbench for exploring how abstract mechanism or protein
signals can perturb an explicit topology model, with safety boundaries kept in
the artifacts instead of left as footnotes.

This project has not solved protein folding. It has built a safe axis-level
falsification stack, a visual contact-map workbench, a first sequence-only
contact-repair pass, and now an audit that freezes the 12-row visual benchmark
as toy/coarse/internal. The current external-safe result has zero unsafe axis
claims, but coverage is reduced. Global fold-class and mechanism-discovery
claims remain locked.

## Current Status

Pharmacotopology now has two related surfaces:

```text
mechanism topology workbench
protein/folding axis falsification stack
visual folding mechanism workbench
contact-topology repair workbench
visual mechanism claim audit
```

The mechanism workbench is a hypothesis simulator. It compares abstract
mechanism vectors against synthetic topology profiles and exposes uncertainty,
collapse cost, and calibration blockers.

The protein/folding stack is a sequence-only safety benchmark. It does not try
to produce a final fold class. It asks which partial axes can be claimed safely:

```text
secondary_structure_axis
architecture_axis
order_axis
environment_axis
```

The visual folding mechanism workbench is a locked 12-row contact-map
visualization layer. It renders native contact targets, sequence-only predicted
contact candidates, overlays, coarse trajectories, closure curves, and failure
cohorts. It is interesting because the failures are visible instead of hidden
inside one score.

The contact-topology repair workbench explains the visual failures with
native-gap clusters after prediction, then tests sequence-only long-range and
beta-pairing repairs. It improves the locked visual benchmark without unlocking
global fold claims.

The visual mechanism audit is deliberately conservative. It says the 12-row
visual/contact-repair benchmark is a toy internal benchmark with coarse native
contact targets, and it reports overfit risk from hardcoded beta-registry
patterns.

Current honest state:

```text
internal 50-row safety benchmark = closed as a development safety milestone
external 100-row unsafe axis claims = closed after quarantine
axis_profile_coverage after repair = 0.55
visual mechanism rows rendered = 12/12
visual mechanism visible partial successes = 5
contact repair visible partial successes = 8
contact repair visible failures = 4
visual_12_is_toy_benchmark = true
contact_repair_overfit_risk_reported = true
mechanism_discovery_claim_allowed = false
global_fold_class_claim_allowed = false
folding_problem_solved = false
```

## Canonical Folding Stack

The active folding stack is:

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
```

Canonical commands:

```bash
python3 scripts/run_fold_axis_adjudication_benchmark.py
python3 scripts/run_fold_axis_profile_benchmark.py
python3 scripts/run_architecture_axis_benchmark.py
python3 scripts/run_external_fold_family_holdout_benchmark.py
python3 scripts/run_external_axis_repair_benchmark.py
python3 scripts/run_visual_folding_mechanism_benchmark.py
python3 scripts/run_contact_topology_repair_benchmark.py
python3 scripts/run_visual_mechanism_audit.py
```

The active generated folding artifacts are:

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
```

Older 10-row and pre-axis 50-row artifacts are preserved under:

```text
first_contact_clean_pharmacotopology_layer_run/archived_legacy/
```

They are development history, not the current benchmark surface.

## Current Metrics

Internal 50-row axis profile:

```text
axis_profile_coverage = 0.86
axis_profile_same_axis_conflict_count = 0
high_confidence_wrong_count_after_axis_scoring = 0
forced_same_axis_conflict_count = 0
global_fold_class_claim_allowed = false
folding_problem_solved = false
```

External 100-row holdout before repair:

```text
axis_profile_coverage = 0.78
axis_profile_same_axis_conflict_count = 22
architecture_axis_same_axis_conflict_count = 14
unsafe_axis_claim_count = 36
high_confidence_wrong_count_after_axis_scoring = 17
```

External 100-row holdout after repair:

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

That coverage loss is intentional. The repair chooses abstention over unsafe
axis claims.

Visual mechanism 12-row workbench:

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

This does not make the model stronger by assertion. It makes the mechanism
more inspectable.

Contact-topology repair 12-row workbench:

```text
baseline_visible_partial_success_count = 5
visible_partial_success_count = 8
visible_failure_count = 4
repaired_mean_contact_map_f1 = 0.263729
long_range_contact_recall_delta = 0.472223
beta_pairing_contact_recall_delta = 0.6541
native_truth_used_before_prediction = false
native_truth_used_before_repair = false
raw_sequence_exposed = false
global_folding_claim_allowed = false
folding_problem_solved = false
```

This repair is still coarse. It explains and improves contact topology on the
locked visual benchmark; it does not infer a real folded structure.

Visual mechanism audit:

```text
visual_12_is_toy_benchmark = true
coarse_native_contacts_only = true
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

## Safety Boundaries

This repository does not:

```text
recommend medication
rank treatments
infer patient state
make clinical claims
design molecules, binders, proteins, drugs, therapies, or dosing logic
claim protein folding is solved
recover global fold-class coverage
export raw sequences in generated benchmark dashboards/CSVs
use truth labels before prediction
```

Safety outputs are explicit:

```text
claim_allowed = false
clinical_use_allowed = false
global_fold_class_claim_allowed = false
folding_problem_solved = false
```

## Repository Shape

```text
src/pharmacotopology/                 simulator, feature extraction, and safety logic
data/                                locked benchmark data and schema examples
scripts/                             canonical runners and historical dev runners
docs/CURRENT_PROJECT_STATUS.md        current truth of the project
docs/PROTEIN_FOLDING_TEST_BOUNDARY.md folding benchmark boundary note
docs/FOLD_AXIS_TRUTH_BOUNDARY.md      orthogonal axis truth boundary
docs/FOLD_AXIS_PROFILE_BOUNDARY.md    axis-safe profile boundary
docs/FOLD_ARCHITECTURE_AXIS_BOUNDARY.md architecture-axis boundary
docs/EXTERNAL_FOLD_FAMILY_HOLDOUT_BOUNDARY.md external holdout boundary
docs/EXTERNAL_AXIS_REPAIR_BOUNDARY.md external-safe repair boundary
docs/VISUAL_FOLDING_MECHANISM_BOUNDARY.md visual mechanism boundary
docs/VISUAL_CONTACT_REPAIR_BOUNDARY.md contact repair boundary
docs/VISUAL_MECHANISM_AUDIT.md visual claim audit
tests/                               pytest coverage for safety and reproducibility
```

Low-level feature modules such as `folding_topology.py`,
`folding_order_aware.py`, `folding_motif_alignment.py`, and
`folding_structure_benchmark.py` are still used as foundations. They are not
part of the public canonical stack, but they are not deleted because the active
pipeline still depends on their primitives.

## Run

Run the mechanism simulator:

```bash
python3 scripts/run_clean_pharmacotopology_layer.py
```

Run the canonical folding stack:

```bash
python3 scripts/run_fold_axis_adjudication_benchmark.py
python3 scripts/run_fold_axis_profile_benchmark.py
python3 scripts/run_architecture_axis_benchmark.py
python3 scripts/run_external_fold_family_holdout_benchmark.py
python3 scripts/run_external_axis_repair_benchmark.py
python3 scripts/run_visual_folding_mechanism_benchmark.py
python3 scripts/run_contact_topology_repair_benchmark.py
python3 scripts/run_visual_mechanism_audit.py
```

Run verification:

```bash
python3 -m pytest
python3 -m compileall src scripts tests
python3 scripts/validate_field_trace.py first_contact_clean_pharmacotopology_layer_run
python3 scripts/measure_field_trace.py first_contact_clean_pharmacotopology_layer_run
```

Export a clean tracked-only zip:

```bash
git archive --format=zip --output pharmacotopology-clean.zip HEAD
```

This avoids packaging `__pycache__/`, `.pytest_cache/`, and other local runtime
state. Finder zips are not accepted as reproducibility evidence.

## Documentation

Start here:

```text
docs/CURRENT_PROJECT_STATUS.md
docs/PROTEIN_FOLDING_TEST_BOUNDARY.md
docs/EXTERNAL_AXIS_REPAIR_BOUNDARY.md
docs/VISUAL_FOLDING_MECHANISM_BOUNDARY.md
docs/VISUAL_CONTACT_REPAIR_BOUNDARY.md
docs/VISUAL_MECHANISM_AUDIT.md
docs/PROTEIN_CENTERED_REFRAME.md
```

The short version:

```text
Safety milestone: yes.
Internal 50-row development benchmark: closed as a safety benchmark.
External unsafe axis claims: closed after quarantine.
Visual mechanism workbench: contact-map evidence is visible for 12 locked rows.
Contact repair workbench: visible partial successes improved from 5 to 8.
Visual mechanism audit: toy benchmark and overfit risk are explicitly reported.
Protein folding: not solved.
Next correct work: improve visual/mechanism evidence quality without unlocking global claims.
```

## License

MIT.
