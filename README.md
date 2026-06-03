# Pharmacotopology

A bounded research workbench for exploring how abstract mechanism or protein
signals can perturb an explicit topology model, with safety boundaries kept in
the artifacts instead of left as footnotes.

This project has not solved protein folding. It has built a safe axis-level
falsification stack. The current external-safe result has zero unsafe axis
claims, but coverage is reduced. Global fold-class claims remain locked.

## Current Status

Pharmacotopology now has two related surfaces:

```text
mechanism topology workbench
protein/folding axis falsification stack
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

Current honest state:

```text
internal 50-row safety benchmark = closed as a development safety milestone
external 100-row unsafe axis claims = closed after quarantine
axis_profile_coverage after repair = 0.55
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
```

Canonical commands:

```bash
python3 scripts/run_fold_axis_adjudication_benchmark.py
python3 scripts/run_fold_axis_profile_benchmark.py
python3 scripts/run_architecture_axis_benchmark.py
python3 scripts/run_external_fold_family_holdout_benchmark.py
python3 scripts/run_external_axis_repair_benchmark.py
```

The active generated folding artifacts are:

```text
real_folding_50_axis_*
real_folding_50_axis_profile_*
real_folding_50_architecture_axis_*
external_fold_family_100_*
external_axis_repair_*
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
state.

## Documentation

Start here:

```text
docs/CURRENT_PROJECT_STATUS.md
docs/PROTEIN_FOLDING_TEST_BOUNDARY.md
docs/EXTERNAL_AXIS_REPAIR_BOUNDARY.md
docs/PROTEIN_CENTERED_REFRAME.md
```

The short version:

```text
Safety milestone: yes.
Internal 50-row development benchmark: closed as a safety benchmark.
External unsafe axis claims: closed after quarantine.
Protein folding: not solved.
Next correct work: improve evidence quality only after the canonical surface stays clean.
```

## License

MIT.
