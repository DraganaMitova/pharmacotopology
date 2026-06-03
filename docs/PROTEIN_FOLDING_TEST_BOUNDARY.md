# Protein Folding Test Boundary

For the active stack and current metrics, start with
`docs/CURRENT_PROJECT_STATUS.md`. This file preserves the longer boundary and
development history behind that current stack.

This branch adds a small protein-folding hypothesis lab to the project. The
goal is not to solve folding. The goal is to create a safe benchmark shell
where topology signatures can be compared against folding-relevant reference
summaries.

The test claim is deliberately narrow:

```text
topology signatures extracted from protein sequences and structures
can be compared against folding-relevant benchmark summaries
under explicit uncertainty and failure labels
```

## Modeled Chain

```text
protein sequence
-> predicted topology signature
-> reference topology signature
-> contact-map proxy similarity
-> fold-class match
-> uncertainty radius
-> evidence readiness
```

The signature dimensions are:

```text
sequence_complexity
secondary_structure_balance
contact_map_closure
hydrophobic_core_closure
loop_disorder_pressure
domain_boundary_stability
long_range_contact_order
conformational_flexibility
knot_or_entanglement_signature
uncertainty_radius
```

## Current Benchmark Status

The default benchmark rows are placeholders. They are useful for exercising the
code path, output schema, tests, and safety boundary, but they are not external
validation.

Real calibration requires reference rows derived from external structure
resources such as contact maps, fold-class annotations, structure ensembles, or
curated folding benchmark sets.

The repository now includes a first locked 10-row external benchmark slice:

```text
data/folding_benchmarks_real_10.locked.json
```

It contains real sequences and source accessions from RCSB PDB/CATH-style
labels plus DisProt disorder references:

```text
2 alpha-rich rows
2 beta-rich rows
2 alpha-beta mixed rows
2 multidomain / boundary-sensitive rows
2 disordered / flexible rows
```

Each row includes:

```text
real sequence
source_database
source_accession
reference_structure_source
reference_label_source
is_external_reference = true
curation_notes
```

The predictor receives the sequence only. Labels, accessions, and reference
topology signatures remain benchmark-side scoring metadata.

## External Benchmark Adapter

The benchmark runner can load externally derived rows:

```bash
python3 scripts/run_folding_topology_benchmark.py --benchmark-file data/folding_benchmarks_real.json --require-external
```

The adapter validates that:

```text
the file contains reference rows
the sequence uses supported amino-acid symbols
the reference source looks external, such as pdb:, afdb:, casp:, cath:, scop:, disprot:, or external:
placeholder/example/template source labels are rejected when --require-external is used
all topology signature dimensions are present
all topology signature values are inside 0..1
no drug-design or clinical-use boundary is opened
```

The repo includes:

```text
data/folding_benchmarks_real.example.json
data/folding_benchmarks_real_10.locked.json
data/folding_benchmarks_real_500.locked.json
```

The example file is a schema template, not evidence. Copy it to
`data/folding_benchmarks_real.json` and replace the row with externally derived
structure summaries before using `--require-external`.

The 10-file is a small locked benchmark slice, not a general validation set. On
the current model it reports:

```text
external rows = 10
fold-class matches = 4
fold-class mismatches = 6
accuracy = 0.4
folding_problem_solved = false
```

That result is frozen as:

```text
real_external_label_benchmark_v0
```

It is real in the sense that the rows contain real sequences and external
accessions. It is still label/prototype-based in the topology-reference channel:

```text
sequence -> predicted topology
external label -> locked broad-class prototype
compare
```

It is not yet a structure-derived folding test.

## Structure-Derived Topology Benchmark

The next benchmark channel is:

```text
structure_derived_topology_benchmark
```

It keeps three evidence channels separate:

```text
A. sequence-derived prediction
B. structure-derived topology evidence
C. external human/database label
```

The report therefore shows:

```text
prediction_vs_structure_score
prediction_vs_label_score
structure_vs_label_agreement
```

For folded PDB rows, the structure evidence is extracted from local PDB
coordinate files using a C-alpha contact graph plus HELIX/SHEET/B-factor
signals. For DisProt rows, the evidence channel is explicitly
`disorder_reference`, not coordinates.

The compact checked-in evidence file is:

```text
data/folding_benchmarks_real_10_structure_evidence.json
```

It records:

```text
8 coordinate_contact_graph rows
2 disorder_reference rows
folding_problem_solved = false
```

To rebuild that evidence file from local PDB files:

```bash
python3 scripts/extract_structure_topology_signatures.py \
  --benchmark-file data/folding_benchmarks_real_10.locked.json \
  --pdb-dir /path/to/pdb_files \
  --output data/folding_benchmarks_real_10_structure_evidence.json
```

To run the structure benchmark:

```bash
python3 scripts/run_structure_folding_topology_benchmark.py
```

The current checked-in structure report says:

```text
prediction_vs_structure_accuracy = 0.3
prediction_vs_label_accuracy = 0.4
structure_vs_label_agreement_rate = 0.9
sequence_order_sensitivity_score = 0.0
composition_only_warning = true
revision_required = true
claim_allowed = false
folding_problem_solved = false
```

That is a useful contradiction/revision signal, not a model success claim.

## Sequence-Order Controls

The structure benchmark also runs internal perturbation controls:

```text
real sequence
composition-preserving shuffled sequence
reversed sequence
local-window shuffled sequence
hydrophobic-cluster destroyed sequence
charge-pattern scrambled sequence
```

Generated control sequences are not written to the output files. The output
records only control type, whether composition was preserved, predicted class
change, and signature delta.

The most important metric is:

```text
sequence_order_sensitivity_score
```

The structure benchmark value is `0.0`, meaning the original recipe did not yet
respond to sequence order changes. That is the point of the test.

## Order-Aware Recipe Layer

The next layer is:

```text
order_aware_folding_topology_benchmark
```

It still receives sequence only. It does not receive CATH labels, PDB classes,
DisProt labels, structure-derived signatures, or reference fold classes during
prediction.

It extracts:

```text
hydrophobic cluster topology
charge pattern topology
proline/glycine breaker distribution
cysteine spacing / bridge potential
windowed local structure pressure
segment boundary contrast
long-range closure potential
predicted contact-prior graph
```

Run it with:

```bash
python3 scripts/run_order_aware_folding_topology_benchmark.py
```

It writes:

```text
real_folding_10_order_aware_report.json
real_folding_10_order_aware_rows.csv
real_folding_10_contact_prior.csv
real_folding_10_control_separation.csv
real_folding_10_order_aware_dashboard.html
```

The current checked-in order-aware report says:

```text
sequence_order_sensitivity_score = 0.278079
real_vs_shuffled_separation_mean = 0.278079
contact_prior_signal_seen = true
recipe_order_blind = false
composition_only_warning = false
prediction_vs_structure_accuracy = 0.2
prediction_vs_label_accuracy = 0.1
revision_required = true
claim_allowed = false
folding_problem_solved = false
```

This means the recipe now reacts to sequence order and contact-prior topology,
but it still has not earned a folding claim. Scaling beyond the 10-row slice
should wait until the order-aware signal remains stable and the
prediction-vs-structure behavior improves without threshold tuning.

## Motif-to-Structure Alignment Layer

The next diagnostic layer is:

```text
motif_to_structure_alignment_benchmark
```

It still receives sequence only during prediction. It does not receive CATH
labels, DisProt labels, PDB classes, structure-derived signatures, or reference
fold classes until after the motif evidence vector has already been produced.

It makes the topology evidence vector the first output:

```text
alpha periodicity
beta alternation
compact core
disorder run
domain boundary
long-range closure
breaker / turn
charge frustration
```

Only after that evidence vector exists does the layer map it to a provisional
broad class, compare it against structure/label channels, and emit per-protein
failure diagnosis. The practical output is not higher accuracy yet; it is a
clearer explanation of why the current recipe fails.

Run it with:

```bash
python3 scripts/run_motif_alignment_benchmark.py
```

It writes:

```text
real_folding_10_motif_alignment_report.json
real_folding_10_motif_alignment_rows.csv
real_folding_10_failure_diagnosis.csv
real_folding_10_evidence_conflicts.csv
real_folding_10_motif_alignment_dashboard.html
```

The current checked-in motif alignment report says:

```text
prediction_vs_structure_accuracy = 0.1
prediction_vs_label_accuracy = 0.1
sequence_order_sensitivity_score = 0.278079
real_vs_shuffled_separation_mean = 0.278079
contact_prior_signal_seen = true
motif_signal_seen = true
evidence_conflict_mean = 0.895054
uncertainty_gating_used = true
forced_prediction_count = 2
abstained_prediction_count = 8
high_confidence_wrong_count = 0
revision_required = true
claim_allowed = false
folding_problem_solved = false
```

This means the layer has become more cautious rather than more accurate. That
is intentional for this revision. A high-conflict motif vector now usually
abstains, and high-confidence wrong predictions are gated to zero on this
locked 10-row slice.

## Hierarchical Folding Decision Gates

The next diagnostic layer is:

```text
hierarchical_folding_decision_gate_benchmark
```

It still receives sequence only during prediction. It uses the same motif
evidence vector, but it no longer asks for a fold class first. It routes the
evidence through staged gates:

```text
Gate 1: disorder / flexibility
Gate 2: compactness / closure
Gate 3: domain segmentation
Gate 4: secondary structure
Gate 5: confidence / abstention
```

The important boundary distinction is:

```text
segmented but compact
segmented and disordered
```

Domain-boundary evidence is treated as folded multidomain only when compact
core and long-range closure are also present. If disorder pressure is high and
compact closure is weak, segmentation is interpreted as flexible/disordered
regional structure instead.

Run it with:

```bash
python3 scripts/run_hierarchical_gate_benchmark.py
```

It writes:

```text
real_folding_10_hierarchical_gate_report.json
real_folding_10_hierarchical_gate_rows.csv
real_folding_10_gate_paths.csv
real_folding_10_gate_failures.csv
real_folding_10_hierarchical_gate_dashboard.html
real_folding_50_hierarchical_gate_report.json
real_folding_50_hierarchical_gate_rows.csv
real_folding_50_gate_paths.csv
real_folding_50_gate_failures.csv
real_folding_50_hierarchical_gate_dashboard.html
real_folding_50_certificate.json
real_folding_50_confusion_matrix.csv
```

The current checked-in hierarchical gate report says:

```text
prediction_vs_structure_accuracy = 0.6
prediction_vs_label_accuracy = 0.6
forced_prediction_count = 6
abstained_prediction_count = 4
high_confidence_wrong_count = 0
evidence_conflict_mean = 0.549305
disorder_gate_accuracy = 1.0
compactness_gate_accuracy = 1.0
segmentation_gate_accuracy = 1.0
secondary_structure_gate_accuracy = 0.428571
flexible_segmentation_false_multidomain_count = 0
false_beta_from_disorder_count = 0
false_mixed_from_alpha_count = 0
hierarchy_changed_raw_prediction_count = 9
revision_required = true
claim_allowed = false
folding_problem_solved = false
```

This is still not a folding solution. It is a better routing layer: it reduces
specific gate-order failures while preserving abstention when secondary
structure evidence is not stable enough.

## Locked 50-Row Stability Test

The hierarchy is now frozen against:

```text
data/folding_benchmarks_real_50.locked.json
```

The 50-row file contains externally sourced sequence rows with this label
stratification:

```text
10 alpha_rich
10 beta_rich
10 alpha_beta_mixed
10 multidomain_boundary
10 disordered_flexible
```

The companion structure evidence file is:

```text
data/folding_benchmarks_real_50_structure_evidence.json
```

It contains 40 coordinate-contact graph rows from PDB files and 10
DisProt-style disorder-reference rows. The predictor still receives sequence
only. Structure-derived classes and external labels are used only after
prediction for scoring and failure review.

Run the frozen hierarchy on 50 with:

```bash
python3 scripts/run_hierarchical_gate_benchmark.py \
  --benchmark-file data/folding_benchmarks_real_50.locked.json \
  --structure-evidence-file data/folding_benchmarks_real_50_structure_evidence.json \
  --report-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_hierarchical_gate_report.json \
  --rows-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_hierarchical_gate_rows.csv \
  --gate-paths-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_gate_paths.csv \
  --gate-failures-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_gate_failures.csv \
  --dashboard-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_hierarchical_gate_dashboard.html \
  --certificate-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_certificate.json \
  --confusion-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_confusion_matrix.csv
```

It writes:

```text
real_folding_50_hierarchical_gate_report.json
real_folding_50_hierarchical_gate_rows.csv
real_folding_50_gate_paths.csv
real_folding_50_gate_failures.csv
real_folding_50_hierarchical_gate_dashboard.html
real_folding_50_certificate.json
real_folding_50_confusion_matrix.csv
```

The current checked-in 50-row stability report says:

```text
benchmark_size = 50
external_rows = 50
prediction_vs_structure_accuracy = 0.34
prediction_vs_label_accuracy = 0.38
forced_prediction_count = 25
abstained_prediction_count = 25
high_confidence_wrong_count = 6
false_beta_from_disorder_count = 0
false_mixed_from_alpha_count = 0
flexible_segmentation_false_multidomain_count = 0
disorder_gate_accuracy = 0.9
compactness_gate_accuracy = 0.82
segmentation_gate_accuracy = 0.9
secondary_structure_gate_accuracy = 0.194444
accuracy_delta_from_10 = -0.26
abstention_delta_from_10 = 21
high_confidence_wrong_delta_from_10 = 6
stability_status = unstable_accuracy_drop
claim_allowed = false
folding_problem_solved = false
```

This is a useful negative result. The targeted failure counters stayed closed,
but the hierarchy did not remain stable under the 50-row structure-derived
stress test. The next revision should be based on the 50-row failure
distribution.

## Protein Regime Routing And Failure Cohorts

The next diagnostic layer is:

```text
protein_regime_routing_failure_cohort_analysis
```

This layer does not retune the hierarchy thresholds. It adds a sequence-only
router before accepting the hierarchy output. The router sees only derived
features from the amino-acid sequence:

```text
sequence complexity
hydrophobic fraction and hydrophobic run pressure
disorder-promoting composition
breaker / turn pressure
charge blockiness
alpha / beta / disorder motif pressure
domain-boundary and long-range-closure motif pressure
repeat and periodicity pressure
```

It does not use:

```text
external labels
PDB coordinates
CATH / SCOP classes
DisProt truth labels
structure-derived fold classes
reference topology signatures
```

The predicted protein regimes are:

```text
compact_single_domain
multidomain_modular
intrinsically_disordered
repeat_like
coiled_coil_or_fibrous
membrane_like
low_complexity_region_rich
small_peptide_or_fragment
ambiguous_regime
```

The point is leakage-safe routing. If a sequence looks membrane-like,
repeat-like, fibrous/coiled, too short, or ambiguous, the current broad
fold-class gates are not trusted to force a class. The model abstains and
records the reason. Structure evidence and labels are used only afterwards to
score whether that abstention prevented a bad claim or simply exposed a
manual-review row.

Run it with:

```bash
python3 scripts/run_regime_analysis_benchmark.py
```

It writes:

```text
real_folding_50_regime_analysis_report.json
real_folding_50_regime_rows.csv
real_folding_50_failure_cohorts.csv
real_folding_50_high_confidence_wrong.csv
real_folding_50_abstention_analysis.csv
real_folding_50_regime_dashboard.html
```

The current checked-in regime report says:

```text
benchmark_size = 50
prediction_vs_structure_accuracy = 0.28
prediction_vs_label_accuracy = 0.28
forced_prediction_count = 14
abstained_prediction_count = 36
high_confidence_wrong_count = 0
regime_accuracy = 0.74
regime_confidence_mean = 0.738114
structure_label_disagreement_count = 17
ambiguous_reference_count = 22
possible_bad_rows_count = 8
hierarchical_high_confidence_wrong_count_before_regime_routing = 6
hierarchical_high_confidence_wrong_prevented_by_regime_routing = 6
dominant_failure_cohort = abstained_unresolved
old failure counters = closed
revision_required = true
claim_allowed = false
folding_problem_solved = false
```

This is not an accuracy improvement. It is a clearer failure map. The dominant
cohort says the largest remaining non-correct family is unresolved abstention
requiring better axis evidence, while the old targeted failure modes remain
closed.

## Orthogonal Fold-Axis Truth Adjudication

The next layer is:

```text
orthogonal_fold_axis_truth_adjudication
```

It splits the collapsed single fold class into four axes:

```text
secondary_structure_axis
architecture_axis
order_axis
environment_axis
```

This layer does not retune the regime router and does not make a new
prediction. It adjudicates the existing sequence-only outputs against
structure/label truth axes after prediction. See
`docs/FOLD_AXIS_TRUTH_BOUNDARY.md` for the full boundary.

Run it with:

```bash
python3 scripts/run_fold_axis_adjudication_benchmark.py
```

The current checked-in axis report says:

```text
single_class_taxonomy_collapse_detected = true
structure_label_disagreement_count = 17
orthogonal_axis_disagreement_count = 10
true_same_axis_conflict_count = 8
structure_label_same_axis_conflict_count = 7
forced_same_axis_conflict_count = 0
forced_order_axis_conflict_count = 0
forced_secondary_axis_conflict_count = 0
folded_domain_mimic_abstained_count = 1
secondary_axis_ambiguity_abstained_count = 2
coverage_loss_from_safety_guards = 3
axis_unscorable_count = 147
high_confidence_wrong_count_after_axis_scoring = 0
artifact_reproducible = true
claim_allowed = false
folding_problem_solved = false
```

## Axis-Safe Profile Coverage Recovery

The next layer is:

```text
fold_axis_profile_coverage_recovery
```

It does not recover global fold-class coverage. It emits a partial profile:

```text
secondary_structure_axis
architecture_axis
order_axis
environment_axis
```

Each axis can be known or unknown independently. See
`docs/FOLD_AXIS_PROFILE_BOUNDARY.md` for the full boundary.

Run it with:

```bash
python3 scripts/run_fold_axis_profile_benchmark.py
```

The current checked-in axis profile report says:

```text
collapsed_class_coverage = 0.28
axis_profile_coverage = 0.86
secondary_axis_coverage = 0.16
architecture_axis_coverage = 0.14
order_axis_coverage = 0.84
environment_axis_coverage = 0.04
safe_axis_recovered_count = 39
unsafe_class_recovery_count = 0
guard_override_count = 0
forced_same_axis_conflict_count = 0
axis_profile_same_axis_conflict_count = 0
high_confidence_wrong_count_after_axis_scoring = 0
global_fold_class_claim_allowed = false
axis_profile_claim_allowed = true
folding_problem_solved = false
```

This is safe coverage recovery, not accuracy chasing. It lets the benchmark say
that a row is ordered, membrane-like, fragment-scoped, or weakly secondary-axis
supported without turning an abstention into a global fold class.

## Architecture-Axis Evidence Adjudication

The next architecture-specific layer is:

```text
architecture_axis_evidence_adjudication
```

It creates a sequence-only architecture evidence packet per row and then scores
that packet against architecture truth only after prediction. It separates
architecture from secondary structure, disorder, membrane topology, and
length-only inference. See `docs/FOLD_ARCHITECTURE_AXIS_BOUNDARY.md` for the
full boundary.

Run it with:

```bash
python3 scripts/run_architecture_axis_benchmark.py
```

The current checked-in architecture-axis report says:

```text
previous_profile_architecture_axis_coverage = 0.14
architecture_axis_coverage = 0.32
architecture_axis_claim_allowed_count = 16
architecture_axis_abstained_count = 34
architecture_axis_same_axis_conflict_count = 0
architecture_axis_safe_claim_count = 16
fragment_scope_detected_count = 7
compact_single_domain_claim_count = 6
multidomain_claim_count = 3
repeat_like_claim_count = 0
architecture_claim_without_secondary_leakage_count = 0
architecture_claim_without_label_leakage_count = 0
global_fold_class_claim_allowed = false
axis_profile_claim_allowed = true
claim_allowed = false
folding_problem_solved = false
```

This batch improves architecture coverage while preserving the same safety
closure: no global fold class recovery, no secondary-structure leakage, no
label leakage, and no architecture same-axis conflicts.

## External Fold-Family Holdout

The next layer is:

```text
external_fold_family_holdout_benchmark
```

This is a falsification layer, not a repair layer. It runs the current safe
axis stack without tuning against a locked 100-row non-overlapping holdout. See
`docs/EXTERNAL_FOLD_FAMILY_HOLDOUT_BOUNDARY.md` for the full boundary.

Run it with:

```bash
python3 scripts/run_external_fold_family_holdout_benchmark.py
```

The current checked-in holdout report says:

```text
holdout_row_count = 100
holdout_unique_sequence_count = 100
holdout_unique_family_count = 23
development_overlap_count = 0
development_sequence_overlap_count = 0
holdout_lock_valid = true
holdout_non_overlap_valid = true

axis_profile_coverage = 0.78
architecture_axis_coverage = 0.17
axis_profile_same_axis_conflict_count = 22
architecture_axis_same_axis_conflict_count = 14
forced_same_axis_conflict_count = 22
high_confidence_wrong_count_after_axis_scoring = 17
unsafe_axis_claim_count = 36
unsafe_class_recovery_count = 0
guard_override_count = 0
family_generalization_status = conflicts_detected
global_fold_class_claim_allowed = false
axis_profile_claim_allowed = true
architecture_axis_claim_allowed = false
claim_allowed = false
folding_problem_solved = false
```

This result shows partial generalization and real holdout failures. The next
repair batch should be selected from `external_fold_family_100_failure_cohorts.csv`.

## External-Safe Axis Conflict Quarantine

The first external repair batch is intentionally a safety calibration batch,
not an accuracy batch:

```bash
python3 scripts/run_external_axis_repair_benchmark.py
```

It writes:

```text
real_folding_50_axis_* artifacts regenerated for reproducibility
real_folding_50_axis_profile_* artifacts regenerated for reproducibility
external_axis_repair_report.json
external_axis_repair_rows.csv
external_axis_repair_conflict_delta.csv
external_axis_repair_abstention_delta.csv
external_axis_repair_quarantine_rows.csv
external_axis_repair_family_summary.csv
external_axis_repair_dashboard.html
external_axis_repair_certificate.json
```

Current result:

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
legacy_axis_artifacts_reproducible = true
legacy_axis_profile_artifacts_reproducible = true
global_fold_class_claim_allowed = false
claim_allowed = false
folding_problem_solved = false
```

The repair closes unsafe external axis overclaims by abstention:

```text
weak disordered_flexible order projection -> mixed_or_uncertain
weak repeat-like hydrophobic periodicity -> unknown architecture
```

It does not convert the affected rows into better-looking labels and does not
recover collapsed-class coverage.

The 500-file is a target shell, not a completed benchmark. It records the
intended proof ladder and current lock blockers until real rows exist.

```bash
python3 scripts/build_real_folding_benchmark_500.py \
  --size 500 \
  --output data/folding_benchmarks_real_500.locked.json \
  --lock
```

The intended stratification is:

```text
100 mostly-alpha / compact domains
100 mostly-beta / long-range-contact domains
100 alpha-beta mixed domains
100 multidomain / boundary-sensitive proteins
100 disordered or flexible proteins
```

Run the real-10 slice with:

```bash
python3 scripts/run_folding_topology_benchmark.py \
  --benchmark-file data/folding_benchmarks_real_10.locked.json \
  --require-external \
  --report-output first_contact_clean_pharmacotopology_layer_run/real_folding_10_report.json \
  --csv-output first_contact_clean_pharmacotopology_layer_run/real_folding_10_rows.csv
```

And visualize it with:

```bash
python3 scripts/render_folding_benchmark_dashboard.py \
  --report first_contact_clean_pharmacotopology_layer_run/real_folding_10_report.json \
  --csv first_contact_clean_pharmacotopology_layer_run/real_folding_10_rows.csv \
  --output first_contact_clean_pharmacotopology_layer_run/real_folding_10_dashboard.html
```

When the rows are real and locked, the benchmark can be run with:

```bash
python3 scripts/run_folding_topology_benchmark.py \
  --benchmark-file data/folding_benchmarks_real_500.locked.json \
  --require-external \
  --report-output first_contact_clean_pharmacotopology_layer_run/real_folding_500_report.json \
  --csv-output first_contact_clean_pharmacotopology_layer_run/real_folding_500_rows.csv
```

And visualized with:

```bash
python3 scripts/render_folding_benchmark_dashboard.py \
  --report first_contact_clean_pharmacotopology_layer_run/real_folding_500_report.json \
  --csv first_contact_clean_pharmacotopology_layer_run/real_folding_500_rows.csv \
  --output first_contact_clean_pharmacotopology_layer_run/real_folding_500_dashboard.html
```

## Output Fields

Each benchmark comparison writes:

```text
protein_id
sequence_length
source_database
source_accession
reference_structure_source
reference_label_source
predicted_topology_signature
reference_topology_signature
predicted_internal_fold_class
contact_map_similarity
fold_class_match
uncertainty_radius
evidence_readiness
failure_reason
is_external_reference
reference_topology_signature_kind
curation_notes
```

The default `failure_reason` is expected to say that an external structure
benchmark is not attached. That is a useful failure: it prevents placeholder
numbers from being mistaken for evidence.

External benchmark runs also write sidecar files beside the report:

```text
*_certificate.json
*_failures.csv
*_confusion_matrix.csv
```

Structure-derived benchmark runs write:

```text
real_folding_10_structure_report.json
real_folding_10_structure_rows.csv
real_folding_10_structure_dashboard.html
real_folding_10_order_controls.csv
real_folding_10_falsification_report.json
real_folding_10_order_aware_report.json
real_folding_10_order_aware_rows.csv
real_folding_10_contact_prior.csv
real_folding_10_control_separation.csv
real_folding_10_order_aware_dashboard.html
real_folding_10_motif_alignment_report.json
real_folding_10_motif_alignment_rows.csv
real_folding_10_failure_diagnosis.csv
real_folding_10_evidence_conflicts.csv
real_folding_10_motif_alignment_dashboard.html
real_folding_10_hierarchical_gate_report.json
real_folding_10_hierarchical_gate_rows.csv
real_folding_10_gate_paths.csv
real_folding_10_gate_failures.csv
real_folding_10_hierarchical_gate_dashboard.html
real_folding_50_hierarchical_gate_report.json
real_folding_50_hierarchical_gate_rows.csv
real_folding_50_gate_paths.csv
real_folding_50_gate_failures.csv
real_folding_50_hierarchical_gate_dashboard.html
real_folding_50_certificate.json
real_folding_50_confusion_matrix.csv
real_folding_50_regime_analysis_report.json
real_folding_50_regime_rows.csv
real_folding_50_failure_cohorts.csv
real_folding_50_high_confidence_wrong.csv
real_folding_50_abstention_analysis.csv
real_folding_50_regime_dashboard.html
real_folding_50_axis_adjudication_report.json
real_folding_50_axis_rows.csv
real_folding_50_axis_conflicts.csv
real_folding_50_axis_manual_review.csv
real_folding_50_axis_confusion_matrices.csv
real_folding_50_axis_dashboard.html
real_folding_50_axis_profile_report.json
real_folding_50_axis_profile_rows.csv
real_folding_50_axis_profile_abstentions.csv
real_folding_50_axis_profile_recovery_candidates.csv
real_folding_50_axis_profile_dashboard.html
real_folding_50_axis_profile_certificate.json
real_folding_50_architecture_axis_report.json
real_folding_50_architecture_axis_rows.csv
real_folding_50_architecture_axis_conflicts.csv
real_folding_50_architecture_axis_abstentions.csv
real_folding_50_architecture_axis_dashboard.html
real_folding_50_architecture_axis_certificate.json
external_fold_family_100_report.json
external_fold_family_100_rows.csv
external_fold_family_100_family_summary.csv
external_fold_family_100_axis_conflicts.csv
external_fold_family_100_abstentions.csv
external_fold_family_100_failure_cohorts.csv
external_fold_family_100_dashboard.html
external_fold_family_100_certificate.json
```

When a benchmark file is loaded, the JSON report also includes:

```text
reference_dataset_validation
```

This records how many rows were loaded, how many had external reference-source
labels, whether external rows were required, and any violations or warnings.

The dashboard makes the following proof surfaces visible:

```text
confusion matrix
similarity distribution
per-class accuracy
worst mismatch table
radar overlay
locked benchmark certificate
dataset hash
commit hash
lock blockers
```

## Boundary

This module does not:

```text
solve protein folding
predict atomic structure
generate protein sequences
design molecules
recommend compounds
make clinical claims
```

It is a benchmark adapter and hypothesis-review shell. Its best use is to make
the next evidence step obvious.
