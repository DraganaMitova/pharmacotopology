# Pharmacotopology

A bounded hypothesis workbench and visual-only research simulator for exploring
how hypothesis-level mechanism vectors might perturb an abstract topology
profile.

This repository is public because the idea is easier to examine when the model,
assumptions, generated artifacts, and safety boundaries are visible. It is not
a clinical tool. It does not recommend medication, infer treatment, map brand
names, mention dosage, or provide real-world prescribing guidance.

## What This Is

Pharmacotopology treats a mechanism class as a bounded perturbation over an
abstract topology state:

```text
mechanism vector
-> topology delta
-> pathology reduction score
-> collapse cost score
-> net topology health score
```

The included prototype uses a synthetic schizophrenia-like topology profile.
That profile is a simulated pressure map only; it is not a diagnosis, a patient
model, or a medication target.

The current layer compares hypothesis-only mechanism classes such as
`d2_antagonist_like`, `nmda_support_like`, and
`glutamate_amplification_stressor_like`. Each vector changes one or more
topology dimensions and carries a separate collapse cost. The ranking asks a
narrow modeling question:

```text
Did this perturbation move the simulated topology closer to bounded review,
and at what modeled collapse cost?
```

The practical use is not clinical decision support. It is structured hypothesis
work: make assumptions explicit, compare them under the same scoring rules, add
evidence weights when evidence exists, and keep uncertainty visible until a
claim can be externally calibrated.

## Protein-Centered Reframe

The project is moving from compound-first pharmacology toward protein-centered
mechanism topology. Chemical compounds are treated only as abstract perturbation
sources. The modeled chain is:

```text
abstract_compound_class
-> protein_family
-> protein_state_shift
-> pathway/network perturbation
-> topology_delta
-> collapse_cost
-> evidence_readiness
```

This is still not a drug-design tool. It does not propose molecules, optimize
compounds, recommend treatments, infer patient state, or make clinical claims.
See [docs/PROTEIN_CENTERED_REFRAME.md](docs/PROTEIN_CENTERED_REFRAME.md).

## Protein Folding Hypothesis Lab

The protein-centered layer now has a small folding benchmark shell. It asks
whether abstract topology signatures extracted from protein sequences can be
compared against folding-relevant reference summaries:

```text
protein sequence
-> predicted topology signature
-> reference topology signature
-> contact-map proxy similarity
-> fold-class match
-> uncertainty radius
-> evidence readiness
```

The default references are placeholders. They exercise the schema, metrics, and
safety boundary; they do not validate any folding claim.

The repository now also includes a first locked 10-row external benchmark slice
at `data/folding_benchmarks_real_10.locked.json`. It contains real protein
sequences and external reference accessions from RCSB PDB/CATH-style labels and
DisProt. It is deliberately small: two alpha-rich rows, two beta-rich rows, two
alpha-beta mixed rows, two multidomain or boundary-sensitive rows, and two
disordered/flexible rows. The predictor receives only the sequence; source
labels and topology-reference prototypes are used after prediction for scoring.
On the current model, this slice reports 4 fold-class matches and 6 mismatches
(`accuracy = 0.4`). That is useful because the failures are visible, not because
the model is strong.

That first slice is frozen as `real_external_label_benchmark_v0`: it compares a
sequence-derived prediction against broad external labels and locked
label-derived prototype signatures. The next layer is separate and more
important: `structure_derived_topology_benchmark`.

```text
sequence -> predicted topology
PDB coordinates / DisProt disorder evidence -> extracted topology evidence
external label -> label agreement channel
```

The structure-derived layer keeps three truth channels apart:

```text
prediction_vs_structure_score
prediction_vs_label_score
structure_vs_label_agreement
```

It also runs internal sequence-order controls without writing generated control
sequences: composition shuffle, reversal, local-window shuffle,
hydrophobic-cluster destruction, and charge-pattern scrambling. The current
structure-derived control pass has `sequence_order_sensitivity_score = 0.0`,
which is the honest warning: the original recipe was still reading composition
more than sequence order.

The next checked-in layer is `order_aware_folding_topology_benchmark`. It adds
sequence-position features: hydrophobic cluster topology, charge pattern
topology, proline/glycine breaker distribution, cysteine spacing, windowed
local structure pressure, segment boundary contrast, long-range closure
potential, and a predicted contact-prior graph. On the same 10 locked rows it
now reports `sequence_order_sensitivity_score = 0.278079`,
`real_vs_shuffled_separation_mean = 0.278079`,
`contact_prior_signal_seen = true`, `recipe_order_blind = false`, and
`composition_only_warning = false`. The prediction accuracy is still weak, so
`revision_required = true` and `folding_problem_solved = false`.

The current follow-up layer is `motif_to_structure_alignment_benchmark`. It
makes the topology evidence vector the first output, then only afterwards maps
that evidence to a broad class for failure diagnosis. It exposes alpha
periodicity, beta alternation, compact-core, disorder-run, domain-boundary,
long-range-closure, breaker/turn, and charge-frustration signals. On the same
10 locked rows it reports `forced_prediction_count = 2`,
`abstained_prediction_count = 8`, and `high_confidence_wrong_count = 0`. This
is not better folding accuracy yet. It is a more useful failure surface: the
model now shows when motif evidence is conflicted enough that it should
abstain instead of making a confident class claim.

The next layer is `hierarchical_folding_decision_gate_benchmark`. It routes the
same sequence-only motif evidence through explicit gates: disorder/flexibility,
compactness/closure, domain segmentation, secondary structure, and confidence.
This separates flexible segmentation from folded multidomain evidence and
requires beta evidence to include pairing support, not only alternation. On the
same 10 locked rows it reports `prediction_vs_structure_accuracy = 0.6`,
`forced_prediction_count = 6`, `abstained_prediction_count = 4`,
`high_confidence_wrong_count = 0`, and all three targeted failure counters at
zero: disorder-as-beta, alpha-as-mixed, and flexible-segmentation-as-multidomain.
It is still not a folding solution.

The hierarchy is now frozen against a 50-row stability benchmark at
`data/folding_benchmarks_real_50.locked.json`. That run is intentionally
harsher: `prediction_vs_structure_accuracy = 0.34`,
`abstained_prediction_count = 25`, and `high_confidence_wrong_count = 6`, so
`stability_status = unstable_accuracy_drop`. The useful part is that the
targeted failure counters remain closed (`false_beta_from_disorder_count = 0`,
`false_mixed_from_alpha_count = 0`,
`flexible_segmentation_false_multidomain_count = 0`). The next revision should
come from this 50-row failure distribution, not from retuning the 10-row set.

That follow-up layer is now `protein_regime_routing_failure_cohort_analysis`.
It does not tune the hierarchy. It detects a sequence-only protein regime
first, such as compact single domain, multidomain modular, intrinsically
disordered, membrane-like, small fragment, repeat-like, fibrous/coiled, or
ambiguous. The router can then abstain before accepting a broad fold-class
claim when the regime needs a gate the current hierarchy does not have.

On the locked 50-row set it reports:

```text
prediction_vs_structure_accuracy = 0.28
prediction_vs_label_accuracy = 0.28
forced_prediction_count = 14
abstained_prediction_count = 36
high_confidence_wrong_count = 0
regime_accuracy = 0.74
hierarchical_high_confidence_wrong_prevented_by_regime_routing = 6
dominant_failure_cohort = abstained_unresolved
old failure counters = closed
claim_allowed = false
folding_problem_solved = false
```

This is a diagnosis layer, not an accuracy win. The main finding is that the
largest remaining non-correct cohort is not a reopened gate-order bug; it is a
set of unresolved abstentions that need better axis evidence before any broader
claim is safe.

The next layer is `orthogonal_fold_axis_truth_adjudication`. It splits the
single fold class into secondary-structure, architecture, order, and
environment axes. That makes lossy benchmark rows visible: architecture is not
secondary structure, membrane regime is not fold class, disorder is not simply
absence of alpha/beta, and fragment evidence is not global fold truth. On the
locked 50-row set it detects `structure_label_disagreement_count = 17`,
`orthogonal_axis_disagreement_count = 10`,
`true_same_axis_conflict_count = 8`, and
`forced_same_axis_conflict_count = 0`.

See [docs/PROTEIN_FOLDING_TEST_BOUNDARY.md](docs/PROTEIN_FOLDING_TEST_BOUNDARY.md).
See [docs/FOLD_AXIS_TRUTH_BOUNDARY.md](docs/FOLD_AXIS_TRUTH_BOUNDARY.md).

External benchmark rows can be loaded separately:

```bash
python3 scripts/run_folding_topology_benchmark.py --benchmark-file data/folding_benchmarks_real.json --require-external
```

The repository includes `data/folding_benchmarks_real.example.json` as a schema
template only. It is not evidence, and `--require-external` rejects placeholder
or template source labels.

There is also a 500-protein locked-benchmark target shell at
`data/folding_benchmarks_real_500.locked.json`. It starts empty on purpose. It
records the intended stratification and lock blockers without pretending that
real rows have been collected.

```text
100 mostly-alpha / compact
100 mostly-beta / long-range-contact
100 alpha-beta mixed
100 multidomain / boundary-sensitive
100 disordered or flexible
```

The proof standard remains: same recipe, same code, locked external rows,
visible failures, then blind repeatability.

## What This Is Not

This project intentionally refuses clinical interpretation:

```text
not medication advice
not treatment ranking
not symptom prediction
not side-effect prediction
not brand-name mapping
not patient-specific inference
```

The useful reading is:

```text
perturbation changes topology dimensions
```

The unsafe reading is:

```text
pill works / pill does not work
```

## Why It Exists

The experiment is modest: make the hidden assumptions of a topology-style
pharmacology metaphor inspectable. A mechanism vector should not get to look
"helpful" just because it suppresses pressure somewhere; it also has to carry
its modeled collapse cost beside the reduction score.

That means the layer can distinguish between:

```text
reduced unstable topology with low collapse cost
reduced unstable topology with high collapse cost
weak topology correction
destabilizing perturbation
```

The numbers are hypotheses, not truth. Future calibration would need external
evidence such as trial outcomes, receptor binding data, side-effect profiles,
symptom scales, relapse data, cognition data, and discontinuation rates. Until
then, this remains a bounded perturbation model.

The repository now exposes that gap directly through calibration readiness
metadata: evidence weight, uncertainty radius, score intervals, blockers, and
next calibration steps. That makes the project useful as research
infrastructure without pretending to be useful as medicine.

## Repository Shape

```text
src/pharmacotopology/                 core simulator and boundary model
data/folding_benchmarks_real.example.json
data/folding_benchmarks_real_10.locked.json
data/folding_benchmarks_real_10_structure_evidence.json
data/folding_benchmarks_real_50.locked.json
data/folding_benchmarks_real_50_structure_evidence.json
data/folding_benchmarks_real_500.locked.json
data/folding_benchmarks_external_fold_family_100.locked.json
scripts/run_clean_pharmacotopology_layer.py
scripts/render_pharmacotopology_dashboard.py
scripts/render_profile_comparison_dashboard.py
scripts/export_pharmacotopology_csv.py
scripts/run_sensitivity_analysis.py
scripts/explore_sensitivity.py
scripts/build_real_folding_benchmark_500.py
scripts/run_folding_topology_benchmark.py
scripts/render_folding_benchmark_dashboard.py
scripts/extract_structure_topology_signatures.py
scripts/run_structure_folding_topology_benchmark.py
scripts/run_order_aware_folding_topology_benchmark.py
scripts/run_motif_alignment_benchmark.py
scripts/run_hierarchical_gate_benchmark.py
scripts/run_regime_analysis_benchmark.py
scripts/run_fold_axis_adjudication_benchmark.py
scripts/run_fold_axis_profile_benchmark.py
scripts/run_architecture_axis_benchmark.py
scripts/run_external_fold_family_holdout_benchmark.py
scripts/run_external_axis_repair_benchmark.py
scripts/validate_field_trace.py
scripts/measure_field_trace.py
tests/                               pytest coverage for ranking and safety
docs/ADDING_PROFILES.md              small guide for adding synthetic profiles
docs/PROTEIN_CENTERED_REFRAME.md     protein-centered mechanism note
docs/PROTEIN_FOLDING_TEST_BOUNDARY.md folding benchmark boundary note
docs/FOLD_AXIS_TRUTH_BOUNDARY.md     orthogonal fold-axis truth boundary
docs/FOLD_AXIS_PROFILE_BOUNDARY.md   axis-safe partial profile boundary
docs/FOLD_ARCHITECTURE_AXIS_BOUNDARY.md architecture-axis evidence boundary
docs/EXTERNAL_FOLD_FAMILY_HOLDOUT_BOUNDARY.md external holdout boundary
docs/EXTERNAL_AXIS_REPAIR_BOUNDARY.md external-safe axis repair boundary
PHARMACOTOPOLOGY_LAYER.md             longer design note
LICENSE                              MIT license
```

The package has no runtime dependencies and targets Python 3.9+.

## Run

Generate a clean run:

```bash
python3 scripts/run_clean_pharmacotopology_layer.py
```

Run against another synthetic topology profile:

```bash
python3 scripts/run_clean_pharmacotopology_layer.py --profile anxiety_like
```

Render the static dashboard:

```bash
python3 scripts/render_pharmacotopology_dashboard.py
```

Render the multi-profile dashboard:

```bash
python3 scripts/render_profile_comparison_dashboard.py
```

Run the protein folding topology benchmark shell:

```bash
python3 scripts/run_folding_topology_benchmark.py
```

Run it against an externally derived benchmark file:

```bash
python3 scripts/run_folding_topology_benchmark.py --benchmark-file data/folding_benchmarks_real.json --require-external
```

Run the first locked 10-row external benchmark slice and render its dashboard:

```bash
python3 scripts/run_folding_topology_benchmark.py --benchmark-file data/folding_benchmarks_real_10.locked.json --require-external --report-output first_contact_clean_pharmacotopology_layer_run/real_folding_10_report.json --csv-output first_contact_clean_pharmacotopology_layer_run/real_folding_10_rows.csv
python3 scripts/render_folding_benchmark_dashboard.py --report first_contact_clean_pharmacotopology_layer_run/real_folding_10_report.json --csv first_contact_clean_pharmacotopology_layer_run/real_folding_10_rows.csv --output first_contact_clean_pharmacotopology_layer_run/real_folding_10_dashboard.html
```

That run also writes:

```text
real_folding_10_certificate.json
real_folding_10_failures.csv
real_folding_10_confusion_matrix.csv
```

Run the structure-derived topology benchmark and order controls:

```bash
python3 scripts/run_structure_folding_topology_benchmark.py
```

If you want to rebuild the compact structure evidence file from local PDB
coordinate files, place files such as `1MBN.pdb` and `2LZM.pdb` in a directory
and run:

```bash
python3 scripts/extract_structure_topology_signatures.py --pdb-dir /path/to/pdb_files
```

The checked-in structure benchmark report currently says:

```text
coordinate rows = 8
disorder-reference rows = 2
prediction_vs_structure_accuracy = 0.3
prediction_vs_label_accuracy = 0.4
structure_vs_label_agreement_rate = 0.9
sequence_order_sensitivity_score = 0.0
revision_required = true
folding_problem_solved = false
```

Run the order-aware recipe and its contact-prior/control-separation package:

```bash
python3 scripts/run_order_aware_folding_topology_benchmark.py
```

The checked-in order-aware report currently says:

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

Run the motif-to-structure alignment layer and failure diagnosis package:

```bash
python3 scripts/run_motif_alignment_benchmark.py
```

The checked-in motif alignment report currently says:

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

Run the hierarchical folding decision gates:

```bash
python3 scripts/run_hierarchical_gate_benchmark.py
```

The checked-in hierarchical gate report currently says:

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

Run the frozen hierarchy against the locked 50-row stability benchmark:

```bash
python3 scripts/run_hierarchical_gate_benchmark.py --benchmark-file data/folding_benchmarks_real_50.locked.json --structure-evidence-file data/folding_benchmarks_real_50_structure_evidence.json --report-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_hierarchical_gate_report.json --rows-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_hierarchical_gate_rows.csv --gate-paths-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_gate_paths.csv --gate-failures-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_gate_failures.csv --dashboard-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_hierarchical_gate_dashboard.html --certificate-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_certificate.json --confusion-output first_contact_clean_pharmacotopology_layer_run/real_folding_50_confusion_matrix.csv
```

The checked-in 50-row stability report currently says:

```text
benchmark_size = 50
external_rows = 50
class_distribution = 10 per broad class
prediction_vs_structure_accuracy = 0.34
prediction_vs_label_accuracy = 0.38
forced_prediction_count = 25
abstained_prediction_count = 25
high_confidence_wrong_count = 6
false_beta_from_disorder_count = 0
false_mixed_from_alpha_count = 0
flexible_segmentation_false_multidomain_count = 0
accuracy_delta_from_10 = -0.26
abstention_delta_from_10 = 21
high_confidence_wrong_delta_from_10 = 6
stability_status = unstable_accuracy_drop
claim_allowed = false
folding_problem_solved = false
```

Run the sequence-only protein-regime routing and failure cohort analysis:

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

The checked-in regime analysis report currently says:

```text
prediction_vs_structure_accuracy = 0.28
prediction_vs_label_accuracy = 0.28
forced_prediction_count = 14
abstained_prediction_count = 36
high_confidence_wrong_count = 0
regime_accuracy = 0.74
structure_label_disagreement_count = 17
possible_bad_rows_count = 8
dominant_failure_cohort = abstained_unresolved
hierarchical_high_confidence_wrong_prevented_by_regime_routing = 6
claim_allowed = false
folding_problem_solved = false
```

Run the orthogonal fold-axis truth adjudication layer:

```bash
python3 scripts/run_fold_axis_adjudication_benchmark.py
```

It writes:

```text
real_folding_50_axis_adjudication_report.json
real_folding_50_axis_rows.csv
real_folding_50_axis_conflicts.csv
real_folding_50_axis_manual_review.csv
real_folding_50_axis_confusion_matrices.csv
real_folding_50_axis_dashboard.html
```

The checked-in axis adjudication report currently says:

```text
single_class_taxonomy_collapse_detected = true
structure_label_disagreement_count = 17
orthogonal_axis_disagreement_count = 10
true_same_axis_conflict_count = 8
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

Run the axis-safe partial profile coverage layer:

```bash
python3 scripts/run_fold_axis_profile_benchmark.py
```

It writes:

```text
real_folding_50_axis_profile_report.json
real_folding_50_axis_profile_rows.csv
real_folding_50_axis_profile_abstentions.csv
real_folding_50_axis_profile_recovery_candidates.csv
real_folding_50_axis_profile_dashboard.html
real_folding_50_axis_profile_certificate.json
```

The checked-in axis profile report currently says:

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

Run the architecture-axis evidence adjudication layer:

```bash
python3 scripts/run_architecture_axis_benchmark.py
```

It writes:

```text
real_folding_50_architecture_axis_report.json
real_folding_50_architecture_axis_rows.csv
real_folding_50_architecture_axis_conflicts.csv
real_folding_50_architecture_axis_abstentions.csv
real_folding_50_architecture_axis_dashboard.html
real_folding_50_architecture_axis_certificate.json
```

The checked-in architecture-axis report currently says:

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

Run the external fold-family holdout falsification benchmark:

```bash
python3 scripts/run_external_fold_family_holdout_benchmark.py
```

It writes:

```text
external_fold_family_100_report.json
external_fold_family_100_rows.csv
external_fold_family_100_family_summary.csv
external_fold_family_100_axis_conflicts.csv
external_fold_family_100_abstentions.csv
external_fold_family_100_failure_cohorts.csv
external_fold_family_100_dashboard.html
external_fold_family_100_certificate.json
```

The checked-in holdout report currently says:

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

This is a falsification result, not a regression fix. The current stack was run
without threshold tuning on non-overlapping fold-family rows, and the failure
cohorts now define the next repair target.

Run the external-safe axis repair benchmark:

```bash
python3 scripts/run_external_axis_repair_benchmark.py
```

It writes:

```text
external_axis_repair_report.json
external_axis_repair_rows.csv
external_axis_repair_conflict_delta.csv
external_axis_repair_abstention_delta.csv
external_axis_repair_quarantine_rows.csv
external_axis_repair_family_summary.csv
external_axis_repair_dashboard.html
external_axis_repair_certificate.json
```

The checked-in repair report currently says:

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
global_fold_class_claim_allowed = false
claim_allowed = false
folding_problem_solved = false
```

This is a safety repair, not an accuracy claim. It refuses weak
`disordered_flexible -> order_axis` projection and refuses `repeat_like`
architecture when the signal is mostly hydrophobic periodicity without enough
recurrence support. The repair lowers coverage by abstention instead of forcing
new labels.

Build the 500-protein locked benchmark target shell:

```bash
python3 scripts/build_real_folding_benchmark_500.py --size 500 --output data/folding_benchmarks_real_500.locked.json --lock
```

After externally sourced rows are attached, run and render the proof dashboard:

```bash
python3 scripts/run_folding_topology_benchmark.py --benchmark-file data/folding_benchmarks_real_500.locked.json --require-external --report-output first_contact_clean_pharmacotopology_layer_run/real_folding_500_report.json --csv-output first_contact_clean_pharmacotopology_layer_run/real_folding_500_rows.csv
python3 scripts/render_folding_benchmark_dashboard.py --report first_contact_clean_pharmacotopology_layer_run/real_folding_500_report.json --csv first_contact_clean_pharmacotopology_layer_run/real_folding_500_rows.csv --output first_contact_clean_pharmacotopology_layer_run/real_folding_500_dashboard.html
```

Default artifacts are written to:

```text
first_contact_clean_pharmacotopology_layer_run/
```

The dashboard is generated at:

```text
first_contact_clean_pharmacotopology_layer_run/pharmacotopology_dashboard.html
```

The multi-profile dashboard is generated at:

```text
first_contact_clean_pharmacotopology_layer_run/multi_profile_dashboard.html
```

It is plain HTML/CSS/SVG and does not require a web server or external packages.

## Profiles

The built-in source profiles are synthetic pressure maps:

```text
schizophrenia_like
depression_like
mania_like
anxiety_like
mixed_state_like
```

They are for hypothesis comparison only. They are not diagnoses, patient
models, or medication targets.

To add another bounded synthetic profile, see [docs/ADDING_PROFILES.md](docs/ADDING_PROFILES.md).

## Export CSV

Export machine-readable rankings and per-dimension deltas:

```bash
python3 scripts/export_pharmacotopology_csv.py
```

Default CSV outputs:

```text
first_contact_clean_pharmacotopology_layer_run/pharmacotopology_rankings.csv
first_contact_clean_pharmacotopology_layer_run/pharmacotopology_deltas.csv
```

The ranking CSV includes score intervals, evidence weight, uncertainty radius,
and evidence readiness labels. The delta CSV includes every topology dimension
plus collapse cost.

## Calibration Readiness Table

The dashboard includes a per-vector calibration table with:

```text
mechanism_id
abstract compound class
protein family
protein mechanism class
protein state shift
pathway/network perturbation
evidence_stage
evidence_weight
uncertainty_radius
evidence_readiness_label
primary evidence sources
evidence refs
calibration blockers
interval kind
pathology_reduction interval
collapse_cost interval
net_topology_health interval
```

The current default vectors remain `uncalibrated_hypothesis`, which is the
honest starting point. The table is meant to show exactly where evidence would
attach later.

The dashboard also includes a synthetic profile comparison table that runs the
same mechanism set against every built-in profile and reports top mechanism,
fit label, net interval, and destabilizing count. It is a ranking robustness
view, not a clinical comparison.

## Multi-Profile Dashboard

Render a dedicated cross-profile dashboard:

```bash
python3 scripts/render_profile_comparison_dashboard.py
```

It includes:

```text
baseline topology maps for every built-in profile
combined ranking table with rank / net score per profile
mechanism/profile net-score heatmap
consistently positive vs highly profile-specific summary
```

## Sensitivity Analysis

Run a local pressure sweep to see whether rankings are stable when source
profile pressures move up or down:

```bash
python3 scripts/run_sensitivity_analysis.py
```

Sweep a smaller set of dimensions:

```bash
python3 scripts/run_sensitivity_analysis.py --profile anxiety_like --dimensions threat_propagation,sleep_instability
```

Default sensitivity outputs:

```text
first_contact_clean_pharmacotopology_layer_run/sensitivity_analysis_report.json
first_contact_clean_pharmacotopology_layer_run/sensitivity_rankings.csv
```

This is a robustness review, not a clinical interpretation.

## Sensitivity Explorer

Explore a mechanism assumption with deterministic Monte Carlo perturbations:

```bash
python3 scripts/explore_sensitivity.py --profile schizophrenia_like --mechanism nmda_support_like --vary collapse_cost 0.05:0.25
```

You can also vary a delta dimension:

```bash
python3 scripts/explore_sensitivity.py --mechanism nmda_support_like --vary delta:cognitive_fragmentation -0.30:-0.05
```

Default explorer outputs:

```text
first_contact_clean_pharmacotopology_layer_run/sensitivity_explorer_report.json
first_contact_clean_pharmacotopology_layer_run/sensitivity_explorer_samples.csv
```

The explorer reports a net-score distribution for the selected mechanism and a
robustness score for every mechanism under sampled perturbations.

## Interactive Wrapper

A small Streamlit or Gradio wrapper would be a natural later layer for local
exploration: profile selection, pressure sliders, ranking preview, CSV export,
and the same safety warnings. It is intentionally not included in the core
package yet, so the current project remains dependency-free and easy to audit.

## Validate

Run the tests:

```bash
python3 -m pytest
```

Validate and measure a generated trace:

```bash
python3 scripts/validate_field_trace.py first_contact_clean_pharmacotopology_layer_run
python3 scripts/measure_field_trace.py first_contact_clean_pharmacotopology_layer_run
```

Expected safety posture:

```text
coherence = 1.0
stop_integrity = 1.0
surface_leakage = 0.0
```

The run report also asserts:

```text
simulation_only = true
hypothesis_numbers_only = true
calibration_status = uncalibrated_hypothesis_workbench
practical_use = bounded_hypothesis_comparison_and_falsification
clinical_use_allowed = false
clinical_advice_created = false
medication_recommendation_created = false
real_patient_inference_created = false
brand_name_mapping_created = false
treatment_claim_created = false
voice_opened = false
schedule_created = false
command_queue_created = false
autonomous_loop_created = false
```

## Output Artifacts

The local tools can create a small set of auditable artifacts:

```text
input_stream.jsonl
memory.jsonl
audit.jsonl
session_records.jsonl
clean_pharmacotopology_layer_report.json
calibration_readiness_report.json
pharmacotopology_rankings.csv
pharmacotopology_deltas.csv
sensitivity_analysis_report.json
sensitivity_rankings.csv
pharmacotopology_dashboard.html
multi_profile_dashboard.html
sensitivity_explorer_report.json
sensitivity_explorer_samples.csv
folding_topology_benchmark_report.json
folding_topology_benchmark.csv
real_folding_500_report.json
real_folding_500_rows.csv
real_folding_500_dashboard.html
real_folding_500_certificate.json
real_folding_500_failures.csv
real_folding_500_confusion_matrix.csv
real_folding_10_report.json
real_folding_10_rows.csv
real_folding_10_dashboard.html
real_folding_10_certificate.json
real_folding_10_failures.csv
real_folding_10_confusion_matrix.csv
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
field_validation.json
field_metrics.json
```

The JSON files are meant to make the review trace auditable. The dashboard is
meant to make the topology ranking, evidence readiness, and uncertainty
intervals easier to inspect without changing the safety boundary.

## Safety Footer

```text
Φ.review does not recommend medication.
Φ.review does not infer treatment.
Φ.review only simulates bounded topology perturbations.
Clinical claims are denied unless externally validated.
```

## Status

This is an early research prototype with a practical non-clinical target:
bounded hypothesis comparison and falsification. The most valuable
contributions right now are careful critique of the modeling assumptions,
evidence-source mapping, calibration datasets, clearer validation boundaries,
and better ways to explain what the scores do and do not mean.

## License

MIT. See [LICENSE](LICENSE).
