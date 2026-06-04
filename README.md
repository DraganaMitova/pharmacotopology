# Pharmacotopology

A bounded research workbench for exploring how abstract mechanism or protein
signals can perturb an explicit topology model, with safety boundaries kept in
the artifacts instead of left as footnotes.

This project has not solved protein folding. It has built a safe axis-level
falsification stack, a visual contact-map workbench, a first sequence-only
contact-repair pass, an audit that freezes the 12-row visual benchmark as
toy/coarse/internal, and a new real-coordinate visual contact benchmark. The
real-coordinate benchmark is now followed by a contact-law threshold search
that rejects the current scalar contact heuristic as a stable law. The current
threshold search is followed by a cooperative folding-nucleus closure search.
That nucleus search recovers long-range native regions better than pair-level
thresholds, but it overgenerates false nuclei. Global fold-class and
mechanism-discovery claims remain locked.
The active physical selection layer now reranks and gates nuclei with coarse
physical-state terms, then rejects the selector law after ablation.
The coupling-preserved nucleus layer adds a locked coupling-constraint channel
and an iterative trace loop. It passes selector survival gates under an
oracle-control coupling file derived from coordinate-native contacts, so it
identifies a missing channel without unlocking mechanism claims.

## Current Status

Pharmacotopology now has several related surfaces:

```text
mechanism topology workbench
protein/folding axis falsification stack
visual folding mechanism workbench
contact-topology repair workbench
visual mechanism claim audit
real-coordinate visual contact benchmark
contact-law threshold falsification
folding-nucleus closure search
competitive nucleus selection
nucleus graph selectivity and decoy falsification
active physical selection and term ablation
coupling-preserved nucleus selector
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

The real-coordinate visual contact benchmark upgrades the proof target. It
stores locked C-alpha coordinate traces for eight RCSB PDB chains and derives
native contact maps from coordinate geometry at runtime. Prediction remains
sequence-only and blind to coordinates until scoring. This is stronger than the
toy locked-contact target, but still only a coarse C-alpha contact benchmark.

The contact-law threshold search tests whether one sequence-only threshold law
survives across held-out coordinate-backed proteins. The current scalar contact
score is rejected. A pair-plus-entropy candidate improves held-out F1, but it
does not survive as a general law because long-range recall collapses.

The folding-nucleus closure search tests cooperative segment closure events
instead of independent residue-pair contacts. It recovers long-range
coordinate-native regions better than the pair-level threshold baseline, but
the false-nucleus rate remains too high for any law claim.

The competitive nucleus selection layer reduces the closure-event flood with a
fixed sequence-only budget and frustration/geometry/overlap audit. It preserves
some long-range recovery, but it still fails false-rate and cluster-precision
survival gates.

The nucleus graph selectivity layer asks whether selected closure-graph cores
beat matched decoys. They do not: matched decoys currently overlap native
contacts more often than the selected graph, so the graph law is rejected.

The physical closure-state evaluator instantiates those graph-selected closures
as coarse sequence-only physical states. Physical scores beat matched decoys by
score, but false nuclei and contact precision still fail, so the physical-state
law is rejected too.

The active physical selection layer moves those coarse physical terms from
audit into selection. Physical reranking improves decoy enrichment and
long-range recall, while viability gates reduce false nuclei and improve
precision. No selector passes all survival gates, and ablation says only
burial gain clearly helps under the current proxy.

The coupling-preserved nucleus selector tests the missing native-selective
constraint channel. Coupling rerank improves false-nucleus rate, precision,
long-range recall, and decoy enrichment. The iterative trace loop selects a
small compatible closure path with zero false nuclei and passes survival gates.
The checked-in coupling file is explicitly an oracle control, not external
evolutionary evidence, so mechanism-discovery claims stay locked.

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
real_coordinate_visual_rows = 8
real_coordinate_mean_contact_map_f1 = 0.051198
real_coordinate_visible_partial_success_count = 3
real_coordinate_visible_failure_count = 5
toy_locked_contact_targets_used = false
coarse_ca_only = true
current_scalar_score_law_rejected = true
best_law_candidate_model = pair_plus_entropy_score
best_law_candidate_loo_mean_test_f1 = 0.233569
law_generalizes = false
nucleus_native_recall = 0.50826
nucleus_long_range_recall = 0.913242
nucleus_false_rate = 0.692721
nucleus_law_survives = false
competitive_selected_events = 686
competitive_event_reduction_ratio = 0.821633
competitive_long_range_recall = 0.588907
competitive_false_nucleus_rate = 0.594592
competitive_contact_cluster_precision = 0.044608
competitive_nucleus_law_survives = false
nucleus_graph_selected_events = 320
nucleus_graph_false_nucleus_rate = 0.609375
nucleus_graph_contact_cluster_precision = 0.043311
nucleus_graph_long_range_recall = 0.372496
nucleus_graph_real_vs_decoy_enrichment_ratio = 0.806452
nucleus_graph_law_survives = false
physical_state_count = 320
physical_real_vs_decoy_enrichment_ratio = 1.171535
physical_false_nucleus_rate = 0.609375
physical_contact_cluster_precision = 0.043311
physical_long_range_contact_recall = 0.372496
physical_state_law_survives = false
physical_rerank_false_nucleus_rate = 0.559375
physical_rerank_cluster_precision = 0.050488
physical_rerank_long_range_recall = 0.405008
physical_rerank_real_vs_decoy_enrichment_ratio = 1.585882
physical_gate_false_nucleus_rate = 0.316338
physical_gate_cluster_precision = 0.069038
physical_gate_long_range_recall = 0.088436
future_frustration_false_nucleus_rate = 0.309455
future_frustration_cluster_precision = 0.071567
future_frustration_long_range_recall = 0.070742
best_physical_term = burial_gain
worst_physical_term = future_frustration
active_physical_selection_survives = false
coupling_constraint_count = 745
coupling_rerank_false_nucleus_rate = 0.34375
coupling_rerank_cluster_precision = 0.072461
coupling_rerank_long_range_recall = 0.565164
coupling_rerank_real_vs_decoy_enrichment_ratio = 1.931956
coupling_trace_loop_selected_event_count = 54
coupling_trace_loop_false_nucleus_rate = 0.0
coupling_trace_loop_cluster_precision = 0.164128
coupling_trace_loop_long_range_recall = 0.397896
coupling_trace_loop_constraint_recall = 0.466536
coupling_trace_loop_real_vs_decoy_enrichment_ratio = 1.693887
coupling_selector_targets_met = true
external_evolutionary_couplings_used = false
oracle_constraint_control = true
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
real-coordinate visual contact benchmark
contact-law threshold falsification
folding-nucleus closure search
competitive nucleus selection
nucleus graph selectivity and decoy falsification
physical closure-state evaluator
active physical selection and term ablation
coupling-preserved nucleus selector
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
python3 scripts/run_real_coordinate_visual_benchmark.py
python3 scripts/run_contact_law_threshold_search.py
python3 scripts/run_folding_nucleus_closure_search.py
python3 scripts/run_competitive_nucleus_selection.py
python3 scripts/run_nucleus_graph_selectivity_benchmark.py
python3 scripts/run_physical_closure_state_benchmark.py
python3 scripts/run_active_physical_selection_benchmark.py
python3 scripts/run_evolutionary_coupling_nucleus_benchmark.py
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
real_coordinate_visual_8_*
real_coordinate_visuals/*/*
contact_law_threshold_*
folding_nucleus_closure_*
competitive_nucleus_selection_*
nucleus_graph_selectivity_*
physical_closure_state_*
active_physical_selection_*
coupling_nucleus_selector_*
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

Real-coordinate visual 8-row benchmark:

```text
real_coordinate_native_contacts_extracted = true
toy_locked_contact_targets_used = false
coarse_ca_only = true
full_atomic_folding_available = false
benchmark_size = 8
contact_map_f1_computed_count = 8
mean_contact_map_f1 = 0.051198
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

Contact-law threshold search:

```text
pair_feature_row_count = 94546
threshold_grid_row_count = 505
holdout_row_count = 32
current_scalar_score_best_global_threshold = 0.52
current_scalar_score_best_global_f1 = 0.150906
current_scalar_score_best_global_micro_f1 = 0.147772
current_scalar_score_threshold_stable = false
current_scalar_score_law_rejected = true
pair_only_best_f1 = 0.178009
best_law_candidate_model = pair_plus_entropy_score
best_law_candidate_loo_mean_test_f1 = 0.233569
best_law_candidate_loo_threshold_std = 0.007071
best_law_candidate_survives = false
law_generalizes = false
native_truth_used_before_feature_generation = false
row_specific_thresholds_forbidden = true
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

The candidate improves held-out F1 but fails the law boundary because it loses
long-range contact recall. The correct statement is law rejected / not yet
found, not folding solved.

Folding-nucleus closure search:

```text
candidate_closure_event_count = 5041
selected_threshold = 0.3
native_nucleus_recall = 0.50826
long_range_contact_recall_after_nucleus = 0.913242
pair_level_mean_long_range_contact_recall = 0.0
long_range_recall_delta_vs_pair_level = 0.913242
false_nucleus_rate = 0.692721
contact_cluster_precision = 0.032261
nucleus_level_long_range_beats_pair_level = true
nucleus_law_survives = false
native_truth_used_before_event_generation = false
row_specific_nucleus_thresholds_forbidden = true
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

The useful finding is that cooperative closure is a better proof object for
long-range native-region recovery. The blocking failure is false closure
overgeneration.

Competitive nucleus selection:

```text
pre_competition_event_count = 3846
post_competition_selected_event_count = 686
event_reduction_ratio = 0.821633
pre_false_nucleus_rate = 0.692721
post_false_nucleus_rate = 0.594592
pre_contact_cluster_precision = 0.032261
post_contact_cluster_precision = 0.044608
pre_long_range_contact_recall = 0.913242
post_long_range_contact_recall = 0.588907
selected_event_count_target_met = true
false_nucleus_rate_target_met = false
contact_cluster_precision_target_met = false
long_range_contact_recall_target_met = true
nucleus_competition_law_survives = false
native_truth_used_before_selection = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

The useful finding is narrower: competition can reduce closure flooding while
keeping some long-range contact recovery. It does not yet reject false nuclei
well enough to claim a mechanism.

Nucleus graph selectivity and decoy falsification:

```text
pre_graph_selected_event_count = 686
post_graph_selected_event_count = 320
pre_false_nucleus_rate = 0.594592
post_false_nucleus_rate = 0.609375
pre_contact_cluster_precision = 0.044608
post_contact_cluster_precision = 0.043311
pre_long_range_contact_recall = 0.588907
post_long_range_contact_recall = 0.372496
rank_enrichment_at_25 = 0.959248
decoy_native_overlap_rate = 0.484375
real_native_positive_rate = 0.390625
real_vs_decoy_enrichment_ratio = 0.806452
nucleus_graph_law_survives = false
native_truth_used_before_graph_selection = false
native_truth_used_before_decoy_matching = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

This is a negative but useful result. The graph selector reduces event count,
but matched decoys beat the selected graph. The missing variable is still not
captured.

Physical closure-state evaluator:

```text
candidate_state_count = 320
state_build_success_count = 320
state_build_failure_count = 0
mean_loop_strain = 0.012782
mean_steric_clash_score = 0.019839
mean_burial_gain = 0.583596
mean_unsatisfied_polar_penalty = 0.126116
mean_future_frustration_score = 0.282609
real_vs_decoy_physical_enrichment_ratio = 1.171535
physical_state_rank_enrichment_at_25 = 1.1264
post_physical_false_nucleus_rate = 0.609375
post_physical_contact_cluster_precision = 0.043311
post_physical_long_range_contact_recall = 0.372496
physical_enrichment_target_met = true
post_physical_false_nucleus_rate_target_met = false
post_physical_contact_cluster_precision_target_met = false
post_physical_long_range_contact_recall_target_met = true
physical_state_law_survives = false
native_truth_used_before_physical_scoring = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

This is the first physical-state crack attempt. It sees a real score difference
against decoys, but that difference is not enough to improve native-contact
truth. The missing variable is probably still deeper: orientation, side-chain
packing, solvent exposure, or kinetic order.

Active physical selection and term ablation:

```text
graph_only_false_nucleus_rate = 0.609375
physical_rerank_false_nucleus_rate = 0.559375
physical_gate_false_nucleus_rate = 0.316338
future_frustration_false_nucleus_rate = 0.309455
graph_only_cluster_precision = 0.043311
physical_rerank_cluster_precision = 0.050488
physical_gate_cluster_precision = 0.069038
future_frustration_cluster_precision = 0.071567
graph_only_long_range_recall = 0.372496
physical_rerank_long_range_recall = 0.405008
physical_gate_long_range_recall = 0.088436
future_frustration_long_range_recall = 0.070742
physical_rerank_real_vs_decoy_enrichment_ratio = 1.585882
best_physical_term = burial_gain
worst_physical_term = future_frustration
physical_terms_with_positive_ablation_effect = burial_gain
active_physical_selection_survives = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

This is an active-selector rejection, not a failed project. Reranking improves
some signals, but the gates trade away too much long-range recall. The current
physical proxy is informative, not sufficient.

Coupling-preserved nucleus selector:

```text
coupling_source_kind = coordinate_oracle_surrogate_for_missing_evolutionary_channel_v1
external_evolutionary_couplings_used = false
coordinate_truth_used_to_build_constraints = true
native_truth_used_before_coupling_selection = true
oracle_constraint_control = true
coupling_selector_targets_met = true
claim_allowed = false
```

This is the first selector surface that kills fake nuclei in the checked-in
coordinate-backed benchmark, but it does so with an oracle-control coupling
file. The next legitimate test is to replace
`data/folding_real_coordinate_visual_8_couplings.locked.json` with external
MSA/DCA couplings and rerun the exact same benchmark.

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
docs/COMPETITIVE_NUCLEUS_SELECTION_BOUNDARY.md competitive nucleus boundary
docs/NUCLEUS_GRAPH_SELECTIVITY_BOUNDARY.md nucleus graph selectivity boundary
docs/PHYSICAL_CLOSURE_STATE_BOUNDARY.md physical closure-state boundary
docs/ACTIVE_PHYSICAL_SELECTION_BOUNDARY.md active physical selection boundary
docs/EVOLUTIONARY_COUPLING_NUCLEUS_BOUNDARY.md coupling selector boundary
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
python3 scripts/run_real_coordinate_visual_benchmark.py
python3 scripts/run_contact_law_threshold_search.py
python3 scripts/run_folding_nucleus_closure_search.py
python3 scripts/run_competitive_nucleus_selection.py
python3 scripts/run_nucleus_graph_selectivity_benchmark.py
python3 scripts/run_physical_closure_state_benchmark.py
python3 scripts/run_active_physical_selection_benchmark.py
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
docs/COMPETITIVE_NUCLEUS_SELECTION_BOUNDARY.md
docs/NUCLEUS_GRAPH_SELECTIVITY_BOUNDARY.md
docs/PHYSICAL_CLOSURE_STATE_BOUNDARY.md
docs/ACTIVE_PHYSICAL_SELECTION_BOUNDARY.md
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
Competitive nucleus selection: event flooding is reduced, but false closures keep the law unclaimed.
Nucleus graph selectivity: selected graph cores fail against matched decoys.
Physical closure-state evaluator: physical score enrichment appears, but native/contact gates fail.
Active physical selection: rerank helps, gates over-prune, selector law rejected.
Protein folding: not solved.
Next correct work: improve visual/mechanism evidence quality without unlocking global claims.
```

## License

MIT.
