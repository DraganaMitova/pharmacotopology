# Current Project Status

This is the current truth of the repository after the coupling-preserved
nucleus selector milestone.

## Latest Safety Baseline

```text
latest recorded safety commit = e8cc5eb Audit visual mechanism claims and freeze toy-contact benchmark
current target = Add coupling-preserved nucleus trace loop
```

The current target adds a locked coupling-constraint channel and an iterative
trace loop. It shows that native-selective constraints can kill fake burial and
pass selector survival gates, but the checked-in coupling file is an oracle
control derived from coordinate-native contacts. Mechanism-discovery claims
remain locked until the same benchmark passes with external MSA/DCA couplings.

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

## Contact-Law Threshold Status

The threshold search computes sequence-only pair features for every residue
pair in the eight coordinate-backed rows, then attaches native labels only for
scoring:

```text
pair_feature_row_count = 94546
native_pair_label_count = 2988
threshold_grid_row_count = 505
holdout_row_count = 32
native_truth_used_before_feature_generation = false
row_specific_thresholds_forbidden = true
```

The current scalar score is rejected as a stable law:

```text
current_scalar_score_best_global_threshold = 0.52
current_scalar_score_best_global_f1 = 0.150906
current_scalar_score_best_global_micro_f1 = 0.147772
current_scalar_score_threshold_std = 0.217715
current_scalar_score_threshold_stable = false
current_scalar_score_law_rejected = true
```

The best candidate is better, but does not survive as a law:

```text
pair_only_best_f1 = 0.178009
best_law_candidate_model = pair_plus_entropy_score
best_law_candidate_loo_mean_test_f1 = 0.233569
best_law_candidate_loo_threshold_std = 0.007071
best_law_candidate_survives = false
law_generalizes = false
candidate_law_failure_count = 4
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

Interpretation:

```text
Pair-plus-entropy improves held-out F1 and reduces false contacts, but it loses
long-range contact recall. The law is not discovered.
```

## Folding-Nucleus Closure Status

The nucleus search generates candidate closure events between sequence
segments, scores those events from sequence-only features, then attaches
coordinate-native labels after event generation:

```text
candidate_closure_event_count = 5041
nucleus_threshold_min = 0.3
selected_threshold = 0.3
accepted_event_count = 3846
native_truth_used_before_event_generation = false
native_label_attached_after_event_generation = true
row_specific_nucleus_thresholds_forbidden = true
```

The useful signal:

```text
native_nucleus_recall = 0.50826
long_range_contact_recall_after_nucleus = 0.913242
pair_level_mean_long_range_contact_recall = 0.0
long_range_recall_delta_vs_pair_level = 0.913242
nucleus_level_long_range_beats_pair_level = true
cooperative_closure_supported = true
```

The blocking failure:

```text
false_nucleus_rate = 0.692721
contact_cluster_precision = 0.032261
trap_event_count = 3056
candidate_law_failure_count = 6
nucleus_law_survives = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

Interpretation:

```text
Cooperative closure is a better object than independent pair contacts for
long-range native-region recovery, but the current closure score overgenerates
false nuclei. No folding law has been found.
```

## Competitive Nucleus Selection Status

The competitive layer uses sequence-only features to rank and filter the
accepted closure events. Native contact labels are attached only after
selection:

```text
pre_competition_event_count = 3846
post_competition_selected_event_count = 686
event_reduction_ratio = 0.821633
native_truth_used_before_selection = false
native_label_attached_after_selection = true
row_specific_nucleus_thresholds_forbidden = true
```

What improved:

```text
pre_long_range_contact_recall = 0.913242
post_long_range_contact_recall = 0.588907
selected_event_count_target_met = true
long_range_contact_recall_target_met = true
competition_reduces_event_flood = true
```

What still fails:

```text
pre_false_nucleus_rate = 0.692721
post_false_nucleus_rate = 0.594592
pre_contact_cluster_precision = 0.032261
post_contact_cluster_precision = 0.044608
false_nucleus_rate_target_met = false
contact_cluster_precision_target_met = false
nucleus_competition_law_survives = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

Interpretation:

```text
Competition reduces the event flood and keeps some long-range contact recovery.
It does not yet reject false nuclei or raise cluster precision enough to claim a
folding mechanism.
```

## Nucleus Graph Selectivity Status

The graph-selectivity layer treats selected closure events as graph edges
between segment nodes, then compares the selected graph core against matched
decoys from the same protein. Native labels are attached only after graph
selection and decoy matching:

```text
pre_graph_selected_event_count = 686
post_graph_selected_event_count = 320
native_truth_used_before_graph_selection = false
native_truth_used_before_decoy_matching = false
native_label_attached_after_graph_selection = true
finder_zip_allowed = false
```

What improved:

```text
post_graph_selected_event_target_met = true
post_long_range_contact_recall_target_met = true
pre_long_range_contact_recall = 0.588907
post_long_range_contact_recall = 0.372496
```

What failed:

```text
pre_false_nucleus_rate = 0.594592
post_false_nucleus_rate = 0.609375
pre_contact_cluster_precision = 0.044608
post_contact_cluster_precision = 0.043311
rank_enrichment_at_25 = 0.959248
decoy_native_overlap_rate = 0.484375
real_native_positive_rate = 0.390625
real_vs_decoy_enrichment_ratio = 0.806452
post_false_nucleus_rate_target_met = false
post_contact_cluster_precision_target_met = false
decoy_enrichment_target_met = false
nucleus_graph_law_survives = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

Interpretation:

```text
The selected graph core is not enriched over matched decoys. This is a useful
negative result: graph scoring reduced count but did not discover the missing
selector law.
```

## Physical Closure-State Status

The physical-state evaluator converts graph-selected closures into coarse
sequence-only state packets. Native labels are attached only after physical
scoring:

```text
candidate_state_count = 320
state_build_success_count = 320
state_build_failure_count = 0
native_truth_used_before_physical_scoring = false
native_label_attached_after_physical_scoring = true
```

What improved:

```text
mean_physical_state_score = 0.41558
mean_decoy_physical_state_score = 0.354731
real_vs_decoy_physical_enrichment_ratio = 1.171535
real_beats_decoy_physical_score_rate = 0.65625
physical_state_rank_enrichment_at_25 = 1.1264
physical_enrichment_target_met = true
post_physical_long_range_contact_recall_target_met = true
```

What still fails:

```text
post_physical_false_nucleus_rate = 0.609375
post_physical_contact_cluster_precision = 0.043311
post_physical_long_range_contact_recall = 0.372496
post_physical_false_nucleus_rate_target_met = false
post_physical_contact_cluster_precision_target_met = false
physical_state_law_survives = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

Interpretation:

```text
Coarse physical scoring sees a real score difference against decoys, but the
difference does not translate into enough native-contact precision. The missing
variable is likely deeper than this coarse state proxy.
```

## Active Physical Selection Status

The active selector compares four sequence-only selection surfaces, then
attaches native labels only after selection:

```text
graph selector only
graph + active physical rerank
graph + physical viability gate
graph + physical viability gate + future-frustration gate
native_truth_used_before_active_selection = false
native_label_attached_after_active_selection = true
```

What improved:

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
physical_rerank_real_vs_decoy_enrichment_ratio = 1.585882
```

What still fails:

```text
physical_gate_long_range_recall = 0.088436
future_frustration_long_range_recall = 0.070742
best_physical_term = burial_gain
worst_physical_term = future_frustration
physical_terms_with_positive_ablation_effect = burial_gain
active_physical_selection_survives = false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
```

Interpretation:

```text
Physical reranking is informative, but active viability gates trade away too
much long-range recovery. The current physical selector is rejected. This is not
a folding mechanism discovery.
```

## Coupling-Preserved Nucleus Status

The coupling selector tests the missing native-selective constraint channel:

```text
direct coupling support
future closure preservation
matched-decoy coupling margin
iterative coupling trace loop
```

The checked-in coupling file is an oracle-control layer:

```text
coupling_constraint_count = 745
external_evolutionary_couplings_used = false
per_constraint_coordinate_truth_used = true
coordinate_truth_used_to_build_constraints = true
native_truth_used_before_coupling_selection = true
oracle_constraint_control = true
claim_mode_validation_passed = false
```

What improved:

```text
physical_rerank_false_nucleus_rate = 0.559375
coupling_rerank_false_nucleus_rate = 0.34375
coupling_trace_loop_false_nucleus_rate = 0.0

physical_rerank_cluster_precision = 0.050488
coupling_rerank_cluster_precision = 0.072461
coupling_trace_loop_cluster_precision = 0.164128

physical_rerank_long_range_recall = 0.405008
coupling_rerank_long_range_recall = 0.565164
coupling_trace_loop_long_range_recall = 0.397896

coupling_trace_loop_constraint_recall = 0.466536
coupling_trace_loop_real_vs_decoy_enrichment_ratio = 1.693887
coupling_selector_targets_met = true
```

What remains locked:

```text
mechanism_discovery_claim_allowed = false
global_folding_claim_allowed = false
folding_problem_solved = false
claim_allowed = false
```

Interpretation:

```text
The missing selector channel is real enough to test: native-selective
constraints plus a trace-building loop remove fake nuclei in the oracle
control. The next hard test is replacing the oracle coupling file with
external evolutionary couplings and rerunning:

EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_V0
```

The current concrete blocker is:

```text
REAL_EXTERNAL_COUPLING_FILE_BUILD_V0
```

This build stage creates the provenance-locked external coupling file, target
manifest, and build log before the selector benchmark runs. The current checked
build artifact is allowed to contain zero usable rows because no raw MSA/DCA
coupling file has been supplied.

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
contact-law threshold falsification
folding-nucleus closure search
competitive nucleus selection
nucleus graph selectivity and decoy falsification
physical closure-state evaluator
active physical selection and term ablation
coupling-preserved nucleus selector
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
scripts/run_contact_law_threshold_search.py
scripts/run_folding_nucleus_closure_search.py
scripts/run_competitive_nucleus_selection.py
scripts/run_nucleus_graph_selectivity_benchmark.py
scripts/run_physical_closure_state_benchmark.py
scripts/run_active_physical_selection_benchmark.py
scripts/run_evolutionary_coupling_nucleus_benchmark.py
scripts/build_real_external_coupling_file_v0.py
scripts/build_external_evolutionary_coupling_trace_loop_v0.py
scripts/run_external_evolutionary_coupling_trace_loop_benchmark.py
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
contact_law_threshold_*
folding_nucleus_closure_*
competitive_nucleus_selection_*
nucleus_graph_selectivity_*
physical_closure_state_*
active_physical_selection_*
coupling_nucleus_selector_*
external_coupling_trace_loop_*
external_coupling_target_manifest_v0.json
external_coupling_build_log_v0.csv
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
current scalar contact score is rejected as a stable threshold law
best threshold candidate does not generalize because long-range recall collapses
cooperative closure recovers long-range native regions but overgenerates traps
competitive nucleus selection reduces events but does not pass false-rate gates
nucleus graph selectivity fails against matched decoys
physical closure-state score enrichment does not pass native/contact gates
active physical selection rerank helps but gates over-prune long-range recall
coupling trace loop is oracle-control until external couplings replace coordinate-derived constraints
external trace-loop V0 requires a supplied provenance-locked MSA/DCA coupling file
legacy feature modules still contain useful low-level primitives
```

The current state is useful as a safe falsification, abstention, visual
mechanism-inspection, and contact-repair stack, not as a predictive folding
system.

## Next Allowed Research Target

Do not add more broad biology claims before the canonical surface stays clean.

The next research target can be:

```text
false-nucleus rejection
cluster-precision improvement
decoy-beating graph selectivity
physical state variables that improve native-contact precision
external MSA/DCA coupling replacement for the oracle-control coupling layer
```

Allowed shape:

```text
improve contact evidence packets
preserve truth-after-prediction scoring
keep global fold-class claims locked
track coverage cost
record uncertainty, failure cohorts, and visible contact-map drift
analyze coordinate-native missed contact cohorts
explain long-range recall collapse under pair-plus-entropy scoring
reduce false segment closures without destroying long-range recall
add steric/orientation/trap filters at the event level, not pair level
keep visual_12_is_toy_benchmark = true
keep contact_repair_overfit_risk_reported = true
keep toy_locked_contact_targets_used = false on real-coordinate surfaces
keep current_scalar_score_law_rejected = true until a stronger law survives
keep nucleus_law_survives = false until false_nucleus_rate is under control
keep nucleus_competition_law_survives = false until false-rate and precision gates pass
keep nucleus_graph_law_survives = false until matched decoys are beaten
keep physical_state_law_survives = false until physical enrichment improves native/contact gates
keep coupling oracle_constraint_control = true from creating mechanism claims
keep external trace-loop V0 claim_allowed = false on the 8-row probe
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
