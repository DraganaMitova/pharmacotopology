# External Evolutionary Coupling Trace Loop Boundary

This batch is:

```text
EXTERNAL_EVOLUTIONARY_COUPLING_TRACE_LOOP_V0
```

It tests whether externally supplied MSA/DCA-style sequence-covariation
couplings preserve any of the anti-fake-nucleus signal that the oracle coupling
control exposed.

It does not solve protein folding, discover a folding mechanism, design
proteins, create molecules, or allow clinical use.

## Frozen Benchmark

V0 keeps the current 8-row coordinate benchmark frozen:

```text
1MBN:A  myoglobin
2LZM:A  T4 lysozyme
1TEN:A  tenascin FN3
1CSP:A  cold shock protein
1TIM:A  triosephosphate isomerase
1PGA:A  protein G B1
4AKE:A  adenylate kinase
1CLL:A  calmodulin
```

Every row must be preregistered in the external coupling file. No row may be
silently dropped or replaced.

## Real File Build V0

Before running the selector benchmark, build the external input surface:

```text
REAL_EXTERNAL_COUPLING_FILE_BUILD_V0
```

The builder script is:

```text
scripts/build_real_external_coupling_file_v0.py
```

It writes:

```text
data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json
first_contact_clean_pharmacotopology_layer_run/external_coupling_target_manifest_v0.json
first_contact_clean_pharmacotopology_layer_run/external_coupling_build_log_v0.csv
```

The builder does not run the selector and does not fabricate couplings. If no
raw external MSA/DCA coupling file is supplied, it emits an empty locked coupling
file and marks all rows with honest build-log rejection reasons.

The build-stage statuses are more specific than the selector-stage statuses:

```text
external_couplings_available
external_couplings_rejected_no_sequence_mapping
external_couplings_rejected_low_msa_depth
external_couplings_rejected_low_coverage
external_couplings_rejected_domain_boundary_ambiguous
external_couplings_rejected_position_mapping_ambiguous
external_couplings_rejected_tool_failed
external_couplings_rejected_coordinate_taint
```

V0 rejects duplicate residue-pair couplings. A real external DCA output must
contain unique `(row_id, i, j)` pairs:

```text
reject_duplicate_coupling_pairs = true
duplicate_count_dropped = 0
```

V0 also requires every external constraint to carry the full provenance surface:

```text
row_id
source_accession
i
j
sequence_separation
normalized_separation
confidence
raw_score
apc_corrected_score
rank
rank_fraction
source_kind
msa_source_kind
msa_sha256
msa_depth
effective_sequence_count
effective_sequence_count_over_length
target_coverage
focus_sequence_mapping_confidence
coordinate_truth_used_to_build_constraint
native_truth_used_before_coupling_selection
structure_model_used
raw_sequence_exposed
```

Missing field means the file is rejected.

## Accepted Sources

Only these source kinds are accepted:

```text
external_msa_dca_plmc_v1
external_msa_dca_ccmpred_v1
external_evcouplings_sequence_covariation_v1
external_pfam_msa_dca_v1
external_uniref_msa_dca_v1
```

These source kinds are rejected:

```text
coordinate_native_contacts
pdb_contact_map
alphafold_distance_map
supervised_structure_contact_predictor
manual_contact_annotation
oracle_from_benchmark_coordinates
oracle_from_native_contacts
```

## Row Statuses

Each preregistered row receives exactly one external-coupling status:

```text
external_couplings_available
external_couplings_rejected_low_depth
external_couplings_rejected_low_coverage
external_couplings_rejected_mapping_ambiguous
external_couplings_rejected_coordinate_taint
```

Low-depth rows are honest biological failures, not project failures. They mean
there was not enough usable covariation signal for the V0 selector.

## Quality Gates

The serious V0 gate is:

```text
target_coverage >= 0.70
focus_sequence_mapping_confidence >= 0.98
effective_sequence_count_over_length >= 5.0
top_L_couplings_available = true
coordinate_truth_used_to_build_constraint = false
native_truth_used_before_coupling_selection = false
structure_model_used = false
raw_sequence_exposed = false
```

If any constraint is coordinate-derived, native-truth-derived, or
structure-model-derived before selection, the row is rejected and the imported
dataset remains oracle-tainted for claim gating.

## Controls

Every real external coupling set is attacked by matched controls:

```text
external_shuffled_same_row_same_separation
external_confidence_permuted
external_cross_row_swapped
external_random_long_range_same_count
external_low_confidence_tail
physical_no_coupling_baseline
oracle_coordinate_positive_control
```

The same-row same-separation control is the critical one: it preserves the easy
long-range separation distribution so the real external signal has to win on
more than distance.

## V0 Success

The V0 report records:

```text
result = no_external_data_built / insufficient_external_signal / external_channel_not_yet_supported / external_channel_supported_in_v0
external_probe_passed = true/false
external_couplings_available_rows = N
external_rows_rejected_low_depth = N
external_rows_rejected_mapping = N
external_real_beats_physical = true/false
external_real_beats_matched_controls = true/false
score_margin_expansion_row_count = N
score_margin_expansion_row_count_margin_vs_matched_controls = N
score_margin_expansion_row_count_margin_vs_adversarial_controls = N
score_margin_expansion_repeated_independent_row_signal_seen = true/false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
claim_allowed = false
```

The best possible V0 interpretation is:

```text
external evolutionary couplings preserve part of the anti-fake-nucleus signal under matched controls
```

If score-margin expansion candidates repeat across independent benchmark rows,
that is recorded as stability evidence for the next selector design. The guarded
score-margin-expanded selector is only interpretable when it beats matched and
adversarial controls without native-truth filtering, and it remains claim-locked
even when that selector-level comparison succeeds.

It is not:

```text
folding solved
mechanism discovered
```

## REAL_EXTERNAL_SEQUENCE_TO_DCA_BUILD_V0 Result

The first real external sequence-family run used PDBe SIFTS-derived Pfam
mappings, InterPro Pfam full alignments, EBI HMMER jackhmmer for the `1PGA:A`
protein G B1 row, and focus-mode `plmc` built locally from source. It did not
use coordinate contacts or native labels to build constraints.

```text
external_couplings_available_rows = 8
external_constraint_count = 1139
external_real_false_nucleus_rate = 0.390272
physical_rerank_false_nucleus_rate = 0.559375
external_real_cluster_precision = 0.082689
physical_rerank_cluster_precision = 0.050488
external_real_long_range_recall = 0.275503
oracle_trace_loop_long_range_recall = 0.397896
external_real_vs_control_enrichment_ratio = 0.997262
external_real_beats_physical = true
external_real_beats_matched_controls = false
external_real_meets_oracle_recall_floor = true
external_margin_gated_false_nucleus_rate = 0.320991
external_margin_gated_cluster_precision = 0.088524
external_margin_gated_long_range_recall = 0.280851
external_margin_gated_vs_control_enrichment_ratio = 0.99413
external_margin_gated_beats_physical = true
external_margin_gated_beats_matched_controls = true
external_margin_gated_meets_oracle_recall_floor = true
claim_allowed = false
```

Interpretation:

```text
real external sequence-family couplings are no longer empty
focus-plmc improves over physical rerank
margin-gated trace-loop beats matched controls on false-rate and precision
confidence/enrichment controls are still too close
1PGA is now covered by a sequence-only HMMER jackhmmer/PF01378 rescue route
claims remain locked
```

## Persistent Trace Recovery Diagnostic

The KNOT-Core admission rule separates evidence from claims and requires
falsification controls before a stronger word is allowed. The local KNOT
particle simulation similarly treats a loop as more meaningful only when it
persists as a supported track instead of appearing as an isolated fragment.

The folding analogue added here is:

```text
coupling_trace_loop_persistent_rank_consistent_cluster_gated
```

This selector keeps the rank-consistent calibrated core, then allows a
lower-cluster recovery event only if it has nearby compatible recovery
candidates with external-coupling support, low blocked-future pressure, and a
minimum persistence score. Native labels are still attached only after event
generation and scoring.

Checked-in V0 result:

```text
external_rank_consistent_cluster_gated_selected_event_count = 23
external_rank_consistent_cluster_gated_false_nucleus_rate = 0.0
external_rank_consistent_cluster_gated_cluster_precision = 0.162044
external_rank_consistent_cluster_gated_long_range_recall = 0.208635
external_rank_consistent_cluster_gated_vs_control_enrichment_ratio = 1.291737

external_persistent_rank_consistent_cluster_gated_selected_event_count = 24
external_persistent_rank_consistent_cluster_gated_false_nucleus_rate = 0.0
external_persistent_rank_consistent_cluster_gated_cluster_precision = 0.167578
external_persistent_rank_consistent_cluster_gated_long_range_recall = 0.231865
external_persistent_rank_consistent_cluster_gated_vs_control_enrichment_ratio = 1.245257
external_persistent_rank_consistent_cluster_gated_vs_adversarial_calibrated_enrichment_ratio = 0.991698
external_persistent_rank_consistent_cluster_gated_real_vs_decoy_coupling_nucleus_enrichment_ratio = 1.591368
external_persistent_rank_consistent_cluster_gated_vs_control_nucleus_score_enrichment_ratio = 1.262676
external_persistent_rank_consistent_cluster_gated_vs_adversarial_calibrated_nucleus_score_enrichment_ratio = 1.049653
external_persistent_rank_consistent_cluster_gated_recovered_event_count = 1
external_persistent_rank_consistent_cluster_gated_recovered_native_long_range_contact_count = 21
external_persistent_rank_consistent_cluster_gated_probe_passed = false
external_persistent_rank_consistent_cluster_gated_selector_score_probe_passed = true
external_persistent_rank_consistent_cluster_gated_hard_selector_score_probe_passed = true
external_persistent_rank_consistent_cluster_gated_recall_frontier_count = 23
external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_count = 2
external_persistent_rank_consistent_cluster_gated_score_margin_expansion_candidate_native_long_range_contact_count = 38
external_persistent_rank_consistent_cluster_gated_score_margin_expansion_claim_allowed = false

external_score_margin_expanded_selected_event_count = 26
external_score_margin_expanded_added_event_count = 2
external_score_margin_expanded_added_native_long_range_contact_count = 38
external_score_margin_expanded_added_false_event_count = 0
external_score_margin_expanded_false_nucleus_rate = 0.0
external_score_margin_expanded_cluster_precision = 0.179688
external_score_margin_expanded_long_range_recall = 0.252423
external_score_margin_expanded_long_range_recall_delta_vs_persistent = 0.020558
external_score_margin_expanded_long_range_recall_margin_vs_matched_controls = 0.117799
external_score_margin_expanded_long_range_recall_margin_vs_adversarial_controls = 0.174647
external_score_margin_expanded_beats_matched_controls = true
external_score_margin_expanded_beats_adversarial_calibrated_controls = true
external_score_margin_expanded_claim_allowed = false

external_boundary_continuity_expanded_selected_event_count = 28
external_boundary_continuity_expanded_added_event_count = 2
external_boundary_continuity_expanded_added_native_long_range_contact_count = 2
external_boundary_continuity_expanded_added_false_event_count = 0
external_boundary_continuity_expanded_false_nucleus_rate = 0.0
external_boundary_continuity_expanded_cluster_precision = 0.174316
external_boundary_continuity_expanded_long_range_recall = 0.259572
external_boundary_continuity_expanded_long_range_recall_delta_vs_score_margin = 0.007149
external_boundary_continuity_expanded_long_range_recall_margin_vs_matched_controls = 0.124948
external_boundary_continuity_expanded_long_range_recall_margin_vs_adversarial_controls = 0.119626
external_boundary_continuity_expanded_beats_matched_controls = true
external_boundary_continuity_expanded_beats_adversarial_calibrated_controls = true
external_boundary_continuity_expanded_claim_allowed = false

external_edge_continuity_expanded_selected_event_count = 30
external_edge_continuity_expanded_added_event_count = 2
external_edge_continuity_expanded_added_native_long_range_contact_count = 3
external_edge_continuity_expanded_added_false_event_count = 0
external_edge_continuity_expanded_false_nucleus_rate = 0.0
external_edge_continuity_expanded_cluster_precision = 0.172363
external_edge_continuity_expanded_long_range_recall = 0.269988
external_edge_continuity_expanded_long_range_recall_delta_vs_boundary_continuity = 0.010416
external_edge_continuity_expanded_long_range_recall_margin_vs_matched_controls = 0.135364
external_edge_continuity_expanded_long_range_recall_margin_vs_adversarial_controls = 0.130042
external_edge_continuity_expanded_beats_matched_controls = true
external_edge_continuity_expanded_beats_adversarial_calibrated_controls = true
external_edge_continuity_expanded_claim_allowed = false

external_pressure_release_expanded_selected_event_count = 32
external_pressure_release_expanded_added_event_count = 2
external_pressure_release_expanded_added_native_long_range_contact_count = 13
external_pressure_release_expanded_added_false_event_count = 0
external_pressure_release_expanded_false_nucleus_rate = 0.0
external_pressure_release_expanded_cluster_precision = 0.169434
external_pressure_release_expanded_long_range_recall = 0.287035
external_pressure_release_expanded_long_range_recall_delta_vs_edge_continuity = 0.017047
external_pressure_release_expanded_long_range_recall_margin_vs_matched_controls = 0.152411
external_pressure_release_expanded_long_range_recall_margin_vs_adversarial_controls = 0.147089
external_pressure_release_expanded_beats_matched_controls = true
external_pressure_release_expanded_beats_adversarial_calibrated_controls = true
external_pressure_release_expanded_claim_allowed = false

external_registry_extension_expanded_selected_event_count = 35
external_registry_extension_expanded_added_event_count = 3
external_registry_extension_expanded_added_native_long_range_contact_count = 13
external_registry_extension_expanded_added_false_event_count = 0
external_registry_extension_expanded_false_nucleus_rate = 0.0
external_registry_extension_expanded_cluster_precision = 0.160862
external_registry_extension_expanded_long_range_recall = 0.292161
external_registry_extension_expanded_long_range_recall_delta_vs_pressure_release = 0.005126
external_registry_extension_expanded_long_range_recall_margin_vs_matched_controls = 0.157537
external_registry_extension_expanded_long_range_recall_margin_vs_adversarial_controls = 0.152215
external_registry_extension_expanded_beats_matched_controls = true
external_registry_extension_expanded_beats_adversarial_calibrated_controls = true
external_registry_extension_expanded_claim_allowed = false
claim_allowed = false
```

Interpretation:

```text
persistence recovers a real long-range trace event without reintroducing fake nuclei
persistence improves precision and long-range recall over the strict rank gate
the added event weakens the decoy-enrichment margin
raw coupling-only adversarial calibrated enrichment remains too close
the selector's full coupling-nucleus score clears matched and adversarial enrichment controls
the score-margin-expanded selector admits two addable native-positive trace events
that add 38 long-range contacts without adding a false nucleus
the score-margin-expanded selector beats matched and adversarial controls on precision and recall
the boundary-continuity selector adds two more long-range contacts without adding a false nucleus
the boundary-continuity selector improves recall but lowers cluster precision from 0.179688 to 0.174316
the edge-continuity selector adds three more long-range contacts without adding a false nucleus
the edge-continuity selector improves recall again but lowers cluster precision to 0.172363
the pressure-release selector adds thirteen more long-range contacts without adding a false nucleus
the pressure-release selector improves recall again but lowers cluster precision to 0.169434
the registry-extension selector adds thirteen more long-range contacts without adding a false nucleus
the registry-extension selector improves recall again but lowers cluster precision to 0.160862
claims remain locked
```

## Biological Warning

An external coupling unsupported by the monomer PDB contact map still counts as
a false positive for this selector benchmark. It may still be biologically
meaningful as an oligomer-interface, ligand-mediated, alternate-conformation,
or allosteric/family-level signal. The constraint audit fields therefore keep
these ideas separate:

```text
native_contact_supported
monomer_coordinate_unsupported
possible_interdomain_or_allosteric_signal
possible_oligomer_interface_signal
benchmark_counts_as_false_positive
```

For this benchmark, unsupported still counts as false. The science note stays
honest about why.
