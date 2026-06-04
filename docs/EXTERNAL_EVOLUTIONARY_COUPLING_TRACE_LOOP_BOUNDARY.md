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
result = insufficient_external_signal / external_channel_not_yet_supported / external_channel_supported_in_v0
external_probe_passed = true/false
external_couplings_available_rows = N
external_rows_rejected_low_depth = N
external_rows_rejected_mapping = N
external_real_beats_physical = true/false
external_real_beats_matched_controls = true/false
mechanism_discovery_claim_allowed = false
folding_problem_solved = false
claim_allowed = false
```

The best possible V0 interpretation is:

```text
external evolutionary couplings preserve part of the anti-fake-nucleus signal under matched controls
```

It is not:

```text
folding solved
mechanism discovered
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
