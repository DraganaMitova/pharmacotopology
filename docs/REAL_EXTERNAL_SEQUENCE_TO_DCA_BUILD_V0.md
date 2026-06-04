# REAL_EXTERNAL_SEQUENCE_TO_DCA_BUILD_V0

This batch is the first real external sequence-family build for the frozen
8-row coordinate benchmark. It does not use coordinate contacts, AlphaFold
distance maps, native labels, or manual coupling fabrication.

## Data Sources

- PDBe SIFTS-derived Pfam mapping API:
  `https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{pdb_id}`
- InterPro Pfam full Stockholm alignments:
  `https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pfam_id}/?annotation=alignment:full`

The local machine initially did not have `plmc`, `ccmpred`, `jackhmmer`, or
`hmmscan`. `plmc` was then built from source in `/private/tmp/plmc-master` and
the best current locked artifact uses focus-mode plmc against real InterPro
Pfam full alignments:

```text
coupling_source_kind = external_msa_dca_plmc_v1
dca_tool_used = true
dca_tool_mode = focus_plmc_unweighted_pfam_full_sample_2000
```

The MI/APC fallback remains useful as a tool-missing fallback, but the locked
score batch below is true plmc pseudolikelihood DCA.

## Builder

```text
scripts/build_real_external_sequence_to_dca_v0.py
```

Default outputs:

```text
data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json
first_contact_clean_pharmacotopology_layer_run/external_sequence_mapping_v0.csv
first_contact_clean_pharmacotopology_layer_run/external_msa_build_log_v0.csv
first_contact_clean_pharmacotopology_layer_run/external_dca_build_log_v0.csv
```

The builder writes `external_evolutionary_couplings_used = true` only when at
least one external coupling constraint is present. Empty builds remain locked as
`external_evolutionary_couplings_used = false`.

## Real-Data Row Outcome

```text
1MBN:A available, 153 accepted constraints, Pfam PF00042
2LZM:A available, 164 accepted constraints, Pfam PF00959
1TEN:A available, 90 accepted constraints, Pfam PF00041
1CSP:A available, 67 accepted constraints, Pfam PF00313
1TIM:A available, 247 accepted constraints, Pfam PF00121
1PGA:A rejected, no PDBe Pfam mapping
4AKE:A available, 214 accepted constraints, Pfam PF00406/PF05191
1CLL:A available, 148 accepted constraints, Pfam PF13499
```

## Benchmark Result

Latest real external focus-plmc trace-loop run:

```text
result = external_channel_not_yet_supported
external_probe_passed = false
external_couplings_available_rows = 7
usable_external_rows = 7
external_constraint_count = 1083
external_real_false_nucleus_rate = 0.365272
physical_rerank_false_nucleus_rate = 0.559375
external_real_cluster_precision = 0.074096
physical_rerank_cluster_precision = 0.050488
external_real_long_range_recall = 0.194621
oracle_trace_loop_long_range_recall = 0.397896
external_real_vs_control_enrichment_ratio = 1.004367
external_real_beats_physical = true
external_real_beats_matched_controls = false
external_real_meets_oracle_recall_floor = false
claim_allowed = false
```

## Interpretation

The focus-plmc channel moves in the right direction versus physical reranking on
false nuclei, cluster precision, and long-range recall. It is not enough. It
still fails matched controls and remains below the oracle recall floor.

What is still missing:

- query-centered homolog search/MSA construction instead of broad Pfam full
  alignments
- a non-Pfam acquisition path for the `1PGA` protein G B1 row
- control-robust confidence calibration; confidence-permuted and low-confidence
  tail controls remain too close
- stronger quality gates for weak or noisy top-L tails

Claims remain locked.
