# REAL_EXTERNAL_SEQUENCE_TO_DCA_BUILD_V0

This batch is the first real external sequence-family build for the frozen
8-row coordinate benchmark. It does not use coordinate contacts, AlphaFold
distance maps, native labels, or manual coupling fabrication.

## Data Sources

- PDBe SIFTS-derived Pfam mapping API:
  `https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{pdb_id}`
- InterPro Pfam full Stockholm alignments:
  `https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pfam_id}/?annotation=alignment:full`
- EBI HMMER jackhmmer API for the `1PGA:A` protein G B1 rescue route:
  `https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/jackhmmer`

The local machine initially did not have `plmc`, `ccmpred`, `jackhmmer`, or
`hmmscan`. `plmc` was then built from source in `/private/tmp/plmc-master` and
the best current locked artifact uses focus-mode plmc against real InterPro
Pfam full alignments plus an EBI HMMER jackhmmer UniProt alignment for the
`1PGA:A` core:

```text
coupling_source_kind = external_msa_dca_plmc_v1
dca_tool_used = true
dca_tool_mode = focus_plmc_unweighted_pfam_full_plus_hmmer_jackhmmer_1pga_iter3
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
1PGA:A available, 56 accepted constraints, EBI HMMER jackhmmer/PF01378
4AKE:A available, 214 accepted constraints, Pfam PF00406/PF05191
1CLL:A available, 148 accepted constraints, Pfam PF13499
```

`1PGA:A` still has no PDBe Pfam mapping in this batch. The accepted route is
sequence-only: EBI HMMER jackhmmer iteration 3 against UniProt produced a
Stockholm alignment with 491 external domain rows; an exact aligned UniProt hit
maps the `1PGA:A` core positions 6-56 with identity 1.0, target coverage
0.910714, and focus mapping confidence 1.0. The plmc focus run had 492 valid
sequences over 51 focus sites, so the row passes the serious depth gate with
effective_sequence_count_over_length = 8.785714.

## Benchmark Result

Latest real external focus-plmc trace-loop run:

```text
result = external_channel_not_yet_supported
external_probe_passed = false
external_couplings_available_rows = 8
usable_external_rows = 8
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
claim_allowed = false
```

## Interpretation

The HMMER-assisted focus-plmc channel moves in the right direction versus
physical reranking on false nuclei, cluster precision, and long-range recall.
Adding `1PGA:A` pushes long-range recall above the 50 percent oracle recall
floor, but it also raises the false nucleus rate and still fails matched
controls. It is not enough.

What is still missing:

- query-centered homolog search/MSA construction for every row, not just the
  `1PGA` rescue
- control-robust confidence calibration; confidence-permuted and low-confidence
  tail controls remain too close
- stronger quality gates for weak or noisy top-L tails

Claims remain locked.
