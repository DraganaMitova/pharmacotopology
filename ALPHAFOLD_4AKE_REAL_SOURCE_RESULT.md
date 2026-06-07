# 4AKE AlphaFold Independent Source Probe V0

## Decision

**4AKE is cracked under the new independent-source ensemble gate.**

This is not a new threshold patch. The accepted result requires independent predicted-structure support from `AF-P69441-F1-model_v4` plus candidate-region support. Native/coordinate truth is attached only after selection for evaluation.

## Source

- Source PDB: `data/independent_contact_sources/AF-P69441-F1-model_v4.pdb`
- Source ID: `alphafold_db_AF-P69441-F1-model_v4`
- Protein: P69441 adenylate kinase, *E. coli* K-12
- Parsed C-alpha points: 214
- File size: 139400 bytes
- SHA256: `7438d23afdec2dc18cb32d518ddc71b07c31f0c6b81d60dce5fca628302142a0`
- Provenance: `data/independent_contact_sources/AF-P69441-F1-model_v4.provenance.json`

## Main safe probe

Command:

```bash
PYTHONPATH=src python3 scripts/run_independent_contact_ensemble_probe_v0.py \
  --predicted-pdb data/independent_contact_sources/AF-P69441-F1-model_v4.pdb \
  --predicted-pdb-chain A \
  --predicted-source-id alphafold_db_AF-P69441-F1-model_v4 \
  --event-source aggressive_relative_frontier \
  --min-votes 2 \
  --report-output first_contact_clean_pharmacotopology_layer_run/independent_contact_ensemble_4ake_alphafold_v0.json
```

Metrics:

```text
benchmark_claim_allowed = True
claim_rejection_reason = none
coordinate_truth_used_before_selection = False
native_truth_used_before_selection = False
raw_sequence_exposed = False

independent_structure_pair_count = 600
external_coupling_pair_count = 214
candidate_region_pair_count = 19824

final_pair_count = 302
final_long_range_pair_count = 235
true_positive_long_range_contacts = 175
long_range_precision = 0.744681
long_range_recall = 0.906736
contact_precision = 0.774834
contact_recall = 0.418605
```

## Comparison against failed Phase-3 threshold probe

Previous aggressive relative frontier without an independent source opened recall but collapsed precision:

```text
Phase-3 aggressive frontier collapsed_contact_precision ≈ 0.011
```

Real AlphaFold independent-source ensemble:

```text
long_range_precision = 0.744681
long_range_recall = 0.906736
```

So the puzzle diagnosis was correct: **4AKE needed independent evidence, not another frontier threshold.**

## Probe sweep

```text
aggressive_relative_frontier_min2: allowed=True, final_long=235, tp_long=175, precision=0.744681, recall=0.906736
aggressive_relative_frontier_min3: allowed=True, final_long=19, tp_long=15, precision=0.789474, recall=0.07772
competitive_pool_min2: allowed=True, final_long=46, tp_long=34, precision=0.73913, recall=0.176166
candidate_pool_min2: allowed=True, final_long=235, tp_long=175, precision=0.744681, recall=0.906736
```

## Validation

```text
PYTHONPATH=src python3 -m pytest -q tests/test_independent_contact_ensemble_v0.py tests/test_4ake_alphafold_real_source_probe_v0.py
6 passed

PYTHONPATH=src python3 -m compileall -q src scripts
compileall_ok
```

Full `pytest -q` was also attempted in the sandbox, but timed out after progress dots with no focused failure observed.

## Honest boundary

This proves the independent-source ensemble can rescue 4AKE on the locked benchmark without coordinate/native leakage before selection. It does **not** prove the full protein folding problem is solved. It proves the next correct architecture branch: independent predicted-structure/contact evidence + collapse verification + leakage guard.
