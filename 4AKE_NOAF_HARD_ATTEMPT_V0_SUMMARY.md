# 4AKE no-AlphaFold hard attempt v0

## Result

This package does **not** claim that 4AKE is solved without AlphaFold yet.

The hard attempt is now automated, bounded, and tested:

- GIF generation stays off by default.
- Default pytest no longer runs the heavy self-consistency regression or visual/GIF test.
- Slow/full tests are opt-in.
- The internet MSA-free path is implemented and bounded.
- AlphaFold-like sources are still rejected for no-AF claims unless explicitly allowed as a positive-control.
- A direct-structure metric was added so a real ESMFold/OmegaFold/SPIRED PDB can be judged without losing recall through the candidate-region intersection filter.

## What happened in this sandbox

`ESM Atlas ESMFold` request was attempted through the new bounded internet runner, but the shell/Python runtime in this sandbox cannot resolve external DNS:

```text
URLError: Temporary failure in name resolution
```

Therefore no real non-AlphaFold MSA-free PDB was produced here, and `folding_problem_solved = false`.

## Why this is still progress

The missing piece is now executable: a real single-sequence structure predictor output can be passed into:

```bash
PYTHONPATH=src python3 scripts/run_msa_free_single_sequence_structure_probe_v0.py \
  --source-accession 4AKE:A \
  --predicted-pdb /path/to/esmfold_or_omegafold_or_spired_4ake.pdb \
  --predicted-source-id esmfold_single_sequence_4ake \
  --predicted-pdb-chain A \
  --include-sequence-physical-priors \
  --out-dir first_contact_clean_pharmacotopology_layer_run/msa_free_single_sequence_structure_v0
```

The probe now reports both:

1. `report`: ensemble candidate-region + coupling + independent-structure votes.
2. `direct_structure_metric`: all contacts implied by the single-sequence predicted PDB.

A no-AlphaFold folding claim is allowed only if the source is not AlphaFold-like and direct precision/recall pass the thresholds.

## Positive-control safety check

The bundled AlphaFold fixture was run only as a blocked positive-control. It proves the new direct metric can detect a good structure, while refusing to count AlphaFold as a no-AF solution:

```text
alphafold_used_by_this_script = true
benchmark_claim_allowed_by_direct_structure = false
folding_problem_solved = false
```

Its direct metric was high:

```text
native_contact_precision = 0.865
native_contact_recall = 0.928444
contact_map_f1 = 0.8956
long_range_contact_recall = 0.906736
```

That is not a no-AF success. It is a safety positive-control.

## Commands validated

Fast/default suite:

```bash
PYTHONPATH=src PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 timeout 120s python3 -m pytest -q -vv
```

Result:

```text
4 passed, 2 skipped in 14.03s
```

Opt-in full non-visual suite:

```bash
PYTHONPATH=src PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 timeout 150s python3 -m pytest -q -vv --run-full-suite
```

Observed after fixes:

```text
5 passed, 1 skipped in 26.02s
```

Hard attempt orchestrator:

```bash
PYTHONPATH=src PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 scripts/run_4ake_noaf_hard_attempt_v0.py \
  --out-dir first_contact_clean_pharmacotopology_layer_run/4ake_noaf_hard_attempt_v0
```

## Internet/source notes

The implemented internet path targets public or key-backed MSA-free single-sequence structure services:

- ESM Atlas foldSequence endpoint: `https://api.esmatlas.com/foldSequence/v1/pdb/`
- BioLM ESMFold endpoint: `https://biolm.ai/api/v3/esmfold/predict/` when `BIOLM_API_KEY` is supplied.
- Hugging Face documents `facebook/esmfold_v1` as an end-to-end folding model that does not require lookup/MSA.

Relevant source pages checked during this attempt:

- `https://huggingface.co/facebook/esmfold_v1`
- `https://zitniklab.hms.harvard.edu/ToolUniverse/_modules/tooluniverse/esmfold_tool.html`
- `https://biolm.ai/models/esmfold/`
- `https://www.rcsb.org/docs/additional-resources/structure-prediction`

## Bottom line

We did not fake victory. The project now has the missing execution slot and the direct structure scoring needed for the next real test. In this sandbox, internet DNS prevented obtaining a real ESMFold/OmegaFold/SPIRED PDB, so the honest result remains:

```text
4AKE no-AlphaFold solved = false
reason = no real non-AlphaFold single-sequence PDB was available in this runtime
```
