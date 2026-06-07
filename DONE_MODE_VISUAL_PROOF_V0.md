# DONE_MODE_VISUAL_PROOF_V0

This repo is now cleaned around one active regression path:

```text
4AKE:A + real AlphaFold independent source + ensemble/collapse + leakage gate + visual proof GIF
```

## What "done mode" means

This does **not** claim universal de novo folding for every possible protein.

It means the system now has one stable interface that can be run for any locked target:

```bash
PYTHONPATH=src python3 scripts/run_done_mode_contact_prediction_v0.py \
  --source-accession 4AKE:A \
  --predicted-pdb data/independent_contact_sources/AF-P69441-F1-model_v4.pdb \
  --predicted-pdb-chain A \
  --predicted-source-id alphafold_db_AF-P69441-F1-model_v4 \
  --render-visual-proof
```

The system either:

1. allows a benchmark claim when independent evidence exists and leakage gates stay closed, or
2. abstains/rejects when the independent source is missing or unsafe.

This is the honest version of "predict every protein": every target can pass through the same predictor contract, but the system does not fake a result when evidence is insufficient.

## 4AKE result

```text
benchmark_claim_allowed = true
claim_rejection_reason = none
coordinate_truth_used_before_selection = false
native_truth_used_before_selection = false
raw_sequence_exposed = false

long_range_precision = 0.744681
long_range_recall    = 0.906736
contact_precision    = 0.774834
contact_recall       = 0.418605
```

## Visual proof outputs

Primary GIF:

```text
first_contact_clean_pharmacotopology_layer_run/visual_proofs/4ake_alphafold_ensemble_visual_proof_v0.gif
```

DONE-mode GIF:

```text
first_contact_clean_pharmacotopology_layer_run/done_mode_v0/4ake_a/4ake_visual_proof.gif
```

The GIF has five frames:

1. candidate region pool
2. external DCA/coupling signal
3. AlphaFold independent structure contacts
4. final ensemble/collapse contacts
5. post-selection evaluation: TP / FP / missed native long-range contacts

Native truth is used only in the final visualization/evaluation frame after selection.

## Abstention example

`1MBN:A` without independent evidence produces:

```text
status = abstained_or_rejected
claim_rejection_reason = missing_independent_structure_source
```

This keeps the system from pretending it can solve a weak-signal target without a second source.

## Active test suite

Only one active test remains:

```text
tests/test_4ake_alphafold_real_source_end_to_end_v0.py
```

It runs the done-mode wrapper, the ensemble probe, and the GIF renderer with a 30-second hard timeout.
