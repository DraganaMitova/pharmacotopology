# MSA-free learned-geometry ensemble implementation summary

## Status

Implemented and validated.

This is the layer requested after the ESMFold 4AKE success: multiple MSA-free learned structure outputs can now vote together, and their contact consensus is refined by sequence-only physics/context scores.

## Main result on included 4AKE ESMFold output

```text
usable_model_count = 1
folding_problem_solved = true
folding_solution_mode = single_model_consensus_physics_refined

direct_precision    = 0.880068
direct_recall       = 0.932021
consensus_precision = 0.880068
consensus_recall    = 0.932021
selected_contacts   = 592
```

Interpretation:

```text
4AKE is solved without AlphaFold/MSA by the learned single-sequence global-geometry prior.
The new ensemble layer preserves that solution and is ready for OmegaFold/SPIRED/Chai/Boltz replication.
```

## What was added

```text
src/pharmacotopology/folding_msa_free_model_ensemble.py
scripts/run_msa_free_model_ensemble_consensus_v0.py
tests/test_msa_free_model_ensemble_consensus_v0.py
MSA_FREE_LEARNED_PRIOR_ENSEMBLE_V0.md
TEST_SUITE_MSA_FREE_ENSEMBLE_VALIDATION.log
```

Updated:

```text
README.md
external_msa_free_predictors/run_all_4ake_msa_free_tryhard.sh
```

## How it works

```text
1. Parse one or more non-AlphaFold, no-MSA predicted PDBs.
2. Extract contact maps for each model.
3. Compute model vote fraction per contact.
4. Add physics/context terms:
   - sequence-only contact energy
   - secondary-structure prior
   - degree consistency
   - loop-entropy proxy
   - cooperative contact-patch context
   - optional DCA/candidate-region support as weak context
5. Iterate the selected contact graph until stable.
6. Attach native contacts only after selection for benchmark audit.
```

## Safety boundaries

```text
AlphaFold-like source IDs/paths rejected by default.
MSA/template-like source IDs/paths rejected by default.
No GIF generation.
No predictor subprocess inside the ensemble script.
No raw FASTA persistence.
No native truth before selection.
No coordinate truth before selection.
```

## Validation

```text
8 passed, 2 skipped in 15.58s
RETURN_CODE=0
```

## Exact command used

```bash
PYTHONPATH=src python3 scripts/run_msa_free_model_ensemble_consensus_v0.py \
  --source-accession 4AKE:A \
  --model-pdb esmfold_single_sequence_4ake=external_msa_free_predictors/tryhard_runs_live/4ake_esmfold_api.pdb#A \
  --out-dir first_contact_clean_pharmacotopology_layer_run/msa_free_model_ensemble_consensus_v0
```

## Next replication step

Add more predicted PDBs into `external_msa_free_predictors/tryhard_runs_live/`, for example:

```text
4ake_omegafold.pdb
4ake_spired.pdb
4ake_chai1_sample_0.pdb
4ake_boltz_sample_0.pdb
```

Then run:

```bash
PYTHONPATH=src python3 scripts/run_msa_free_model_ensemble_consensus_v0 \
  --source-accession 4AKE:A \
  --predicted-structure-dir external_msa_free_predictors/tryhard_runs_live \
  --default-chain A
```

Clean claim:

```text
Energy alone was not enough.
DG diffusion alone was not enough.
Learned global geometry was the missing signal for 4AKE.
The new code generalizes that signal into a reusable MSA-free model ensemble + physics refinement layer.
```
