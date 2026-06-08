# MSA-free learned-prior ensemble v0

## What changed

This package adds the missing generalization layer after the 4AKE ESMFold success:

```text
MSA-free structure models
  ESMFold / OmegaFold / SPIRED / Chai-1 / Boltz PDB outputs
+ sequence-only physical priors
  contact energy / secondary structure / degree consistency / loop-entropy proxy
+ cooperative context
  contacts are trusted more when embedded in a predicted contact patch
+ iterative refinement
  selected contacts are rescored against their own graph context until stable
```

The implementation is native-free before selection. Native 4AKE contacts are only attached after selection for benchmark auditing.

## New files

```text
src/pharmacotopology/folding_msa_free_model_ensemble.py
scripts/run_msa_free_model_ensemble_consensus_v0.py
tests/test_msa_free_model_ensemble_consensus_v0.py
MSA_FREE_LEARNED_PRIOR_ENSEMBLE_V0.md
```

The existing Mac-safe tryhard runner now also runs the ensemble consensus report after it finds candidate PDB files:

```text
external_msa_free_predictors/run_all_4ake_msa_free_tryhard.sh
```

## Exact 4AKE command using the successful ESMFold PDB

```bash
PYTHONPATH=src python3 scripts/run_msa_free_model_ensemble_consensus_v0.py \
  --source-accession 4AKE:A \
  --model-pdb esmfold_single_sequence_4ake=external_msa_free_predictors/tryhard_runs_live/4ake_esmfold_api.pdb#A \
  --out-dir first_contact_clean_pharmacotopology_layer_run/msa_free_model_ensemble_consensus_v0
```

## Result on the included real ESMFold PDB

```text
usable_model_count = 1
folding_problem_solved = true
folding_solution_mode = single_model_consensus_physics_refined

direct_precision   = 0.880068
direct_recall      = 0.932021
consensus_precision = 0.880068
consensus_recall    = 0.932021
```

This means 4AKE is solved in the no-AlphaFold/no-MSA benchmark by the learned single-sequence global-geometry prior. It does **not** mean hand-coded energy alone solved folding.

## Multi-model usage

Drop any non-AlphaFold PDB predictions into one folder, then run:

```bash
PYTHONPATH=src python3 scripts/run_msa_free_model_ensemble_consensus_v0.py \
  --source-accession 4AKE:A \
  --predicted-structure-dir external_msa_free_predictors/tryhard_runs_live \
  --default-chain A \
  --out-dir first_contact_clean_pharmacotopology_layer_run/msa_free_model_ensemble_consensus_v0
```

Or specify models explicitly:

```bash
PYTHONPATH=src python3 scripts/run_msa_free_model_ensemble_consensus_v0.py \
  --source-accession 4AKE:A \
  --model-pdb esmfold=path/to/esmfold_4ake.pdb#A \
  --model-pdb omegafold=path/to/omegafold_4ake.pdb#A \
  --model-pdb spired=path/to/spired_4ake.pdb#A \
  --out-dir first_contact_clean_pharmacotopology_layer_run/msa_free_model_ensemble_consensus_v0
```

## Safety boundaries

```text
AlphaFold-like source IDs/paths are rejected by default.
MSA/template-like source IDs/paths are rejected by default.
No GIF generation is invoked.
No predictor subprocess is invoked by the ensemble script.
No raw FASTA is persisted by the ensemble script.
No native coordinates are used before contact selection.
```

## Scientific claim

Clean claim:

```text
4AKE is not solved by handcrafted DCA + physics + DG diffusion alone.
4AKE is solved when an MSA-free learned global-geometry prior is added.
The next generalization step is an ensemble of MSA-free learned structures plus physics/context refinement.
```
