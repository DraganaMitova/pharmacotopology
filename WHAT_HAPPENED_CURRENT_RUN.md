# What happened in the current try-hard run

## Direct answer

The project did not hang. The full no-AF/MSA-free try-hard runner executed to completion with hard timeouts.

The sandbox did not produce a real ESMFold/OmegaFold/SPIRED PDB because the container cannot resolve external predictor hosts and local ESMFold is missing `openfold`.

Therefore this run did not biologically solve 4AKE; it produced a safe, bounded, no-fake-success verdict with zero candidate PDBs tested.

## Current validation commands

```bash
cd /mnt/data/pharmacotopology_4ake_tryhard_live_solution_full
PYTHONPATH=src PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
PHARMACOTOPOLOGY_TEST_TIMEOUT_SECONDS=70 \
timeout 180s python3 -m pytest -q -vv
```

Result:

```text
5 passed, 2 skipped in 15.06s
```

Try-hard runner command:

```bash
cd external_msa_free_predictors
ESMFOLD_API_TIMEOUT_SECONDS=30 \
LOCAL_ESMFOLD_TIMEOUT_SECONDS=30 \
OMEGAFOLD_TIMEOUT_SECONDS=30 \
RUN_ROOT=../current_validation/live_tryhard_runs \
timeout 140s ./run_all_4ake_msa_free_tryhard.sh ..
```

Result:

```text
ESMFold API: no PDB, DNS failure for api.esmatlas.com from sandbox container
Local ESMFold v1: no PDB, ModuleNotFoundError: openfold
OmegaFold: no PDB in this sandbox
Candidate PDBs found: 0
Probe reports tested: 0
folding_problem_solved: false
```

## Current final verdict JSON

```json
{
  "all": [],
  "best": [],
  "folding_problem_solved": false,
  "kind": "tryhard_final_verdict_v1",
  "solved_count": 0,
  "tested_reports": 0
}
```

## What is now fixed for the real next run

- `external_msa_free_predictors/esmfold_api.sh` works with no arguments.
- ESMFold API runner posts raw sequence as `text/plain`.
- Every external predictor path has hard timeouts.
- Runner scans all produced/dropped `.pdb` files and tests them automatically.
- No GIF generation is used.
- AlphaFold-like source IDs remain rejected by default in the no-AF probe.

## Real machine command

On a machine with internet, run:

```bash
cd external_msa_free_predictors
./run_all_4ake_msa_free_tryhard.sh /path/to/pharmacotopology_4ake_tryhard_live_solution_full
```

If an external PDB is already produced by Colab/GPU, drop it into:

```text
external_msa_free_predictors/tryhard_runs/
```

and rerun the same command. The final verdict appears at:

```text
external_msa_free_predictors/tryhard_runs/final_tryhard_verdict.json
```
