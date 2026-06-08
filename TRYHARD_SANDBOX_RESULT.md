# Try-hard sandbox execution result

## What was actually executed here

Command:

```bash
cd external_msa_free_predictors
ESMFOLD_API_TIMEOUT_SECONDS=30 \
LOCAL_ESMFOLD_TIMEOUT_SECONDS=30 \
OMEGAFOLD_TIMEOUT_SECONDS=30 \
RUN_ROOT=../tryhard_sandbox_validation/live_tryhard_runs \
./run_all_4ake_msa_free_tryhard.sh ..
```

Result:

```text
ESMFold API: failed because sandbox DNS cannot resolve api.esmatlas.com
Local ESMFold v1: failed because openfold dependency is not installed
OmegaFold: failed because sandbox DNS cannot resolve github.com/model source
Candidate PDBs found: 0
Probe reports tested: 0
folding_problem_solved: false
```

This is not a 4AKE biological failure. It is an execution-environment failure: no external predictor PDB could be generated inside this sandbox.

## Validation

Command:

```bash
PYTHONPATH=src PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
PHARMACOTOPOLOGY_TEST_TIMEOUT_SECONDS=70 \
timeout 180s python3 -m pytest -q -vv
```

Result:

```text
5 passed, 2 skipped in 15.56s
```

## Important fixes in this package

- ESMFold API runner now works with no args: `./esmfold_api.sh`
- ESMFold API uses raw `text/plain` POST body.
- Full all-source runner has hard timeouts.
- No GIF generation is triggered.
- AlphaFold-like sources are still rejected by the project probe unless explicitly overridden.
- Any PDB generated externally can be dropped into `external_msa_free_predictors/tryhard_runs/` and tested automatically.

## The real test command on a machine with internet/GPU

```bash
cd external_msa_free_predictors
./run_all_4ake_msa_free_tryhard.sh /path/to/pharmacotopology_4ake_tryhard_live_solution_full
```

Final verdict will be in:

```text
external_msa_free_predictors/tryhard_runs/final_tryhard_verdict.json
```
