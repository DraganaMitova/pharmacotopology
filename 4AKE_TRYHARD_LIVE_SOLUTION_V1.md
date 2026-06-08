# 4AKE try-hard live solution package v1

This is the package that refuses to close the problem prematurely.

## What changed

The previous package targeted ESMFold v1 only and had a weak runner interface. This version fixes that:

- `external_msa_free_predictors/esmfold_api.sh` now works with no arguments.
- The ESMFold API call uses raw `text/plain` POST body, matching current API examples.
- Every external call is bounded by `timeout`.
- The all-sources runner tries ESMFold API, local ESMFold v1, OmegaFold, and externally dropped PDBs.
- It then runs the pharmacotopology no-AF probe for every produced PDB.
- It produces one final verdict JSON.

## Main command

```bash
cd external_msa_free_predictors
./run_all_4ake_msa_free_tryhard.sh /path/to/pharmacotopology
```

## Where to drop manual results

Put any externally generated PDB under:

```text
external_msa_free_predictors/tryhard_runs/
```

Then rerun the same command. It will scan and test all PDBs.

## Success threshold

```text
precision >= 0.70 and recall >= 0.70
```

## Sandbox result

The included sandbox validation proves the runner does not hang, but this sandbox cannot resolve DNS for `api.esmatlas.com`, `github.com`, `huggingface.co`, or `files.rcsb.org`. Therefore the sandbox could not generate the learned PDB. The package is built to run on an internet/GPU machine.
