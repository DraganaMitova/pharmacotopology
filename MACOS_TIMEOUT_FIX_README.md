# macOS try-hard runner fix

The previous runner used Linux `timeout`. macOS does not ship that command, so ESMFold/OmegaFold were never actually called on the Mac run.

This package replaces Linux `timeout` with `external_msa_free_predictors/run_with_timeout.py`, a cross-platform Python timeout wrapper.

Run:

```bash
cd external_msa_free_predictors
chmod +x *.sh *.py
rm -rf tryhard_runs_live
RUN_ROOT="$PWD/tryhard_runs_live" \
ESMFOLD_API_TIMEOUT_SECONDS=600 \
LOCAL_ESMFOLD_TIMEOUT_SECONDS=1200 \
OMEGAFOLD_TIMEOUT_SECONDS=1200 \
./run_all_4ake_msa_free_tryhard.sh ..
```

If ESMFold API succeeds, you should see a candidate PDB and tested_count > 0.
