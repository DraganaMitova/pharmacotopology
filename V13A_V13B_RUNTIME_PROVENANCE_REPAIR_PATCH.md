# V13a/V13b runtime provenance repair patch

This patch is intentionally narrow. It does **not** tune role rules, selector thresholds, biological scoring, or the V12 locked contract logic.

## What was actually broken

The failed V13a/V13b runs did not produce trajectories. The replica logs show file/path failures before MD could run, for example missing target PDB files. That means the downstream V9 selector correctly reported `no replica trajectories found`.

So the current failure class is:

```text
replica_runtime_input_failure -> no trajectories -> postprocess has nothing to read
```

not:

```text
biological_transfer_failure
```

## What this patch changes

### Modified

- `scripts/run_openmm_tmd_replicas_v0.py`
  - Adds an input preflight before launching replica workers.
  - Writes `input_preflight.json` or `input_preflight_failed.json`.
  - Fails early if required files are missing.
  - Fails early if 1UBQ is pointed at an obvious calmodulin target file.
  - Adds `--input-preflight-only`.
  - Adds `--allow-target-pdb-mismatch` for deliberate debugging only.
  - Records `target_pdb` and `input_preflight` inside the final certificate.

- `scripts/run_openmm_dca_closure_md_v0.py`
  - Adds child-runner file checks for benchmark, coupling, target PDB, and resume PDB.
  - This prevents hidden `FileNotFoundError` tracebacks inside every replica.

### Added

- `scripts/download_v13_independent_target_pdbs.sh`
  - Downloads raw AlphaFold target PDBs.
  - Calls the segment preparer.

- `scripts/prepare_v13_independent_target_pdbs.py`
  - Trims raw AlphaFold models to the exact benchmark sequence segment.
  - Renumbers the segment to 1..N so the runner sees the correct residue indexing.
  - Writes provenance JSON next to each prepared target PDB.

- `scripts/run_v13a_1ubq_repair_fixed.sh`
  - Clean runnable wrapper for 1UBQ V13a repair.
  - Uses the prepared mature ubiquitin segment target PDB.

- `scripts/run_v13b_1cll_calmodulin_fixed.sh`
  - Clean runnable wrapper for 1CLL V13b.
  - Uses the prepared calmodulin segment target PDB.

- `scripts/run_v13_v9_postprocess_fixed.sh`
  - Runs V9 postprocessing only after fixed run directories exist.

- `scripts/check_v13_runtime_provenance_v0.py`
  - Reads a run directory and reports whether trajectories, certificate, preflight, and V10 selection outputs exist.

## How to run

From repo root:

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts"

bash "$REPO_ROOT/scripts/download_v13_independent_target_pdbs.sh"

CLEAN=1 bash "$REPO_ROOT/scripts/run_v13a_1ubq_repair_fixed.sh"
CLEAN=1 bash "$REPO_ROOT/scripts/run_v13b_1cll_calmodulin_fixed.sh"

python3 "$REPO_ROOT/scripts/check_v13_runtime_provenance_v0.py" \
  "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_REPAIR_FIXED"

python3 "$REPO_ROOT/scripts/check_v13_runtime_provenance_v0.py" \
  "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13b_1CLL_CALMODULIN_FIXED"

bash "$REPO_ROOT/scripts/run_v13_v9_postprocess_fixed.sh"
```

## Honest interpretation rules

If `successful_replicas = 0`, do not interpret biology.

If `trajectory_count = 0`, V9 cannot postprocess.

If input preflight says blocked, call it:

```text
preflight_abstain_missing_or_invalid_runtime_input
```

Only if input preflight is ready and trajectories exist should you interpret V5/V10 transfer results.

Even then, if `effective_strict_count = 0` remains, keep:

```text
claim_allowed = false
physics_interpretation_allowed = false
```
