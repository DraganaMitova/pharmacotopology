# Clean pharmacotopology repo — V13 runtime provenance repair

This zip is meant to be used as a clean replacement repo, not as a tiny patch overlay.

## What was kept

- `src/`
- `scripts/`
- `tests/`
- `data/`
- `docs/`
- `external_msa/`
- `external_msa_free_predictors/`
- core project files: `README.md`, `pyproject.toml`, `LICENSE`, `.gitignore`, `sitecustomize.py`

## What was removed

- `.venv/`
- `.git/`
- `.pytest_cache/`, `__pycache__/`, macOS `__MACOSX` / `.DS_Store`
- generated run output folders such as `first_contact_clean_pharmacotopology_layer_run/`
- top-level generated logs and old one-off run artifacts

## V13 repair files included

- `scripts/run_openmm_tmd_replicas_v0.py`
  - same role-aware selector/audit structure preserved
  - adds early input preflight so missing/wrong PDBs stop before launching doomed replicas
- `scripts/run_openmm_dca_closure_md_v0.py`
  - adds clear required-file checks
- `scripts/prepare_v13_independent_target_pdbs.py`
- `scripts/download_v13_independent_target_pdbs.sh`
- `scripts/run_v13a_1ubq_repair_fixed.sh`
- `scripts/run_v13b_1cll_calmodulin_fixed.sh`
- `scripts/run_v13_v9_postprocess_fixed.sh`
- `scripts/check_v13_runtime_provenance_v0.py`

## How to install cleanly on your Mac

Rename your old folder first, do not delete it immediately:

```bash
mv "/Users/draganamitova/My Projects/pharmacotopology" \
   "/Users/draganamitova/My Projects/pharmacotopology_OLD_DIRTY_BACKUP"
```

Then unzip this zip into:

```text
/Users/draganamitova/My Projects/
```

You should end up with:

```text
/Users/draganamitova/My Projects/pharmacotopology
```

Then run:

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts"

python3 -m compileall -q src scripts
```

## V13a / V13b run order

First prepare the independent target PDB files:

```bash
bash "$REPO_ROOT/scripts/download_v13_independent_target_pdbs.sh"
```

Then run 1UBQ:

```bash
CLEAN=1 bash "$REPO_ROOT/scripts/run_v13a_1ubq_repair_fixed.sh"
```

Then check provenance:

```bash
python3 "$REPO_ROOT/scripts/check_v13_runtime_provenance_v0.py" \
  "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_REPAIR_FIXED"
```

Then run 1CLL:

```bash
CLEAN=1 bash "$REPO_ROOT/scripts/run_v13b_1cll_calmodulin_fixed.sh"
```

Then check provenance:

```bash
python3 "$REPO_ROOT/scripts/check_v13_runtime_provenance_v0.py" \
  "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13b_1CLL_CALMODULIN_FIXED"
```

Only after trajectories exist, run V9 postprocess:

```bash
bash "$REPO_ROOT/scripts/run_v13_v9_postprocess_fixed.sh"
```

## Important interpretation rule

If `input_preflight_failed.json` appears, the result is not biological. It is a runtime/provenance/input failure.

If replicas produce no trajectories, V9 cannot diagnose biology because there is no trajectory evidence to postprocess.
