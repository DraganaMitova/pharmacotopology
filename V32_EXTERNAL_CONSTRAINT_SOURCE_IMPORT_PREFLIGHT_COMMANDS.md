# V32 External Constraint Source Import Preflight Commands

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 -m compileall -q src scripts
python3 -m pytest -q tests/test_v31_constraint_backed_operator_readout_preflight.py tests/test_v32_external_constraint_source_import_preflight.py

bash "$REPO_ROOT/scripts/run_v32_external_constraint_source_import_preflight.sh"
python3 "$REPO_ROOT/scripts/print_v32_external_constraint_source_import_preflight.py"

python3 "$REPO_ROOT/scripts/export_locked_runtime_certificates_v0.py"
```

If V32 writes a template manifest, fill this file with real external source rows and rerun V32:

```text
data/external_constraints/v32_external_constraint_source_import_manifest.json
```

Do not use files from `first_contact_clean_pharmacotopology_layer_run/` as external evidence.
