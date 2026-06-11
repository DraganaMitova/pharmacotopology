# KcsA real external source import — V32 → V33 door

Run these from inside the project folder. No git. No MD. No full suite.

```bash
export REPO_ROOT="$(pwd)"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"
python3 scripts/build_kcsa_rcsb_external_constraint_import_v0.py
python3 scripts/run_v32_external_constraint_source_import_preflight_v0.py
python3 scripts/print_v32_external_constraint_source_import_preflight.py
```

Expected after the builder succeeds:

```text
preflight_status: V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED
selected_V33_target: KcsA
selected_V33_panel: V33_CONSTRAINT_BACKED_OPERATOR_READOUT
claim_allowed: False
new_MD_allowed: False
new_MD_recommended: False
```

Boundary:

```text
This imports external RCSB PDB 1BL8 coordinate-derived pore/filter and assembly/interface constraints.
It tests source-import/readout gating only.
It is not a de novo KcsA folding benchmark and not a universal folding claim.
```
