# V33 constraint-backed operator readout commands

Use this only after V32 says:

```text
preflight_status: V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED
selected_V33_target: KcsA
selected_V33_panel: V33_CONSTRAINT_BACKED_OPERATOR_READOUT
claim_allowed: False
new_MD_allowed: False
```

## Paste commands

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 scripts/run_v32_external_constraint_source_import_preflight_v0.py
python3 scripts/print_v32_external_constraint_source_import_preflight.py
python3 scripts/run_v33_constraint_backed_operator_readout_v0.py
python3 scripts/print_v33_constraint_backed_operator_readout.py
python3 -m pytest -q tests/test_v33_constraint_backed_operator_readout.py
```

Expected V33 status:

```text
readout_status: V33_CONSTRAINT_BACKED_OPERATOR_READOUT_PASSED_CLAIM_DISABLED
selected_V33_target: KcsA
claim_allowed: False
new_MD_allowed: False
positive_folding_evidence_found: False
folding_problem_solved: False
```

Boundary: this is an operator readout, not a de novo folding solution.
