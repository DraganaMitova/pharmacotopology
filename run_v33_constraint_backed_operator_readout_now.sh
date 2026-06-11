#!/usr/bin/env bash
set -euo pipefail
export REPO_ROOT="${REPO_ROOT:-$(pwd)}"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 scripts/run_v32_external_constraint_source_import_preflight_v0.py
python3 scripts/print_v32_external_constraint_source_import_preflight.py
python3 scripts/run_v33_constraint_backed_operator_readout_v0.py
python3 scripts/print_v33_constraint_backed_operator_readout.py
python3 -m pytest -q tests/test_v33_constraint_backed_operator_readout.py
