#!/usr/bin/env bash
set -euo pipefail
export REPO_ROOT="$(pwd)"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"
python3 scripts/build_kcsa_rcsb_external_constraint_import_v0.py
python3 scripts/run_v32_external_constraint_source_import_preflight_v0.py
python3 scripts/print_v32_external_constraint_source_import_preflight.py
