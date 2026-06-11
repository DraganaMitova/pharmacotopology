#!/usr/bin/env bash
set -euo pipefail

: "${REPO_ROOT:=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 "$REPO_ROOT/scripts/run_v32_external_constraint_source_import_preflight_v0.py"
