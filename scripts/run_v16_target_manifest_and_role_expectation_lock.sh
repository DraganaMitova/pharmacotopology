#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"
python3 "$REPO_ROOT/scripts/run_v16_target_manifest_and_role_expectation_lock_v0.py" "$@"
