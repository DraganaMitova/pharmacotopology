#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"
python3 "$REPO_ROOT/scripts/run_v15_4ake_dynamic_grammar_bridge_v0.py" "$@"
