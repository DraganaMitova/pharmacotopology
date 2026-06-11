#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
python3 "$REPO_ROOT/scripts/run_v19_kcsa_membrane_pore_evidence_lock_v0.py" "$@"
