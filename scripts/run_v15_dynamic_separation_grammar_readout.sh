#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

python3 "$REPO_ROOT/scripts/run_v15_dynamic_separation_grammar_readout_v0.py" \
  --out-dir "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V15_DYNAMIC_SEPARATION_GRAMMAR_READOUT"
