#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src:${PYTHONPATH:-}"

python3 "$REPO_ROOT/scripts/run_v14_unified_protein_esperanto_grammar_readout_v0.py" \
  --out-dir "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V14_UNIFIED_PROTEIN_ESPERANTO_GRAMMAR_READOUT"
