#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts"
CLEAN=1 REPLICAS=1 STEPS=1000 TARGET_OPEN_STEPS=100 REPORTER_INTERVAL_STEPS=100 \
  bash "$REPO_ROOT/scripts/run_v13a_1ubq_repair_fixed.sh"
python3 "$REPO_ROOT/scripts/check_v13_runtime_provenance_v0.py" \
  "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_REPAIR_FIXED"
