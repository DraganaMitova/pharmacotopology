#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts"

python3 "$REPO_ROOT/scripts/run_garage_role_aware_rescue_selector_v9.py" \
  --source-run-dir "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_REPAIR_FIXED" \
  --benchmark-file "$REPO_ROOT/data/folding_real_coordinate_holdout_1ubq.locked.json" \
  --out-dir "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/comparison_v9/1UBQ_FIXED"

python3 "$REPO_ROOT/scripts/run_garage_role_aware_rescue_selector_v9.py" \
  --source-run-dir "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13b_1CLL_CALMODULIN_FIXED" \
  --benchmark-file "$REPO_ROOT/data/folding_real_coordinate_visual_8.locked.json" \
  --out-dir "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/comparison_v9/1CLL_FIXED"
