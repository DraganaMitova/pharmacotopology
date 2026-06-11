#!/usr/bin/env bash
set -euo pipefail

export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 scripts/build_kcsa_rcsb_external_constraint_import_v0.py
python3 scripts/run_v32_external_constraint_source_import_preflight_v0.py
python3 scripts/run_v33_constraint_backed_operator_readout_v0.py
python3 scripts/run_v33_negative_controls_v0.py
python3 scripts/run_v34_kcsa_discriminative_content_controls_v0.py
python3 scripts/print_v33_negative_controls.py
python3 scripts/print_v34_kcsa_discriminative_content_controls.py
python3 -m pytest -q tests/test_v33_constraint_backed_operator_readout.py tests/test_v33_negative_controls.py tests/test_v34_kcsa_discriminative_content_controls.py
