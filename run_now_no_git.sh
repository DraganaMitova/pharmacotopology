#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")" || exit 1
python3 scripts/run_clean_start_sanity_v0.py
python3 scripts/run_v32_external_constraint_source_import_preflight_v0.py
python3 scripts/print_v32_external_constraint_source_import_preflight.py
