#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PYTHONPATH=src PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 scripts/run_4ake_iterative_distance_geometry_diffusion_v0.py \
  --out-dir first_contact_clean_pharmacotopology_layer_run/4ake_iterative_distance_geometry_diffusion_v0
