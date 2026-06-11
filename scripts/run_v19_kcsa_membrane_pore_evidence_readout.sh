#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

bash "$REPO_ROOT/scripts/download_v16_pressure_target_structures.sh"
bash "$REPO_ROOT/scripts/run_v16_transfer_data_preflight.sh"
bash "$REPO_ROOT/scripts/run_v16_zero_md_role_transfer_readout.sh"
python3 "$REPO_ROOT/scripts/run_v19_kcsa_membrane_pore_evidence_readout_v0.py" "$@"
