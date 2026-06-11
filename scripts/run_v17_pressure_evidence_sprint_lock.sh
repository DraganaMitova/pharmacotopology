#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}"

# One-command sprint: refresh public target/context files, run V16 preflight,
# zero-MD role classification, V16 gap lock, then V17 evidence sprint decision.
# No MD or threshold tuning is executed here.
bash "$REPO_ROOT/scripts/download_v16_pressure_target_structures.sh"
bash "$REPO_ROOT/scripts/run_v16_transfer_data_preflight.sh"
bash "$REPO_ROOT/scripts/run_v16_zero_md_role_transfer_readout.sh"
bash "$REPO_ROOT/scripts/run_v16_pressure_evidence_gap_lock.sh"
python3 "$REPO_ROOT/scripts/run_v17_pressure_evidence_sprint_lock_v0.py" "$@"
