#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

bash "$REPO_ROOT/scripts/run_v21_pressure_evidence_panel_summary.sh"
python3 "$REPO_ROOT/scripts/run_v22_external_evidence_and_annotation_acquisition_panel_v0.py" "$@"
