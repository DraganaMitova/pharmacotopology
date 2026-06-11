#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
RUN_ROOT="$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run"

if [ ! -f "$RUN_ROOT/V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_DECISION/v29_mechanism_operator_panel_summary_and_md_readiness_certificate.json" ]; then
  bash "$REPO_ROOT/scripts/run_v29_mechanism_operator_panel_summary_and_md_readiness_decision.sh"
fi

python3 "$REPO_ROOT/scripts/run_v30_external_constraint_and_coupling_acquisition_sprint_v0.py" "$@"
