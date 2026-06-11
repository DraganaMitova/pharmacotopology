#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
RUN_ROOT="$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run"

if [ ! -f "$RUN_ROOT/V25_FAST_MECHANISM_EVIDENCE_SPRINT/v25_fast_mechanism_evidence_sprint_certificate.json" ]; then
  bash "$REPO_ROOT/scripts/run_v25_fast_mechanism_evidence_sprint.sh"
fi
if [ ! -f "$RUN_ROOT/V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST/v28_xcl1_state_condition_evidence_contrast_certificate.json" ]; then
  bash "$REPO_ROOT/scripts/run_v28_xcl1_state_condition_evidence_contrast_test.sh"
fi

python3 "$REPO_ROOT/scripts/run_v29_mechanism_operator_panel_summary_and_md_readiness_decision_v0.py" "$@"
