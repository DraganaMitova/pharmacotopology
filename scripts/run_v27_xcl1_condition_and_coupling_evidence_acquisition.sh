#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
RUN_ROOT="$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run"

# Keep acquisition fast; refresh zero-MD prerequisites only if missing.
if [ ! -f "$RUN_ROOT/V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT/v20_xcl1_state_specific_evidence_readout_certificate.json" ]; then
  bash "$REPO_ROOT/scripts/run_v20_xcl1_state_specific_evidence_readout.sh"
fi
if [ ! -f "$RUN_ROOT/V25_FAST_MECHANISM_EVIDENCE_SPRINT/v25_fast_mechanism_evidence_sprint_certificate.json" ]; then
  bash "$REPO_ROOT/scripts/run_v25_fast_mechanism_evidence_sprint.sh"
fi
if [ ! -f "$RUN_ROOT/V26_XCL1_STATE_SEPARATION_OPERATOR_TEST/v26_xcl1_state_separation_operator_test_certificate.json" ]; then
  bash "$REPO_ROOT/scripts/run_v26_xcl1_state_separation_operator_test.sh"
fi
if [ ! -f "$RUN_ROOT/V16_TRANSFER_DATA_PREFLIGHT/v16_transfer_data_preflight_certificate.json" ]; then
  bash "$REPO_ROOT/scripts/run_v16_transfer_data_preflight.sh"
fi

python3 "$REPO_ROOT/scripts/run_v27_xcl1_condition_and_coupling_evidence_acquisition_v0.py" "$@"
