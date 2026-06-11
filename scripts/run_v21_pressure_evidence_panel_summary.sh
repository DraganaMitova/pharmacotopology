#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

bash "$REPO_ROOT/scripts/download_v16_pressure_target_structures.sh"
bash "$REPO_ROOT/scripts/run_v18_p53_tad_mdm2_partner_induced_evidence_test.sh"
bash "$REPO_ROOT/scripts/run_v18_p53_partner_induced_evidence_lock.sh"
bash "$REPO_ROOT/scripts/run_v18b_p53_isolated_tad_abstain_context.sh"
bash "$REPO_ROOT/scripts/run_v19_kcsa_membrane_pore_evidence_readout.sh"
bash "$REPO_ROOT/scripts/run_v19_kcsa_membrane_pore_evidence_lock.sh"
bash "$REPO_ROOT/scripts/run_v20_xcl1_state_specific_evidence_readout.sh"
python3 "$REPO_ROOT/scripts/run_v21_pressure_evidence_panel_summary_v0.py" "$@"
