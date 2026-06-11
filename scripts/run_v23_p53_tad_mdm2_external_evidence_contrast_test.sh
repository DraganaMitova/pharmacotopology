#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

bash "$REPO_ROOT/scripts/run_v22_external_evidence_and_annotation_acquisition_panel.sh"
bash "$REPO_ROOT/scripts/run_v18_p53_tad_mdm2_partner_induced_evidence_test.sh"
bash "$REPO_ROOT/scripts/run_v18b_p53_isolated_tad_abstain_context.sh"
python3 "$REPO_ROOT/scripts/run_v23_p53_tad_mdm2_external_evidence_contrast_test_v0.py" "$@"
