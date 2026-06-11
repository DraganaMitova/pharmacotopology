#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
PDB="$REPO_ROOT/data/v16_pressure_targets/p53_TAD_MDM2/1YCR.pdb"
if [ ! -f "$PDB" ]; then
  bash "$REPO_ROOT/scripts/download_v16_pressure_target_structures.sh"
fi
python3 "$REPO_ROOT/scripts/run_v18_p53_tad_mdm2_partner_induced_evidence_test_v0.py" "$@"
