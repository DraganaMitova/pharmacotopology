#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

SRC="$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13b_1CLL_CALMODULIN_FIXED"
OUT="$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13b_1CLL_ADAPTIVE_CHEMICAL_POLICY_READOUT"

if [[ ! -d "$SRC" ]]; then
  echo "missing source run directory: $SRC" >&2
  echo "first run: CLEAN=1 REPLICAS=10 MAX_PARALLEL=2 CPU_THREADS=4 bash \"$REPO_ROOT/scripts/run_v13b_1cll_calmodulin_fixed.sh\"" >&2
  exit 2
fi

traj_count=$(find "$SRC" -name "openmm_dca_restrained_trajectory.pdb" | wc -l | tr -d ' ')
if [[ "$traj_count" == "0" ]]; then
  echo "no replica trajectories found in: $SRC" >&2
  exit 2
fi

rm -rf "$OUT"

python3 "$REPO_ROOT/scripts/run_v13a_purpose_gate_readout_v0.py" \
  --source-run-dir "$SRC" \
  --out-dir "$OUT" \
  --source-accession "1CLL:A" \
  --benchmark-file "$REPO_ROOT/data/folding_real_coordinate_visual_8.locked.json" \
  --external-coupling-file "$REPO_ROOT/data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json" \
  --anchor-profile-file "$REPO_ROOT/data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json" \
  --audit-pairs-json "$REPO_ROOT/data/audit_1cll_v13b.json" \
  --topology-mode "interdomain" \
  --domain-boundaries "1-75,80-148" \
  --chemical-policy "adaptive_soft_guard" \
  --frequency-grid "0.50:0.98:0.01" \
  --chemical-grid "0.00:0.50:0.05"
