#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src:${PYTHONPATH:-}"

python3 "$REPO_ROOT/scripts/run_v13b_hierarchical_purpose_topology_readout_v0.py" \
  --source-run-dir "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13b_1CLL_CALMODULIN_FIXED" \
  --out-dir "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13b_1CLL_HIERARCHICAL_PURPOSE_TOPOLOGY_READOUT" \
  --source-accession "1CLL:A" \
  --benchmark-file "$REPO_ROOT/data/folding_real_coordinate_visual_8.locked.json" \
  --external-coupling-file "$REPO_ROOT/data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json" \
  --anchor-profile-file "$REPO_ROOT/data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json" \
  --audit-pairs-json "$REPO_ROOT/data/audit_1cll_v13b.json" \
  --domain-boundaries "1-75,80-148" \
  --frequency-grid "0.50:0.98:0.01" \
  --chemical-policy "adaptive_soft_guard"
