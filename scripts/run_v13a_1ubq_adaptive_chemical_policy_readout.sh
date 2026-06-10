#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts${PYTHONPATH:+:$PYTHONPATH}"

SOURCE_RUN_DIR="${SOURCE_RUN_DIR:-$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_REPAIR_FIXED}"
OUT_DIR="${OUT_DIR:-$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_ADAPTIVE_CHEMICAL_POLICY_READOUT}"
FREQUENCY_GRID="${FREQUENCY_GRID:-0.50:0.98:0.01}"
CHEMICAL_GRID="${CHEMICAL_GRID:-0.00:0.50:0.05}"
VOTE_THRESHOLD="${VOTE_THRESHOLD:-7}"
BALANCED_DCA_THRESHOLD="${BALANCED_DCA_THRESHOLD:-0.80}"
LEGACY_CHEMICAL_REFERENCE_THRESHOLD="${LEGACY_CHEMICAL_REFERENCE_THRESHOLD:-0.50}"

python3 "$REPO_ROOT/scripts/run_v13a_purpose_gate_readout_v0.py" \
  --source-run-dir "$SOURCE_RUN_DIR" \
  --out-dir "$OUT_DIR" \
  --source-accession "1UBQ:A" \
  --benchmark-file "$REPO_ROOT/data/folding_real_coordinate_holdout_1ubq.locked.json" \
  --external-coupling-file "$REPO_ROOT/data/folding_real_coordinate_holdout_1ubq_external_couplings.v0.locked.json" \
  --anchor-profile-file "$REPO_ROOT/data/folding_real_coordinate_holdout_1ubq_external_couplings.v0.locked.json" \
  --audit-pairs-json "$REPO_ROOT/data/audit_1ubq_debug.json" \
  --topology-mode "none" \
  --domain-boundaries "" \
  --frequency-grid "$FREQUENCY_GRID" \
  --chemical-grid "$CHEMICAL_GRID" \
  --chemical-policy "adaptive_soft_guard" \
  --legacy-chemical-reference-threshold "$LEGACY_CHEMICAL_REFERENCE_THRESHOLD" \
  --vote-threshold "$VOTE_THRESHOLD" \
  --balanced-dca-threshold "$BALANCED_DCA_THRESHOLD"
