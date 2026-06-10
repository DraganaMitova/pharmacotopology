#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts"

TARGET_PDB="$REPO_ROOT/data/independent_contact_sources/AF-P0DP23-F1-model_v4_1CLL_calmodulin_segment.pdb"
OUT_DIR="$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13b_1CLL_CALMODULIN_FIXED"

if [[ ! -f "$TARGET_PDB" ]]; then
  echo "missing prepared 1CLL target PDB: $TARGET_PDB" >&2
  echo "first run: bash \"$REPO_ROOT/scripts/download_v13_independent_target_pdbs.sh\"" >&2
  exit 2
fi

if [[ -d "$OUT_DIR" ]]; then
  if [[ "${CLEAN:-0}" != "1" ]]; then
    echo "output dir already exists: $OUT_DIR" >&2
    echo "rerun with CLEAN=1 to replace it:" >&2
    echo "  CLEAN=1 bash \"$REPO_ROOT/scripts/run_v13b_1cll_calmodulin_fixed.sh\"" >&2
    exit 2
  fi
  rm -rf "$OUT_DIR"
fi

COMMON_ARGS=(
  --source-accession "1CLL:A"
  --benchmark-file "$REPO_ROOT/data/folding_real_coordinate_visual_8.locked.json"
  --external-coupling-file "$REPO_ROOT/data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
  --anchor-profile-file "$REPO_ROOT/data/folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
  --audit-pairs-json "$REPO_ROOT/data/audit_1cll_v13b.json"
  --target-pdb "$TARGET_PDB"
  --out-dir "$OUT_DIR"
  --replicas "${REPLICAS:-10}"
  --max-parallel "${MAX_PARALLEL:-2}"
  --topology-mode "interdomain"
  --domain-boundaries "1-75,80-148"
  --runtime-v10-role-aware-selector
  --platform "${OPENMM_PLATFORM:-cpu}"
)

python3 "$REPO_ROOT/scripts/run_openmm_tmd_replicas_v0.py" "${COMMON_ARGS[@]}" --input-preflight-only
python3 "$REPO_ROOT/scripts/run_openmm_tmd_replicas_v0.py" "${COMMON_ARGS[@]}"
