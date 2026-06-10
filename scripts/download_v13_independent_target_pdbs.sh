#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"

RAW_DIR="$REPO_ROOT/data/independent_contact_sources/raw"
mkdir -p "$RAW_DIR"

UBQ_AF_ID="${UBQ_AF_ID:-P62988}"
CLL_AF_ID="${CLL_AF_ID:-P0DP23}"

fetch_af() {
  local af_id="$1"
  local out="$2"
  local url="https://alphafold.ebi.ac.uk/files/AF-${af_id}-F1-model_v4.pdb"
  if [[ -s "$out" ]]; then
    echo "already present: $out"
    return 0
  fi
  echo "downloading: $url"
  curl -fL "$url" -o "$out"
}

fetch_af "$UBQ_AF_ID" "$RAW_DIR/AF-${UBQ_AF_ID}-F1-model_v4.pdb"
fetch_af "$CLL_AF_ID" "$RAW_DIR/AF-${CLL_AF_ID}-F1-model_v4.pdb"

# If custom IDs are used, pass matching raw file paths into the segment preparer.
export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts"
python3 "$REPO_ROOT/scripts/prepare_v13_independent_target_pdbs.py" \
  --repo-root "$REPO_ROOT" \
  --ubq-raw-pdb "data/independent_contact_sources/raw/AF-${UBQ_AF_ID}-F1-model_v4.pdb" \
  --cll-raw-pdb "data/independent_contact_sources/raw/AF-${CLL_AF_ID}-F1-model_v4.pdb"

echo "ready target PDBs:"
ls -lh \
  "$REPO_ROOT/data/independent_contact_sources/AF-${UBQ_AF_ID}-F1-model_v4_1UBQ_mature_segment.pdb" \
  "$REPO_ROOT/data/independent_contact_sources/AF-${CLL_AF_ID}-F1-model_v4_1CLL_calmodulin_segment.pdb" 2>/dev/null || true
