#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"

RAW_DIR="$REPO_ROOT/data/independent_contact_sources/raw"
OUT_DIR="$REPO_ROOT/data/independent_contact_sources"
mkdir -p "$RAW_DIR" "$OUT_DIR"

# 1UBQ PDB metadata names P62988, but AlphaFold DB may not expose AF-P62988.
# SIFTS/RCSB maps the same mature ubiquitin sequence to polyubiquitin accessions.
# Default order: use a real AFDB polyubiquitin model and trim/renumber the exact 76-aa segment.
if [[ -n "${UBQ_AF_ID:-}" ]]; then
  UBQ_CANDIDATES=("$UBQ_AF_ID")
else
  UBQ_CANDIDATES=("P0CG47" "P0CG48" "P62988")
fi

if [[ -n "${CLL_AF_ID:-}" ]]; then
  CLL_CANDIDATES=("$CLL_AF_ID")
else
  CLL_CANDIDATES=("P0DP23")
fi

fetch_first_available() {
  local label="$1"
  shift
  local candidates=("$@")
  local id out url
  for id in "${candidates[@]}"; do
    out="$RAW_DIR/AF-${id}-F1-model_v4.pdb"
    url="https://alphafold.ebi.ac.uk/files/AF-${id}-F1-model_v4.pdb"
    if [[ -s "$out" ]]; then
      echo "$label already present: $out" >&2
      printf '%s\n' "$id"
      return 0
    fi
    echo "$label trying: $url" >&2
    if curl -fL "$url" -o "$out"; then
      printf '%s\n' "$id"
      return 0
    fi
    rm -f "$out"
  done
  echo "failed to download $label AlphaFold target from candidates: ${candidates[*]}" >&2
  return 1
}

UBQ_USED_ID="$(fetch_first_available "1UBQ" "${UBQ_CANDIDATES[@]}")"
CLL_USED_ID="$(fetch_first_available "1CLL" "${CLL_CANDIDATES[@]}")"

echo "1UBQ raw AF id used: $UBQ_USED_ID"
echo "1CLL raw AF id used: $CLL_USED_ID"

export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts"
python3 "$REPO_ROOT/scripts/prepare_v13_independent_target_pdbs.py" \
  --repo-root "$REPO_ROOT" \
  --ubq-raw-pdb "data/independent_contact_sources/raw/AF-${UBQ_USED_ID}-F1-model_v4.pdb" \
  --ubq-out-pdb "data/independent_contact_sources/1UBQ_independent_target_segment.pdb" \
  --cll-raw-pdb "data/independent_contact_sources/raw/AF-${CLL_USED_ID}-F1-model_v4.pdb" \
  --cll-out-pdb "data/independent_contact_sources/1CLL_independent_target_segment.pdb"

echo "ready target PDBs:"
ls -lh \
  "$OUT_DIR/1UBQ_independent_target_segment.pdb" \
  "$OUT_DIR/1CLL_independent_target_segment.pdb"
