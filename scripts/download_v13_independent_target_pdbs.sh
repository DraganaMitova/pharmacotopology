#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"

RAW_DIR="$REPO_ROOT/data/independent_contact_sources/raw"
OUT_DIR="$REPO_ROOT/data/independent_contact_sources"
mkdir -p "$RAW_DIR" "$OUT_DIR"

# Strict honest target preparation for V13 runtime tests.
#
# This script does NOT synthesize targets, does NOT use native/debug fallback,
# and does NOT accept partial-length targets. A target is accepted only if the
# downloaded real structure contains the full benchmark sequence and can be
# trimmed/renumbered exactly.
#
# 1UBQ transfer target policy:
#   Use independent real ubiquitin structures with full 76-residue sequence.
#   Defaults avoid 1UBI because its RCSB title/provenance is synthetic.
#   Defaults also avoid 3ONS because it commonly resolves only 72/76 residues.
#   First candidates: 1D3Z (human ubiquitin NMR), 5DK8 (human ubiquitin crystal).
#
# 1CLL transfer target policy:
#   AFDB URLs for common CALM ids can return 404, so strict mode now tries
#   non-native real RCSB calmodulin targets first. Native 1CLL is explicitly
#   forbidden as fallback. Each candidate must contain the exact full benchmark
#   148-aa sequence before it is accepted.
#   First candidates: 1CFC, 1DMO, 1WRZ, 1YR5, 3IF7, 1LIN, 1CTR.
#
# Useful commands:
#   ONLY=1UBQ bash scripts/download_v13_independent_target_pdbs.sh
#   ONLY=1CLL bash scripts/download_v13_independent_target_pdbs.sh
#   UBQ_PDB_IDS="1D3Z 5DK8" bash scripts/download_v13_independent_target_pdbs.sh
#   CLL_RCSB_IDS="1CFC 1DMO 1WRZ 1YR5" bash scripts/download_v13_independent_target_pdbs.sh

ONLY="${ONLY:-all}"
if [[ "$ONLY" != "all" && "$ONLY" != "1UBQ" && "$ONLY" != "1CLL" ]]; then
  echo "invalid ONLY=$ONLY ; use all, 1UBQ, or 1CLL" >&2
  exit 2
fi

fetch_rcsb() {
  local label="$1"
  local id="$2"
  local out url
  id="$(printf '%s' "$id" | tr '[:lower:]' '[:upper:]')"
  if [[ "$id" == "1CLL" ]]; then
    echo "$label refused native RCSB fallback candidate: $id" >&2
    return 1
  fi
  out="$RAW_DIR/RCSB-${id}.pdb"
  url="https://files.rcsb.org/download/${id}.pdb"
  if [[ -s "$out" ]]; then
    echo "$label already present: $out" >&2
    printf '%s\n' "$out"
    return 0
  fi
  echo "$label trying RCSB: $url" >&2
  if curl -fL "$url" -o "$out"; then
    printf '%s\n' "$out"
    return 0
  fi
  rm -f "$out"
  return 1
}

fetch_af() {
  local label="$1"
  local id="$2"
  local out url
  out="$RAW_DIR/AF-${id}-F1-model_v4.pdb"
  url="https://alphafold.ebi.ac.uk/files/AF-${id}-F1-model_v4.pdb"
  if [[ -s "$out" ]]; then
    echo "$label already present: $out" >&2
    printf '%s\n' "$out"
    return 0
  fi
  echo "$label trying AFDB: $url" >&2
  if curl -fL "$url" -o "$out"; then
    printf '%s\n' "$out"
    return 0
  fi
  rm -f "$out"
  return 1
}

export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts"

prepare_ubq_from_raw_exact() {
  local raw="$1"
  python3 "$REPO_ROOT/scripts/prepare_v13_independent_target_pdbs.py" \
    --repo-root "$REPO_ROOT" \
    --only 1UBQ \
    --ubq-raw-pdb "${raw#$REPO_ROOT/}" \
    --ubq-out-pdb "data/independent_contact_sources/1UBQ_independent_target_segment.pdb"
}

prepare_cll_from_raw_exact() {
  local raw="$1"
  python3 "$REPO_ROOT/scripts/prepare_v13_independent_target_pdbs.py" \
    --repo-root "$REPO_ROOT" \
    --only 1CLL \
    --cll-raw-pdb "${raw#$REPO_ROOT/}" \
    --cll-out-pdb "data/independent_contact_sources/1CLL_independent_target_segment.pdb"
}

if [[ "$ONLY" == "all" || "$ONLY" == "1UBQ" ]]; then
  rm -f "$OUT_DIR/1UBQ_independent_target_segment.pdb" "$OUT_DIR/1UBQ_independent_target_segment.pdb.provenance.json"
  UBQ_PREPARED=0
  read -r -a UBQ_IDS <<< "${UBQ_PDB_IDS:-1D3Z 5DK8}"
  for id in "${UBQ_IDS[@]}"; do
    raw="$(fetch_rcsb "1UBQ independent real target" "$id")" || continue
    echo "1UBQ preparing exact full-length target from RCSB $id" >&2
    if prepare_ubq_from_raw_exact "$raw"; then
      UBQ_PREPARED=1
      echo "1UBQ exact independent target accepted from RCSB $id" >&2
      break
    fi
    echo "1UBQ exact preparation failed for RCSB $id; trying next exact candidate" >&2
  done
  if [[ "$UBQ_PREPARED" != "1" ]]; then
    echo "failed to prepare strict 1UBQ independent target from candidates: ${UBQ_IDS[*]}" >&2
    echo "No partial/native/synthetic fallback was used. Add another real full-length candidate via UBQ_PDB_IDS." >&2
    exit 1
  fi
fi

if [[ "$ONLY" == "all" || "$ONLY" == "1CLL" ]]; then
  rm -f "$OUT_DIR/1CLL_independent_target_segment.pdb" "$OUT_DIR/1CLL_independent_target_segment.pdb.provenance.json"
  CLL_PREPARED=0

  # Prefer non-native real RCSB structures. All are still strict exact-sequence gated by prepare_v13_independent_target_pdbs.py.
  read -r -a CLL_RCSB_IDS_ARRAY <<< "${CLL_RCSB_IDS:-1CFC 1DMO 1WRZ 1YR5 3IF7 1LIN 1CTR}"
  for id in "${CLL_RCSB_IDS_ARRAY[@]}"; do
    raw="$(fetch_rcsb "1CLL independent non-native real target" "$id")" || continue
    echo "1CLL preparing exact full-length target from non-native RCSB $id" >&2
    if prepare_cll_from_raw_exact "$raw"; then
      CLL_PREPARED=1
      echo "1CLL exact independent target accepted from non-native RCSB $id" >&2
      break
    fi
    echo "1CLL exact preparation failed for RCSB $id; trying next exact candidate" >&2
  done

  # Optional AFDB fallback remains strict, but is not relied on because common current URLs may 404.
  if [[ "$CLL_PREPARED" != "1" ]]; then
    read -r -a CLL_AF_IDS_ARRAY <<< "${CLL_AF_IDS:-P0DP23 P62158}"
    for id in "${CLL_AF_IDS_ARRAY[@]}"; do
      raw="$(fetch_af "1CLL independent AFDB target" "$id")" || continue
      echo "1CLL preparing exact target from AFDB $id" >&2
      if prepare_cll_from_raw_exact "$raw"; then
        CLL_PREPARED=1
        echo "1CLL exact independent target accepted from AFDB $id" >&2
        break
      fi
      echo "1CLL exact preparation failed for AFDB $id; trying next exact candidate" >&2
    done
  fi

  if [[ "$CLL_PREPARED" != "1" ]]; then
    echo "failed to prepare strict 1CLL independent target." >&2
    echo "Tried non-native RCSB candidates: ${CLL_RCSB_IDS_ARRAY[*]}" >&2
    echo "Tried AFDB candidates: ${CLL_AF_IDS:-P0DP23 P62158}" >&2
    echo "No native RCSB 1CLL fallback, partial target, synthetic target, or generated target was used." >&2
    exit 1
  fi
fi

echo "ready strict target PDBs:"
ls -lh "$OUT_DIR"/*independent_target_segment.pdb
