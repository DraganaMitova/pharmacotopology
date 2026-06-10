#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"

RAW_DIR="$REPO_ROOT/data/independent_contact_sources/raw"
OUT_DIR="$REPO_ROOT/data/independent_contact_sources"
mkdir -p "$RAW_DIR" "$OUT_DIR"

# This script prepares target PDBs needed by the V13 runtime tests.
# It intentionally does NOT tune role rules and every prepared target writes
# provenance with claim_allowed=false.
#
# Important 1UBQ note:
# AFDB currently returns 404 for the usual mature ubiquitin accessions in this setup.
# So V13a defaults to an independent, highly similar experimental ubiquitin structure
# from RCSB (3ONS first, then 1UBI). This is a target/provenance repair, not a claim.
#
# Useful overrides:
#   ONLY=1UBQ bash scripts/download_v13_independent_target_pdbs.sh
#   UBQ_PDB_ID=3ONS bash scripts/download_v13_independent_target_pdbs.sh
#   ALLOW_NATIVE_TARGET_DEBUG=1 bash scripts/download_v13_independent_target_pdbs.sh
#
# ALLOW_NATIVE_TARGET_DEBUG=1 may fall back to bundled 1UBQ coordinates only to test
# the runtime plumbing. That mode is NOT an independent biological transfer result.

ONLY="${ONLY:-all}"
if [[ "$ONLY" != "all" && "$ONLY" != "1UBQ" && "$ONLY" != "1CLL" ]]; then
  echo "invalid ONLY=$ONLY ; use all, 1UBQ, or 1CLL" >&2
  exit 2
fi

fetch_rcsb_first_available() {
  local label="$1"
  shift
  local candidates=("$@")
  local id out url
  for id in "${candidates[@]}"; do
    id="$(printf '%s' "$id" | tr '[:lower:]' '[:upper:]')"
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
  done
  return 1
}

fetch_af_first_available() {
  local label="$1"
  shift
  local candidates=("$@")
  local id out url
  for id in "${candidates[@]}"; do
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
  done
  return 1
}

UBQ_RAW=""
if [[ "$ONLY" == "all" || "$ONLY" == "1UBQ" ]]; then
  if [[ -n "${UBQ_PDB_ID:-}" ]]; then
    UBQ_RAW="$(fetch_rcsb_first_available "1UBQ-like target" "$UBQ_PDB_ID")" || {
      echo "failed to download requested UBQ_PDB_ID=$UBQ_PDB_ID" >&2
      exit 1
    }
  else
    UBQ_RAW="$(fetch_rcsb_first_available "1UBQ-like target" 3ONS 1UBI)" || true
  fi

  if [[ -z "$UBQ_RAW" && "${ALLOW_NATIVE_TARGET_DEBUG:-0}" == "1" ]]; then
    if [[ -s "$REPO_ROOT/data/rcsb_pdb/1UBQ.pdb" ]]; then
      UBQ_RAW="$REPO_ROOT/data/rcsb_pdb/1UBQ.pdb"
      echo "1UBQ DEBUG fallback using bundled native PDB: $UBQ_RAW" >&2
      echo "WARNING: this is runtime plumbing only, not independent transfer evidence." >&2
    else
      UBQ_RAW="$(fetch_rcsb_first_available "1UBQ native debug target" 1UBQ)" || true
      echo "WARNING: downloaded 1UBQ native debug target; not independent transfer evidence." >&2
    fi
  fi

  if [[ -z "$UBQ_RAW" ]]; then
    echo "failed to prepare 1UBQ-like target. Try: UBQ_PDB_ID=3ONS or ALLOW_NATIVE_TARGET_DEBUG=1" >&2
    exit 1
  fi
fi

CLL_RAW=""
if [[ "$ONLY" == "all" || "$ONLY" == "1CLL" ]]; then
  if [[ -n "${CLL_PDB_ID:-}" ]]; then
    CLL_RAW="$(fetch_rcsb_first_available "1CLL target" "$CLL_PDB_ID")" || {
      echo "failed to download requested CLL_PDB_ID=$CLL_PDB_ID" >&2
      exit 1
    }
  else
    CLL_RAW="$(fetch_af_first_available "1CLL target" P0DP23 P62158)" || true
    if [[ -z "$CLL_RAW" ]]; then
      CLL_RAW="$(fetch_rcsb_first_available "1CLL target" 1CLL)" || true
    fi
  fi
  if [[ -z "$CLL_RAW" ]]; then
    echo "failed to prepare 1CLL target from AFDB/RCSB candidates" >&2
    exit 1
  fi
fi

export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts"

if [[ "$ONLY" == "1UBQ" ]]; then
  python3 "$REPO_ROOT/scripts/prepare_v13_independent_target_pdbs.py" \
    --repo-root "$REPO_ROOT" \
    --only 1UBQ \
    --ubq-raw-pdb "${UBQ_RAW#$REPO_ROOT/}" \
    --ubq-out-pdb "data/independent_contact_sources/1UBQ_independent_target_segment.pdb"
elif [[ "$ONLY" == "1CLL" ]]; then
  python3 "$REPO_ROOT/scripts/prepare_v13_independent_target_pdbs.py" \
    --repo-root "$REPO_ROOT" \
    --only 1CLL \
    --cll-raw-pdb "${CLL_RAW#$REPO_ROOT/}" \
    --cll-out-pdb "data/independent_contact_sources/1CLL_independent_target_segment.pdb"
else
  python3 "$REPO_ROOT/scripts/prepare_v13_independent_target_pdbs.py" \
    --repo-root "$REPO_ROOT" \
    --only all \
    --ubq-raw-pdb "${UBQ_RAW#$REPO_ROOT/}" \
    --ubq-out-pdb "data/independent_contact_sources/1UBQ_independent_target_segment.pdb" \
    --cll-raw-pdb "${CLL_RAW#$REPO_ROOT/}" \
    --cll-out-pdb "data/independent_contact_sources/1CLL_independent_target_segment.pdb"
fi

echo "ready target PDBs:"
ls -lh "$OUT_DIR"/*independent_target_segment.pdb
