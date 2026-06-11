#!/usr/bin/env bash
set -euo pipefail

# Downloads only public RCSB target/context structures needed for V16 data preflight.
# This script does not run MD, does not tune grammar, and does not permit claims.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${REPO_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}"
TARGET_ROOT="$REPO_ROOT/data/v16_pressure_targets"

mkdir -p "$TARGET_ROOT/p53_TAD_MDM2" "$TARGET_ROOT/KcsA" "$TARGET_ROOT/XCL1_lymphotactin"

write_provenance() {
  local out_json="$1"
  local pdb_id="$2"
  local target_id="$3"
  local role_use="$4"
  local url="https://files.rcsb.org/download/${pdb_id}.pdb"
  python3 - <<PY
import json, datetime
from pathlib import Path
out = Path(r"$out_json")
out.write_text(json.dumps({
  "kind": "V16_PUBLIC_RCSB_STRUCTURE_PROVENANCE_v0",
  "target_id": "$target_id",
  "pdb_id": "$pdb_id",
  "source": "RCSB PDB public download",
  "source_url": "$url",
  "role_use": "$role_use",
  "usage_boundary": "target_context_only_not_native_selection_not_folding_solved_claim",
  "claim_allowed": False,
  "downloaded_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()
}, indent=2, sort_keys=True), encoding="utf-8")
PY
}

download_pdb() {
  local target_id="$1"
  local pdb_id="$2"
  local role_use="$3"
  local out_dir="$TARGET_ROOT/$target_id"
  local out_pdb="$out_dir/${pdb_id}.pdb"
  local out_prov="$out_dir/${pdb_id}.provenance.json"
  local url="https://files.rcsb.org/download/${pdb_id}.pdb"

  mkdir -p "$out_dir"
  echo "Downloading $pdb_id for $target_id -> $out_pdb"
  curl -fL --retry 3 --connect-timeout 20 "$url" -o "$out_pdb"
  write_provenance "$out_prov" "$pdb_id" "$target_id" "$role_use"
}

download_pdb "p53_TAD_MDM2" "1YCR" "partner_induced_interface_or_helix_context"
download_pdb "KcsA" "1K4C" "membrane_pore_oligomer_context"
download_pdb "XCL1_lymphotactin" "2HDM" "state_A_chemokine_monomer_context"
download_pdb "XCL1_lymphotactin" "2JP1" "state_B_beta_sandwich_or_alternative_state_context"

echo "V16 target/context structure download complete."
echo "Next: bash \"$REPO_ROOT/scripts/run_v16_transfer_data_preflight.sh\""
