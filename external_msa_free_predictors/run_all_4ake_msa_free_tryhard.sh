#!/usr/bin/env bash
# One command to try every no-MSA/non-AlphaFold predictor source and run the project probe on every produced PDB.
# Usage: ./run_all_4ake_msa_free_tryhard.sh /path/to/pharmacotopology_project
set -euo pipefail
PROJECT=${1:-..}
ROOT=$(cd "$(dirname "$0")" && pwd)
RUN_ROOT=${RUN_ROOT:-$ROOT/tryhard_runs}
TIMEOUT_ESMFOLD_API=${ESMFOLD_API_TIMEOUT_SECONDS:-300}
TIMEOUT_LOCAL_ESMFOLD=${LOCAL_ESMFOLD_TIMEOUT_SECONDS:-900}
mkdir -p "$RUN_ROOT"
LOG="$RUN_ROOT/run_all_4ake_msa_free_tryhard.log"
: > "$LOG"
log(){ echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*" | tee -a "$LOG"; }

log "PROJECT=$PROJECT"
log "RUN_ROOT=$RUN_ROOT"
log "NO GIF. NO ALPHAFOLD. NO MSA. HARD TIMEOUTS."

# 1. ESM Atlas / ESMFold API, raw text/plain.
log "Trying ESMFold API raw text/plain"
(
  cd "$ROOT"
  timeout "$TIMEOUT_ESMFOLD_API" ./esmfold_api.sh 4ake.fasta "$RUN_ROOT/4ake_esmfold_api.pdb"
) >> "$LOG" 2>&1 || log "ESMFold API did not produce a usable PDB in this environment"

# 2. Local ESMFold v1 if dependencies/weights exist.
log "Trying local ESMFold v1"
(
  cd "$ROOT"
  timeout "$TIMEOUT_LOCAL_ESMFOLD" ./run_local_esmfold_v1.py 4ake.fasta "$RUN_ROOT/4ake_esmfold_local.pdb"
) >> "$LOG" 2>&1 || log "Local ESMFold v1 unavailable or failed"

# 3. OmegaFold fallback.
log "Trying OmegaFold safe wrapper"
(
  cd "$ROOT"
  ./run_omegafold_safe.sh 4ake.fasta "$RUN_ROOT/omegafold_out"
) >> "$LOG" 2>&1 || log "OmegaFold unavailable or failed"

# 4. Accept externally dropped PDBs too.
log "Scanning for candidate PDB files"
python3 "$ROOT/summarize_candidate_pdbs.py" "$RUN_ROOT" | tee -a "$LOG"

# 5. Run probe on all non-empty candidate PDBs.
log "Running pharmacotopology probe on every candidate PDB"
python3 "$ROOT/run_probe_for_candidate_pdbs.py" \
  --project "$PROJECT" \
  --candidate-root "$RUN_ROOT" \
  --out-root "$RUN_ROOT/probe_results" \
  --min-votes 1 | tee -a "$LOG"

log "Aggregating verdict"
python3 "$ROOT/summarize_probe_reports.py" "$RUN_ROOT/probe_results" | tee "$RUN_ROOT/final_tryhard_verdict.json" | tee -a "$LOG"
log "DONE"
