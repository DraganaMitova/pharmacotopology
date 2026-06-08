#!/usr/bin/env bash
# Usage: ./run_omegafold_safe.sh [fasta] [out_dir]
# Tries an installed omegafold first. If absent and network exists, tries source install.
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
run_timeout(){ python3 "$ROOT/run_with_timeout.py" "$@"; }
FASTA=${1:-4ake.fasta}
OUTDIR=${2:-omegafold_out}
TIMEOUT_SECONDS=${OMEGAFOLD_TIMEOUT_SECONDS:-900}
mkdir -p "$OUTDIR"
LOG="$OUTDIR/omegafold_safe.log"
REPORT="$OUTDIR/omegafold_safe_report.json"
{
  echo "kind=omegafold_safe_runner_v1"
  echo "fasta=$FASTA"
  echo "outdir=$OUTDIR"
  echo "timeout_seconds=$TIMEOUT_SECONDS"
  date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
} > "$LOG"
if command -v omegafold >/dev/null 2>&1; then
  echo "using=installed_omegafold" >> "$LOG"
  run_timeout "$TIMEOUT_SECONDS" omegafold "$FASTA" "$OUTDIR" >> "$LOG" 2>&1 || true
elif [ -d OmegaFold ]; then
  echo "using=existing_OmegaFold_checkout" >> "$LOG"
  run_timeout "$TIMEOUT_SECONDS" python3 OmegaFold/main.py "$FASTA" "$OUTDIR" >> "$LOG" 2>&1 || true
else
  echo "using=git_clone_attempt" >> "$LOG"
  run_timeout 180 git clone https://github.com/HeliXonProtein/OmegaFold >> "$LOG" 2>&1 || true
  if [ -d OmegaFold ]; then
    (cd OmegaFold && python3 setup.py install >> "$LOG" 2>&1) || true
    run_timeout "$TIMEOUT_SECONDS" python3 OmegaFold/main.py "$FASTA" "$OUTDIR" >> "$LOG" 2>&1 || true
  fi
fi
python3 - <<PY > "$REPORT"
import json, pathlib, time
out=pathlib.Path("$OUTDIR")
pdbs=sorted(str(p) for p in out.rglob("*.pdb"))
print(json.dumps({"kind":"omegafold_safe_runner_v1","outdir":"$OUTDIR","pdb_count":len(pdbs),"pdbs":pdbs,"status":"success" if pdbs else "no_pdb_produced"}, indent=2))
PY
cat "$REPORT"
