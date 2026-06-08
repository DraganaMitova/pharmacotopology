#!/usr/bin/env bash
# No-arg default: ./esmfold_api.sh
# Optional: ./esmfold_api.sh input.fasta output.pdb
set -euo pipefail
FASTA=${1:-4ake.fasta}
OUT=${2:-4ake_esmfold_api.pdb}
TIMEOUT=${ESMFOLD_API_TIMEOUT_SECONDS:-300}
python3 "$(dirname "$0")/esmfold_api.py" "$FASTA" "$OUT" --timeout-seconds "$TIMEOUT"
