#!/usr/bin/env bash
set -euo pipefail

DATABASE=""
OUT_DIR="external_msa/4ake_jackhmmer"
QUERY="external_msa_free_predictors/4ake.fasta"
ITER=5
E_VALUE="1e-4"
CPU=2

while [[ $# -gt 0 ]]; do
  case "$1" in
    --database) DATABASE="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --query) QUERY="$2"; shift 2 ;;
    --iter) ITER="$2"; shift 2 ;;
    --evalue) E_VALUE="$2"; shift 2 ;;
    --cpu) CPU="$2"; shift 2 ;;
    -h|--help)
      cat <<'EOF'
Build a raw MSA for 4AKE with jackhmmer.

Required local tools:
  - jackhmmer from HMMER
  - esl-reformat from Easel/HMMER miniapps for aligned FASTA conversion

Required database:
  A local protein FASTA database such as UniProt Swiss-Prot/TrEMBL or UniRef.

Example:
  bash scripts/build_4ake_msa_with_jackhmmer_v0.sh \
    --database /data/uniprot/uniprot_sprot.fasta \
    --out-dir external_msa/4ake_jackhmmer \
    --cpu 4

Outputs:
  4ake_jackhmmer.sto  - Stockholm MSA
  4ake_jackhmmer.afa  - aligned FASTA MSA
  4ake_jackhmmer.tbl  - hit table
  4ake_jackhmmer.log  - command log
EOF
      exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

if [[ -z "$DATABASE" ]]; then
  echo "ERROR: --database /path/to/protein_database.fasta is required" >&2
  exit 2
fi
if [[ ! -f "$QUERY" ]]; then
  echo "ERROR: query FASTA not found: $QUERY" >&2
  exit 2
fi
if [[ ! -f "$DATABASE" ]]; then
  echo "ERROR: database FASTA not found: $DATABASE" >&2
  exit 2
fi
if ! command -v jackhmmer >/dev/null 2>&1; then
  echo "ERROR: jackhmmer not found. Install HMMER first." >&2
  exit 2
fi

mkdir -p "$OUT_DIR"
LOG="$OUT_DIR/4ake_jackhmmer.log"
STO="$OUT_DIR/4ake_jackhmmer.sto"
TBL="$OUT_DIR/4ake_jackhmmer.tbl"
AFA="$OUT_DIR/4ake_jackhmmer.afa"

echo "Running jackhmmer..." | tee "$LOG"
echo "query=$QUERY" | tee -a "$LOG"
echo "database=$DATABASE" | tee -a "$LOG"

jackhmmer \
  --cpu "$CPU" \
  -N "$ITER" \
  -E "$E_VALUE" \
  --incE "$E_VALUE" \
  -A "$STO" \
  --tblout "$TBL" \
  "$QUERY" "$DATABASE" 2>&1 | tee -a "$LOG"

if command -v esl-reformat >/dev/null 2>&1; then
  esl-reformat afa "$STO" > "$AFA"
else
  echo "WARNING: esl-reformat not found; keeping Stockholm only: $STO" | tee -a "$LOG"
fi

if [[ -f "$AFA" ]]; then
  MSA_INPUT="$AFA"
else
  MSA_INPUT="$STO"
fi

echo "Running local MI benchmark with MSA: $MSA_INPUT" | tee -a "$LOG"

PYTHONPATH=src python3 scripts/run_true_local_msa_coevolution_v0.py \
  --source-accession 4AKE:A \
  --msa "$MSA_INPUT" \
  --out-dir first_contact_clean_pharmacotopology_layer_run/true_local_msa_coevolution_v0
