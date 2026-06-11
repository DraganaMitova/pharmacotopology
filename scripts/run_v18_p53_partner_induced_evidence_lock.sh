#!/usr/bin/env bash
set -euo pipefail
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
V18_CERT="$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST/v18_p53_tad_mdm2_partner_induced_evidence_certificate.json"
if [ ! -f "$V18_CERT" ]; then
  bash "$REPO_ROOT/scripts/run_v18_p53_tad_mdm2_partner_induced_evidence_test.sh"
fi
python3 "$REPO_ROOT/scripts/run_v18_p53_partner_induced_evidence_lock_v0.py" "$@"
