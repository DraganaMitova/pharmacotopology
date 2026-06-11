#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
RUN_ROOT="$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run"

# Keep this sprint fast: only refresh upstream zero-MD artifacts if they are missing.
if [ ! -f "$RUN_ROOT/V15_DYNAMIC_SEPARATION_GRAMMAR_READOUT/v15_dynamic_separation_grammar_readout_certificate.json" ]; then
  bash "$REPO_ROOT/scripts/run_v15_dynamic_separation_grammar_readout.sh"
fi

if [ ! -f "$RUN_ROOT/V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST/v23_p53_tad_mdm2_external_evidence_contrast_certificate.json" ]; then
  bash "$REPO_ROOT/scripts/run_v23_p53_tad_mdm2_external_evidence_contrast_test.sh"
fi

if [ ! -f "$RUN_ROOT/V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION/v24_kcsa_external_annotation_and_assembly_acquisition_certificate.json" ]; then
  bash "$REPO_ROOT/scripts/run_v24_kcsa_external_annotation_and_assembly_acquisition.sh"
fi

if [ ! -f "$RUN_ROOT/V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT/v20_xcl1_state_specific_evidence_readout_certificate.json" ]; then
  bash "$REPO_ROOT/scripts/run_v20_xcl1_state_specific_evidence_readout.sh"
fi

python3 "$REPO_ROOT/scripts/run_v25_fast_mechanism_evidence_sprint_v0.py"
