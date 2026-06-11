#!/usr/bin/env bash
set -euo pipefail

: "${REPO_ROOT:=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"
BASE="$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run"

must_exist () {
  if [ ! -f "$1" ]; then
    echo "STOP: missing required artifact:" >&2
    echo "$1" >&2
    exit 1
  fi
}

# Recover only small locked source summaries. This does NOT restore raw trajectories.
python3 "$REPO_ROOT/scripts/recover_v15_locked_source_snapshots_v0.py"

# Rebuild V15 from recovered source summaries.
bash "$REPO_ROOT/scripts/run_v15_4ake_dynamic_grammar_bridge.sh" \
  --role-cert "$BASE/V15_4AKE_BALANCED_CANDIDATE_READOUT/v15_4ake_balanced_candidate_readout_certificate.json"
bash "$REPO_ROOT/scripts/run_v15_dynamic_separation_grammar_readout.sh"
bash "$REPO_ROOT/scripts/run_v15_dynamic_role_grammar_panel_lock.sh"
must_exist "$BASE/V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCKED/v15_dynamic_role_grammar_panel_locked_certificate.json"

# V16-V30 chain.
bash "$REPO_ROOT/scripts/download_v16_pressure_target_structures.sh"
bash "$REPO_ROOT/scripts/run_v16_target_manifest_and_role_expectation_lock.sh"
bash "$REPO_ROOT/scripts/run_v16_transfer_data_preflight.sh"
bash "$REPO_ROOT/scripts/run_v16_zero_md_role_transfer_readout.sh"
bash "$REPO_ROOT/scripts/run_v16_pressure_evidence_gap_lock.sh"
bash "$REPO_ROOT/scripts/run_v17_pressure_evidence_sprint_lock.sh"
bash "$REPO_ROOT/scripts/run_v18_p53_tad_mdm2_partner_induced_evidence_test.sh"
bash "$REPO_ROOT/scripts/run_v18_p53_partner_induced_evidence_lock.sh"
bash "$REPO_ROOT/scripts/run_v18b_p53_isolated_tad_abstain_context.sh"
bash "$REPO_ROOT/scripts/run_v19_kcsa_membrane_pore_evidence_readout.sh"
bash "$REPO_ROOT/scripts/run_v19_kcsa_membrane_pore_evidence_lock.sh"
bash "$REPO_ROOT/scripts/run_v20_xcl1_state_specific_evidence_readout.sh"
bash "$REPO_ROOT/scripts/run_v21_pressure_evidence_panel_summary.sh"
bash "$REPO_ROOT/scripts/run_v22_external_evidence_and_annotation_acquisition_panel.sh"
bash "$REPO_ROOT/scripts/run_v23_p53_tad_mdm2_external_evidence_contrast_test.sh"
bash "$REPO_ROOT/scripts/run_v24_kcsa_external_annotation_and_assembly_acquisition.sh"
bash "$REPO_ROOT/scripts/run_v25_fast_mechanism_evidence_sprint.sh"
bash "$REPO_ROOT/scripts/run_v26_xcl1_state_separation_operator_test.sh"
bash "$REPO_ROOT/scripts/run_v27_xcl1_condition_and_coupling_evidence_acquisition.sh"
bash "$REPO_ROOT/scripts/run_v28_xcl1_state_condition_evidence_contrast_test.sh"
bash "$REPO_ROOT/scripts/run_v29_mechanism_operator_panel_summary_and_md_readiness_decision.sh"
bash "$REPO_ROOT/scripts/run_v30_external_constraint_and_coupling_acquisition_sprint.sh"

# Export small locked certs for future restore.
python3 "$REPO_ROOT/scripts/export_locked_runtime_certificates_v0.py" || true
python3 "$REPO_ROOT/scripts/print_v30_external_constraint_and_coupling_acquisition_sprint.py"
