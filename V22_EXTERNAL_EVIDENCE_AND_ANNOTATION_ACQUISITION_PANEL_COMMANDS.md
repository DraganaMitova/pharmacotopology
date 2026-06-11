# V22 commands

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"
```

```bash
python3 -m compileall -q src scripts
python3 -m pytest -q \
  tests/test_v13a_purpose_gate_readout.py \
  tests/test_v13b_hierarchical_purpose_topology_readout.py \
  tests/test_v14_unified_protein_esperanto_grammar_readout.py \
  tests/test_v15_dynamic_separation_grammar_readout.py \
  tests/test_v15_4ake_dynamic_grammar_bridge.py \
  tests/test_v15_dynamic_role_grammar_panel_lock.py \
  tests/test_v15_4ake_balanced_candidate_readout.py \
  tests/test_v16_target_manifest_and_role_expectation_lock.py \
  tests/test_v16_transfer_data_preflight.py \
  tests/test_v16_zero_md_role_transfer_readout.py \
  tests/test_v16_pressure_evidence_gap_lock.py \
  tests/test_v17_pressure_evidence_sprint_lock.py \
  tests/test_v18_p53_tad_mdm2_partner_induced_evidence_test.py \
  tests/test_v18_p53_partner_induced_evidence_lock.py \
  tests/test_v18b_p53_isolated_tad_abstain_context.py \
  tests/test_v19_kcsa_membrane_pore_evidence_readout.py \
  tests/test_v19_kcsa_membrane_pore_evidence_lock.py \
  tests/test_v20_xcl1_state_specific_evidence_readout.py \
  tests/test_v21_pressure_evidence_panel_summary.py \
  tests/test_v22_external_evidence_and_annotation_acquisition_panel.py
```

```bash
bash "$REPO_ROOT/scripts/run_v22_external_evidence_and_annotation_acquisition_panel.sh"
python3 "$REPO_ROOT/scripts/print_v22_external_evidence_and_annotation_acquisition_summary.py"
```
