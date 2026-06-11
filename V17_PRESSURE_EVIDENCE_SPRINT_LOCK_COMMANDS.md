# V17 pressure evidence sprint commands

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"
```

Run tests:

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
  tests/test_v17_pressure_evidence_sprint_lock.py
```

Run the one-shot sprint:

```bash
bash "$REPO_ROOT/scripts/run_v17_pressure_evidence_sprint_lock.sh"
python3 "$REPO_ROOT/scripts/print_v17_pressure_evidence_sprint_lock_summary.py"
```

Generated artifacts:

```text
first_contact_clean_pharmacotopology_layer_run/V17_PRESSURE_EVIDENCE_SPRINT_LOCK/v17_pressure_evidence_manifest.json
first_contact_clean_pharmacotopology_layer_run/V17_PRESSURE_EVIDENCE_SPRINT_LOCK/v17_pressure_evidence_preflight.json
first_contact_clean_pharmacotopology_layer_run/V17_PRESSURE_EVIDENCE_SPRINT_LOCK/v17_zero_md_evidence_readout.json
first_contact_clean_pharmacotopology_layer_run/V17_PRESSURE_EVIDENCE_SPRINT_LOCK/v17_next_target_decision.json
first_contact_clean_pharmacotopology_layer_run/V17_PRESSURE_EVIDENCE_SPRINT_LOCK/v17_pressure_evidence_sprint_lock_certificate.json
```
