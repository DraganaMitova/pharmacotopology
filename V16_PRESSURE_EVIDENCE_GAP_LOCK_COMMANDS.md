# V16 Pressure Evidence Gap Lock Commands

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
  tests/test_v16_pressure_evidence_gap_lock.py
```

Re-run V16 material and zero-MD role transfer:

```bash
bash "$REPO_ROOT/scripts/download_v16_pressure_target_structures.sh"
bash "$REPO_ROOT/scripts/run_v16_transfer_data_preflight.sh"
python3 "$REPO_ROOT/scripts/print_v16_transfer_data_preflight_summary.py"

bash "$REPO_ROOT/scripts/run_v16_zero_md_role_transfer_readout.sh"
python3 "$REPO_ROOT/scripts/print_v16_zero_md_role_transfer_summary.py"
```

Run the gap lock:

```bash
bash "$REPO_ROOT/scripts/run_v16_pressure_evidence_gap_lock.sh"
python3 "$REPO_ROOT/scripts/print_v16_pressure_evidence_gap_lock_summary.py"
```

Expected high-level lock:

```text
V16_PRESSURE_EVIDENCE_GAP_LOCKED
role_classification_passed_targets = ['KcsA', 'XCL1_lymphotactin', 'p53_TAD_MDM2']
positive_folding_evidence_targets = []
claim_allowed = false
new_md_executed = false
```
