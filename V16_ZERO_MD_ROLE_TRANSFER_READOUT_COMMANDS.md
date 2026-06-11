# V16 Zero-MD Role Transfer Readout Commands

This stage is for pressure-class transfer testing after the V16 target manifest and data preflight are locked.
It does not run MD, does not tune thresholds, does not use native metrics for selection, and does not permit claims.

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
  tests/test_v16_zero_md_role_transfer_readout.py
```

Run or re-run manifest lock:

```bash
bash "$REPO_ROOT/scripts/run_v16_target_manifest_and_role_expectation_lock.sh"
python3 "$REPO_ROOT/scripts/print_v16_target_manifest_and_role_expectation_lock_summary.py"
```

Download public target/context structures:

```bash
bash "$REPO_ROOT/scripts/download_v16_pressure_target_structures.sh"
```

Run V16 data preflight:

```bash
bash "$REPO_ROOT/scripts/run_v16_transfer_data_preflight.sh"
python3 "$REPO_ROOT/scripts/print_v16_transfer_data_preflight_summary.py"
```

Expected after the KcsA chain-context repair and downloaded structures:

```text
V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT
ready_targets = ['p53_TAD_MDM2', 'KcsA', 'XCL1_lymphotactin']
blocked_targets = []
claim_allowed = False
new_md_executed = False
```

Run zero-MD role-transfer readout:

```bash
bash "$REPO_ROOT/scripts/run_v16_zero_md_role_transfer_readout.sh"
python3 "$REPO_ROOT/scripts/print_v16_zero_md_role_transfer_summary.py"
```

Expected interpretation:

```text
p53_TAD_MDM2 = partner-induced complex role context; no autonomous p53 TAD fold claim
KcsA = membrane/pore context detected; no soluble compact-core misclassification; no tetramer/whole-fold claim
XCL1 = two-state metamorphic context detected; no forced single-fold claim
claim_allowed = false
```
