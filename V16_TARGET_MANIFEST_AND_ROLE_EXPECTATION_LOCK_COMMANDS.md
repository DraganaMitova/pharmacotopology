# V16 Target Manifest and Role-Expectation Lock Commands

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"
```

## Tests

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
  tests/test_v16_target_manifest_and_role_expectation_lock.py
```

## Lock V16 manifest

This assumes V15 has already been locked as `V15_THREE_PROTEIN_DYNAMIC_GRAMMAR_PANEL_LOCKED`.

```bash
bash "$REPO_ROOT/scripts/run_v16_target_manifest_and_role_expectation_lock.sh"
python3 "$REPO_ROOT/scripts/print_v16_target_manifest_and_role_expectation_lock_summary.py"
```

Expected:

```text
lock_status = V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCKED
claim_allowed = False
data_preflight_executed = False
new_md_executed = False
```
