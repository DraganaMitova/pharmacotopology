# V15 Three-Protein Panel Run Commands

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
  tests/test_v15_4ake_balanced_candidate_readout.py
```

Run 4AKE balanced candidate readout:

```bash
bash "$REPO_ROOT/scripts/run_v15_4ake_balanced_candidate_readout.sh"
python3 "$REPO_ROOT/scripts/print_v15_4ake_balanced_candidate_summary.py"
```

Bridge 4AKE into V15:

```bash
bash "$REPO_ROOT/scripts/run_v15_4ake_dynamic_grammar_bridge.sh" \
  --role-cert "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V15_4AKE_BALANCED_CANDIDATE_READOUT/v15_4ake_balanced_candidate_readout_certificate.json"
python3 "$REPO_ROOT/scripts/print_v15_4ake_dynamic_grammar_bridge_summary.py"
```

Run all-three V15 panel:

```bash
bash "$REPO_ROOT/scripts/run_v15_dynamic_separation_grammar_readout.sh"
python3 "$REPO_ROOT/scripts/print_v15_dynamic_separation_grammar_summary.py"
```

Lock panel:

```bash
bash "$REPO_ROOT/scripts/run_v15_dynamic_role_grammar_panel_lock.sh"
python3 "$REPO_ROOT/scripts/print_v15_dynamic_role_grammar_panel_lock_summary.py"
```
