# V15 lock + 4AKE bridge commands

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
  tests/test_v15_dynamic_role_grammar_panel_lock.py
```

Run V15 if needed:

```bash
bash "$REPO_ROOT/scripts/run_v15_dynamic_separation_grammar_readout.sh"
python3 "$REPO_ROOT/scripts/print_v15_dynamic_separation_grammar_summary.py"
```

Lock V15 partial panel:

```bash
bash "$REPO_ROOT/scripts/run_v15_dynamic_role_grammar_panel_lock.sh"
python3 "$REPO_ROOT/scripts/print_v15_dynamic_role_grammar_panel_lock_summary.py"
```

Build 4AKE bridge:

```bash
bash "$REPO_ROOT/scripts/run_v15_4ake_dynamic_grammar_bridge.sh"
python3 "$REPO_ROOT/scripts/print_v15_4ake_dynamic_grammar_bridge_summary.py"
```

Rerun V15 after bridge:

```bash
bash "$REPO_ROOT/scripts/run_v15_dynamic_separation_grammar_readout.sh"
python3 "$REPO_ROOT/scripts/print_v15_dynamic_separation_grammar_summary.py"
```
