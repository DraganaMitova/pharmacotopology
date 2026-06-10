# V15 Dynamic Separation Grammar Run Commands

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"
```

```bash
python3 - <<'PY'
import run_v15_dynamic_separation_grammar_readout_v0
print("v15_import_ok")
PY

python3 -m compileall -q src scripts
python3 -m pytest -q \
  tests/test_v13a_purpose_gate_readout.py \
  tests/test_v13b_hierarchical_purpose_topology_readout.py \
  tests/test_v14_unified_protein_esperanto_grammar_readout.py \
  tests/test_v15_dynamic_separation_grammar_readout.py
```

```bash
bash "$REPO_ROOT/scripts/run_v15_dynamic_separation_grammar_readout.sh"
python3 "$REPO_ROOT/scripts/print_v15_dynamic_separation_grammar_summary.py"
```
