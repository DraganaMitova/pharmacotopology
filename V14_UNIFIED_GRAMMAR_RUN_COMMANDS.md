# V14 Unified Grammar Run Commands

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 -m compileall -q src scripts
python3 -m pytest -q \
  tests/test_v13a_purpose_gate_readout.py \
  tests/test_v13b_hierarchical_purpose_topology_readout.py \
  tests/test_v14_unified_protein_esperanto_grammar_readout.py

bash "$REPO_ROOT/scripts/run_v14_unified_protein_esperanto_grammar_readout.sh"
python3 "$REPO_ROOT/scripts/print_v14_unified_protein_esperanto_grammar_summary.py"
```

Outputs:

```text
first_contact_clean_pharmacotopology_layer_run/V14_UNIFIED_PROTEIN_ESPERANTO_GRAMMAR_READOUT/
  v14_unified_protein_esperanto_grammar_readout_certificate.json
  v14_unified_protein_esperanto_grammar_table.csv
  V14_UNIFIED_PROTEIN_ESPERANTO_GRAMMAR_REPORT.md
```
