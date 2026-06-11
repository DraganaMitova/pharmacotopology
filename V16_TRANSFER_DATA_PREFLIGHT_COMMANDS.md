# V16 Transfer Data Preflight Commands

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"
```

## 1. Tests

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
  tests/test_v16_transfer_data_preflight.py
```

## 2. Run preflight before downloading material

```bash
bash "$REPO_ROOT/scripts/run_v16_transfer_data_preflight.sh"
python3 "$REPO_ROOT/scripts/print_v16_transfer_data_preflight_summary.py"
```

This should honestly report missing material if no public target structures have been downloaded yet.

## 3. Download public target/context structures

```bash
bash "$REPO_ROOT/scripts/download_v16_pressure_target_structures.sh"
```

This downloads public RCSB PDB files and writes provenance JSON files. It does not run MD.

## 4. Re-run data preflight

```bash
bash "$REPO_ROOT/scripts/run_v16_transfer_data_preflight.sh"
python3 "$REPO_ROOT/scripts/print_v16_transfer_data_preflight_summary.py"
```

Expected after download:

```text
V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT
```

## Boundary

This is not evidence selection yet. It only tells whether V16 has enough clean material to run the next zero-MD role-transfer readout.
