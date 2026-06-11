# V27 XCL1 Condition/Coupling Acquisition Commands

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 -m compileall -q src scripts
python3 -m pytest -q tests/test_v27_xcl1_condition_and_coupling_evidence_acquisition.py

bash "$REPO_ROOT/scripts/run_v27_xcl1_condition_and_coupling_evidence_acquisition.sh"
python3 "$REPO_ROOT/scripts/print_v27_xcl1_condition_and_coupling_evidence_acquisition_summary.py"
```
