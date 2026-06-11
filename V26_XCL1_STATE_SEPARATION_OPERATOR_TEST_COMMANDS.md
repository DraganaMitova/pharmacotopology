# V26 XCL1 State-Separation Operator Test Commands

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 -m compileall -q src scripts
python3 -m pytest -q tests/test_v26_xcl1_state_separation_operator_test.py

bash "$REPO_ROOT/scripts/run_v26_xcl1_state_separation_operator_test.sh"
python3 "$REPO_ROOT/scripts/print_v26_xcl1_state_separation_operator_test_summary.py"
```
