# V29 Mechanism Operator Panel Summary and MD Readiness Decision

V29 summarizes the operator panel after V25/V28 and decides whether any target is ready for MD.
It does not run MD, does not tune thresholds, and does not upgrade pressure evidence into folding evidence.

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 -m compileall -q src scripts
python3 -m pytest -q tests/test_v29_mechanism_operator_panel_summary_and_md_readiness_decision.py

bash "$REPO_ROOT/scripts/run_v29_mechanism_operator_panel_summary_and_md_readiness_decision.sh"
python3 "$REPO_ROOT/scripts/print_v29_mechanism_operator_panel_summary_and_md_readiness_decision.py"
```

Expected high-level result:

```text
summary_status = V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_LOCKED
positive_folding_evidence_targets = []
new_MD_allowed = false
new_MD_recommended = false
selected_next_panel = V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT
```
