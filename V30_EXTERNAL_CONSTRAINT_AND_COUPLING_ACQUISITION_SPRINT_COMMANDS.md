# V30 External Constraint and Coupling Acquisition Sprint

Purpose: convert the V29 MD-readiness block into a practical constraint/coupling acquisition sprint.

Boundary:
- no MD
- no membrane MD
- no threshold tuning
- no native-metric selection
- no folding claim
- no synthetic couplings

Run:

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 -m compileall -q src scripts
python3 -m pytest -q tests/test_v30_external_constraint_and_coupling_acquisition_sprint.py

bash "$REPO_ROOT/scripts/run_v30_external_constraint_and_coupling_acquisition_sprint.sh"
python3 "$REPO_ROOT/scripts/print_v30_external_constraint_and_coupling_acquisition_sprint.py"
```

Expected high-level status:

```text
sprint_status = V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT_LOCKED
selected_next_panel = V31_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_AND_PREFLIGHT_SPRINT
selected_V31_targets = ['XCL1_lymphotactin', 'KcsA']
new_MD_allowed = false
new_MD_recommended = false
claim_allowed = false
positive_folding_evidence_targets = []
```
