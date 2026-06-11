# V33 Negative Controls Commands

Run this after V33 has passed with KcsA selected.

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 scripts/build_kcsa_rcsb_external_constraint_import_v0.py
python3 scripts/run_v32_external_constraint_source_import_preflight_v0.py
python3 scripts/run_v33_constraint_backed_operator_readout_v0.py
python3 scripts/run_v33_negative_controls_v0.py
python3 scripts/print_v33_negative_controls.py
python3 -m pytest -q tests/test_v33_constraint_backed_operator_readout.py tests/test_v33_negative_controls.py
```

Expected:

```text
negative_control_status: V33_NEGATIVE_CONTROLS_PASSED_CLAIM_DISABLED
passed_control_count: 9 / 9
claim_allowed: False
new_MD_allowed: False
positive_folding_evidence_found: False
folding_problem_solved: False
```

Meaning: the V32/V33 evidence grammar blocks internal-runtime source poisoning, annotation-only promotion, missing KcsA buckets, and wrong-target readout. This still does not prove de novo folding.
