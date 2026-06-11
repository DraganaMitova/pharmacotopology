# V34 KcsA Discriminative Content Controls Commands

Run this after V33 and V33 negative controls.

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 scripts/build_kcsa_rcsb_external_constraint_import_v0.py
python3 scripts/run_v32_external_constraint_source_import_preflight_v0.py
python3 scripts/run_v33_constraint_backed_operator_readout_v0.py
python3 scripts/run_v33_negative_controls_v0.py
python3 scripts/run_v34_kcsa_discriminative_content_controls_v0.py
python3 scripts/print_v34_kcsa_discriminative_content_controls.py
python3 -m pytest -q tests/test_v33_constraint_backed_operator_readout.py tests/test_v33_negative_controls.py tests/test_v34_kcsa_discriminative_content_controls.py
```

Expected:

```text
control_status: V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_PASSED_CLAIM_DISABLED
passed_control_count: 8 / 8
claim_allowed: False
new_MD_allowed: False
positive_folding_evidence_found: False
folding_problem_solved: False
```

Meaning: V33 is no longer only a file-count/readout check. The imported KcsA source must contain discriminative pore/filter identity, K+ identity, TVGYG identity, and tetramer/interface identity. Damaged versions must fail. Still no de-novo folding claim.
