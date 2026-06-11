# V31 Constraint-Backed Operator Readout Preflight

Purpose: separate true external constraint/coupling evidence from annotation-only context and generated internal reports before any V32 constraint-backed readout.

This panel does not run MD, does not tune thresholds, and does not allow a folding claim.

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 -m compileall -q src scripts
python3 -m pytest -q tests/test_v31_constraint_backed_operator_readout_preflight.py

bash "$REPO_ROOT/scripts/run_v31_constraint_backed_operator_readout_preflight.sh"
python3 "$REPO_ROOT/scripts/print_v31_constraint_backed_operator_readout_preflight.py"
python3 "$REPO_ROOT/scripts/export_locked_runtime_certificates_v0.py"
```

Expected honest outcomes:

- If XCL1 or KcsA has real external coupling/constraint files: V31 may select an automatic V32 target.
- If XCL1/KcsA only have annotations or generated internal reports: V31 clean-abstains and blocks constraint-backed V32.

Core policy:

```text
real_external_constraint_or_coupling -> allowed_for_constraint_backed_operator_readout
real_external_alignment_source -> allowed_for_constraint_derivation_preflight_only
annotation_only_external_context -> allowed_for_role_context_only
external_structure_source -> allowed_for_structure_context_or_validation_only
generated_internal_report -> allowed_for_audit_only
unusable/unverified -> excluded
```

No generated runtime report is allowed to count as external evidence.
