# V13 strict honest run commands

This repo variant removes the earlier partial/native/synthetic fallbacks. Official V13 target preparation accepts only a real downloaded target structure containing an exact full-length benchmark sequence segment.

## Compile and unit tests

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts"

python3 -m compileall -q src scripts
PYTHONPATH=src:scripts python3 -m pytest -q
```

Expected package test result in this zip build:

```text
44 passed, 2 skipped
```

## Prepare strict 1UBQ transfer target

```bash
ONLY=1UBQ bash "$REPO_ROOT/scripts/download_v13_independent_target_pdbs.sh"
```

Default candidates:

```text
1D3Z, then 5DK8
```

The script will not use 1UBI, 3ONS partial, native 1UBQ, or any generated/synthetic fallback.

## Smoke run

```bash
CLEAN=1 REPLICAS=1 STEPS=1000 TARGET_OPEN_STEPS=100 REPORTER_INTERVAL_STEPS=100 \
  bash "$REPO_ROOT/scripts/run_v13a_1ubq_repair_fixed.sh"

python3 "$REPO_ROOT/scripts/check_v13_runtime_provenance_v0.py" \
  "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_REPAIR_FIXED"
```

You need:

```text
certificate_present = true
trajectory_count > 0
successful_replicas > 0
target_sequence_exact_match = true
target_prepared_ca_count = 76
target_prepared_coverage_ratio = 1.0
```

## Official V13a 1UBQ run

```bash
CLEAN=1 REPLICAS=10 MAX_PARALLEL=2 \
  bash "$REPO_ROOT/scripts/run_v13a_1ubq_repair_fixed.sh"

python3 "$REPO_ROOT/scripts/check_v13_runtime_provenance_v0.py" \
  "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_REPAIR_FIXED"
```

## Print result summary

```bash
python3 - <<'PY'
import json, os
root=os.environ["REPO_ROOT"]
p=f"{root}/first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_REPAIR_FIXED/openmm_tmd_replicas_v0_certificate.json"
o=json.load(open(p))
pre=o.get("preflight",{}).get("v5_requirements",{})
print("successful_replicas =", o.get("successful_replicas"))
print("preflight_v5_ready =", o.get("preflight",{}).get("v5_ready"))
print("failed_checks =", pre.get("failed_checks"))
print("effective_strict_count =", pre.get("effective_strict_count"))
print("effective_balanced_count =", pre.get("effective_balanced_count"))
print("effective_rescue_count =", pre.get("effective_rescue_count"))
print("audit_pair_count =", pre.get("audit_pair_count"))
print("runtime_v10_selected_pair_count =", o.get("runtime_v10_selected_pair_count"))
print("selected_pair_count =", o.get("selected_pair_count"))
print("failure_types =", o.get("success_evaluation",{}).get("failure_types"))
print("claim_allowed =", o.get("claim_allowed"))
PY
```
