# V13 Runtime Test Now

This clean repo version fixes the V13 target-preparation problem that caused
`successful_replicas=0/10` and `trajectory=null` before the mechanism could be tested.

## Why V2 failed

AlphaFold DB returned 404 for the mature ubiquitin IDs used by the downloader.
That made the prepared target PDB missing, so the V13 wrappers refused to run.

## What V3 does

- V13a 1UBQ now defaults to a very similar independent RCSB ubiquitin structure:
  - first `3ONS`, then `1UBI`
- If RCSB is unavailable and you only need to test runtime plumbing, you may use:
  - `ALLOW_NATIVE_TARGET_DEBUG=1`
  - This can fall back to bundled `data/rcsb_pdb/1UBQ.pdb`.
  - That mode is **not** an independent biological transfer result.
- V13b 1CLL tries AFDB first, then RCSB `1CLL`.

## Run first: 1UBQ only

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts"

python3 -m compileall -q src scripts

ONLY=1UBQ bash "$REPO_ROOT/scripts/download_v13_independent_target_pdbs.sh"
```

If RCSB also fails, use debug runtime-only fallback:

```bash
ONLY=1UBQ ALLOW_NATIVE_TARGET_DEBUG=1 bash "$REPO_ROOT/scripts/download_v13_independent_target_pdbs.sh"
```

Then smoke run:

```bash
CLEAN=1 REPLICAS=1 STEPS=1000 TARGET_OPEN_STEPS=100 REPORTER_INTERVAL_STEPS=100 \
  bash "$REPO_ROOT/scripts/run_v13a_1ubq_repair_fixed.sh"

python3 "$REPO_ROOT/scripts/check_v13_runtime_provenance_v0.py" \
  "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_REPAIR_FIXED"
```

Then official run:

```bash
CLEAN=1 REPLICAS=10 MAX_PARALLEL=2 \
  bash "$REPO_ROOT/scripts/run_v13a_1ubq_repair_fixed.sh"

python3 "$REPO_ROOT/scripts/check_v13_runtime_provenance_v0.py" \
  "$REPO_ROOT/first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_REPAIR_FIXED"
```

## Honest interpretation

The result is not valid unless there are trajectories:

```text
certificate_present = true
trajectory_count > 0
successful_replicas > 0
```

If those are true, then inspect V5/V10:

```bash
python3 - <<'PY'
import json, os
root=os.environ["REPO_ROOT"]
p=f"{root}/first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_REPAIR_FIXED/openmm_tmd_replicas_v0_certificate.json"
o=json.load(open(p))
print("successful_replicas =", o.get("successful_replicas"))
print("preflight_v5_ready =", o.get("preflight",{}).get("v5_ready"))
print("failed_checks =", o.get("preflight",{}).get("v5_requirements",{}).get("failed_checks"))
print("effective_strict_count =", o.get("preflight",{}).get("v5_requirements",{}).get("effective_strict_count"))
print("effective_balanced_count =", o.get("preflight",{}).get("v5_requirements",{}).get("effective_balanced_count"))
print("effective_rescue_count =", o.get("preflight",{}).get("v5_requirements",{}).get("effective_rescue_count"))
print("audit_pair_count =", o.get("preflight",{}).get("v5_requirements",{}).get("audit_pair_count"))
print("runtime_v10_selected_pair_count =", o.get("runtime_v10_selected_pair_count"))
print("failure_types =", o.get("success_evaluation",{}).get("failure_types"))
print("claim_allowed =", o.get("claim_allowed"))
PY
```
