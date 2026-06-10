# V13a 1UBQ Purpose-Gate Readout Commands

This repo is a full clean repo, not an overlay patch.
It includes strict-honest V13 target provenance plus postprocess-only purpose-gate readout.

Official target policy:
- no synthetic target fallback
- no partial target fallback
- no native fallback
- no hardcoded 0.62
- no native-precision threshold selection
- claim_allowed remains false for readout

## Run tests

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

python3 -m compileall -q src scripts
python3 -m pytest -q tests/test_v13a_purpose_gate_readout.py
```

## Prepare strict independent 1UBQ target

```bash
ONLY=1UBQ bash "$REPO_ROOT/scripts/download_v13_independent_target_pdbs.sh"
```

## Fast smoke runtime, only if trajectories do not already exist

```bash
CLEAN=1 REPLICAS=1 STEPS=1000 TARGET_OPEN_STEPS=100 REPORTER_INTERVAL_STEPS=100 \
  bash "$REPO_ROOT/scripts/run_v13a_1ubq_repair_fixed.sh"
```

## Official V13a runtime, only if you need to regenerate trajectories

```bash
CLEAN=1 REPLICAS=10 MAX_PARALLEL=2 CPU_THREADS=4 \
  bash "$REPO_ROOT/scripts/run_v13a_1ubq_repair_fixed.sh"
```

## Purpose-gate readout, postprocess-only

```bash
bash "$REPO_ROOT/scripts/run_v13a_1ubq_purpose_gate_readout.sh"
python3 "$REPO_ROOT/scripts/print_v13a_purpose_gate_summary.py"
```
