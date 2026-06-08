# CALIBRATED_TOPOLOGY_VERIFIER_V0

This batch tests the next hypothesis: a generic global verifier is not enough; the verifier must be calibrated toward real fold topology.

## Implemented

- Optional PyRosetta/DFIRE backend probe.
- Dependency-light fallback when PyRosetta/DFIRE are unavailable.
- Leave-one-target-out transparent structural-feature discriminator.
- Features: geometry score, DCA/restraint confidence, residue-pair potential, local patch support, sequence separation/contact order, C-alpha distance closeness, degree pressure, long-range flag.
- Calibrated global selector: statistical potential surrogate + degree coherence + loop/patch coherence + restraint satisfaction + compactness/collision + contact order criterion.
- Native/contact truth is excluded for the target row before selection. Other locked rows may be used only for calibration, so this cannot claim a universal physical law.

## Backend status in this sandbox

PyRosetta: unavailable.
DFIRE: unavailable.
Backend used: transparent leave-one-target-out structural-feature discriminator.

## Focused 4AKE result

Command:

```bash
PYTHONPATH=src timeout 75 python scripts/run_calibrated_topology_verifier_v0.py \
  --source-accession 4AKE:A \
  --md-steps 24 \
  --candidate-count 2 \
  --max-direct 48 \
  --max-sequence-closure 48 \
  --max-geodesic 48
```

Result:

```text
mean_native_contact_precision_after_audit = 0.331818
mean_native_contact_recall_after_audit = 0.261181
mean_long_range_contact_recall_after_audit = 0.0
mean_contact_map_f1_after_audit = 0.292292
folding_problem_solved = false
universal_physical_law_claim_allowed = false
claim_rejection_reason = calibrated_verifier_claim_rejected_for_rows:4AKE:A
```

Interpretation:

The calibrated verifier improves 4AKE precision compared with the generic verifier, but it deletes/collapses the long-range topology. It is therefore not a replacement for the learned global geometry prior.

## Validation

```text
focused tests: 3 passed
compileall: OK
full_suite_run = false
gif_generation_used = false
predictor_subprocesses_used = false
atomistic_md_engine_used = false
final CPU cleanup matched_count = 0
```
