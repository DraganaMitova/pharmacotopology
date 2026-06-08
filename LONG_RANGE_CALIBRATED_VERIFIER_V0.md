# LONG_RANGE_CALIBRATED_VERIFIER_V0

This batch tests the hypothesis that the previous calibrated verifier erased
long-range topology because it was trained mostly to reject decoy long-range
contacts. The new layer builds a leave-one-target-out long-range contact
potential from reference coordinate rows, then uses that potential during
structure-level contact selection.

## What changed

- Added `src/pharmacotopology/folding_long_range_calibrated_verifier.py`.
- Added `scripts/run_long_range_calibrated_verifier_v0.py`.
- Added `scripts/compute_long_range_potential_v0.py`.
- Added `tests/test_long_range_calibrated_verifier_v0.py`.
- Added locked artifacts under
  `first_contact_clean_pharmacotopology_layer_run/long_range_calibrated_verifier_v0/`.

## Boundary

- No ESMFold/AlphaFold/template geometry is used before selection.
- Target native contacts/coordinates are excluded from target calibration.
- Native truth is attached only after contact selection for audit.
- The long-range potential is PDB-derived from other locked rows, therefore the
  result cannot claim a universal physical law.

## Result summary

For the focused 4AKE run, the long-range calibrated verifier did not solve the
folding problem. It avoided the previous total long-range wipeout, but the
rescued long-range signal remained too weak and too noisy.

Claim gate remains closed.
