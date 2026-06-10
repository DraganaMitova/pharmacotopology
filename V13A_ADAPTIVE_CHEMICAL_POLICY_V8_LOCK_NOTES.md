# V13a Adaptive Chemical Policy V8 Lock Notes

This clean repo keeps the V13a 1UBQ adaptive chemical policy readout claim-safe.

Changes relative to V7:

- Replaced misleading `dca_enrichment_pass` selector check with explicit `dca_absolute_support_pass`.
- Kept `dca_background_enrichment_ratio` and `dca_background_enrichment_pass` as claim-lock diagnostics only.
- Added `dca_pass_semantics = absolute_support_for_selection_background_enrichment_for_claim_lock_only`.
- Added `claim_lock_check` to the certificate.
- The adaptive readout may find a purpose-fit core while the claim lock remains blocked if background enrichment is not passed.
- `claim_allowed`, `physics_interpretation_allowed`, and `biological_transfer_claim_allowed` remain false.

No new MD is required if `V13a_1UBQ_REPAIR_FIXED` trajectories already exist.
