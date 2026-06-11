# V33 Negative Controls Patch Note

This patch adds the claim-disabled V33 negative-control battery and repairs one V32 KcsA bucket-gating weakness caught by the controls.

## What changed

- Added `scripts/run_v33_negative_controls_v0.py`.
- Added `scripts/print_v33_negative_controls.py`.
- Added `tests/test_v33_negative_controls.py`.
- Added a V32 regression test: interface-only KcsA evidence must not satisfy the pore/filter bucket.
- Tightened the V32 KcsA readiness check so a generic `contact` word in an assembly/interface row cannot count as pore/filter evidence.

## Scientific boundary

This is still not a folding claim. Passing V33 negative controls means the V32/V33 evidence grammar blocks:

- internal runtime source poisoning,
- annotation-only promotion,
- missing pore/filter bucket,
- missing assembly/interface bucket,
- empty operator buckets,
- wrong-target V33 scope misuse.

`claim_allowed`, `new_MD_allowed`, and `folding_problem_solved` must remain false.
