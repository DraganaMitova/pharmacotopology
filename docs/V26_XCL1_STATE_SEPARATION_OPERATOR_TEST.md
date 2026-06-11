# V26 XCL1 State-Separation Operator Test

V26 is the first locked mechanism-operator test selected by V25. It does not run MD and does not claim fold switching or universal folding. It tests the repeated `state_separation` operator on XCL1/lymphotactin.

Pass requires:

- state A detected;
- state B detected;
- state-specific buckets preserved;
- no mixed-state contact pooling;
- no mixed-state fake core;
- no forced single native state;
- no fold-switch claim;
- `claim_allowed = false`.

The next step after a pass is condition/coupling acquisition, not MD.
MD remains forbidden until state conditions and state-specific external constraints are locked and false-win risk is low.
