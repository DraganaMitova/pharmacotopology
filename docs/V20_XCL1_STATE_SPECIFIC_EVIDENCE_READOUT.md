# V20 XCL1 State-Specific Evidence Readout

V20 is a zero-MD pressure-evidence sprint for XCL1/lymphotactin as a metamorphic switch object.

It does not try to choose one canonical native state. It keeps state A and state B in separate evidence buckets and forbids mixed-state fake-core selection.

Allowed outcomes:

- `V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED`
- `V20_XCL1_CLEAN_ABSTAIN_STATE_SPECIFIC_EVIDENCE_INSUFFICIENT_CLAIM_DISABLED`

The readout never permits:

- forcing a single fold,
- pooling state A and state B contacts into one fake core,
- claiming fold-switch mechanism recovery without state-specific support,
- fixed residue cutoff selection,
- native-metric selection,
- claim upgrade.

`claim_allowed` remains `false`.
