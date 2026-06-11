# V30 External Constraint and Coupling Acquisition Sprint

V30 is a practical acquisition layer after V29. It scans for local target-specific external coupling/constraint files, preserves the no-MD gate, and selects the next constraint import/preflight panel.

It does not synthesize missing couplings. If KcsA or XCL1 constraints are missing, V30 reports them as missing and keeps MD disabled.

Main outputs:

- `v30_external_constraint_manifest.json`
- `v30_local_constraint_preflight.json`
- `v30_coupling_and_constraint_availability_scan.json`
- `v30_next_acquisition_decision.json`
- `v30_external_constraint_target_rows.json`
- `v30_external_constraint_target_rows.csv`

Interpretation:

V30 does not test a new fold. It identifies which external constraints/couplings must be imported or built before a deeper mechanism test or any MD.
