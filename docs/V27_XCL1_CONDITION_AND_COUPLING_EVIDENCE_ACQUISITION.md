# V27 XCL1 Condition and Coupling Evidence Acquisition

V27 is a zero-MD acquisition sprint after V26 selected XCL1/lymphotactin for the next mechanism-evidence step.

It locks state-condition labels that are already supported by the separated state contexts, scans for local state-specific coupling/constraint files, and chooses the next V28 panel. It does not synthesize missing couplings, does not run MD, does not force a single fold, and does not claim fold switching.

Main boundary:

- claim_allowed = false
- positive_folding_evidence_found = false
- new_md_executed = false
- no fixed residue cutoff
- no native metric selection
- no mixed-state pooling
- no single-native-state assumption
