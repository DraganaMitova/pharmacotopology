# V16 Pressure Evidence Gap Lock

This stage freezes the evidence gaps after the V16 zero-MD pressure-role transfer readout.

It is **not** a folding/contact evidence claim. It does not run MD, tune thresholds, use native metrics for selection, or upgrade role classification into folding evidence.

## Locked distinction

- `role_classification_passed_targets`: targets whose pressure role was recognized by the locked grammar.
- `positive_folding_evidence_targets`: must remain empty at this stage.

V16 zero-MD pressure-role transfer checks whether the grammar avoids wrong-object readings:

- p53/MDM2 is not read as autonomous p53 TAD compact fold.
- KcsA is not read as soluble single-domain compact core.
- XCL1 is not forced into one false native fold.

## Gap lock target requirements

### p53_TAD_MDM2

Missing before a real evidence test:

- isolated TAD disorder or clean-abstain context
- partner-bound complex context
- interface/contact evidence
- isolated vs MDM2-bound state labels
- external couplings if available
- leakage guard separating autonomous fold from partner-induced fold

Next evidence test:

`partner_induced_interface_or_helix_readout`

### KcsA

Missing before a real evidence test:

- membrane topology annotation
- pore/selectivity-filter annotation
- oligomer/tetramer assembly context
- chain/interface identity
- transmembrane helix roles
- external couplings if available
- guard against soluble-core misread

Next evidence test:

`membrane_pore_role_evidence_readout`

### XCL1_lymphotactin

Missing before a real evidence test:

- state A / state B labels
- monomer/dimer context
- condition labels if available
- state-specific structural/contact evidence
- state-specific external couplings or constraints if available
- leakage guard preventing mixed-state fake core

Next evidence test:

`state_specific_role_separation_readout`

## Claim boundary

`claim_allowed=false` remains mandatory. V17 may become the first pressure evidence panel only after these gaps are supplied or explicitly clean-abstained.
