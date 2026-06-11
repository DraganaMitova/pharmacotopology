# V25 Fast Mechanism Evidence Sprint

V25 changes mode from target-by-target micro-locks to one fast mechanism sprint.
It does not run MD, does not tune thresholds, does not upgrade any folding claim,
and keeps `claim_allowed=false`.

## Outputs

- `v25_kcsa_coupling_interface_readout.json`
- `v25_xcl1_state_specific_readout.json`
- `v25_unified_mechanism_operator_table.json`
- `v25_next_mechanism_test_decision.json`
- `v25_fast_mechanism_evidence_sprint_certificate.json`

## Mechanism operators

The operator table uses a shared operator alphabet:

- `role_detect`
- `context_detect`
- `evidence_assign`
- `reachability_check`
- `replica_support`
- `chemical_confidence`
- `topology_context`
- `state_separation`
- `interface_context`
- `clean_abstain`
- `pollution_guard`

Not every protein uses every operator. The point is that different protein
purposes are read through the same bounded operator alphabet, not a single
fixed residue cutoff or single fixed chemical threshold.

## Claim boundary

V25 is evidence for a repeated role/context/operator grammar. It is not a
universal protein-folding claim. `positive_folding_evidence_targets` remains
empty and `claim_allowed` remains false.
MD is only allowed later if evidence readiness and false-win risk allow it.
