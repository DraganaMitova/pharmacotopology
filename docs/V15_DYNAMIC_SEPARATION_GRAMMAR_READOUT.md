# V15 Dynamic Separation Grammar Readout

V15 retires fixed sequence-distance cutoffs as evidence gates.

Sequence separation is still reported as context, but it does not decide whether a residue pair is valid. A pair is interpreted by the role-aware evidence grammar:

- external DCA/coupling support
- geometry reachability from trajectories
- replica persistence
- protein purpose and topology/domain relation
- adaptive chemical policy
- noise/pollution guards
- claim lock or clean abstain

For a single-domain compact object, selected balanced-core evidence is read as compact-core evidence. For a multi-domain composite object, intradomain domain-core and interdomain hinge evidence are separate roles.

This is postprocess-only and claim-safe. It does not rerun molecular dynamics and keeps `claim_allowed = false`.
