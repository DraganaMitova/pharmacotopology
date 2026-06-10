# V13b 1CLL hierarchical purpose topology readout

This is a postprocess-only readout. It does not rerun OpenMM and does not
change the trajectories.

The purpose is to prevent a multi-domain protein from being interpreted as a
single flat interdomain-only object. A multi-domain protein is treated as a
composition of:

- N-domain compact core evidence
- C-domain compact core evidence
- interdomain hinge/closure evidence
- local support
- medium support
- diagnostic shell

For 1CLL, a pair such as 97-133 is no longer considered a topology failure just
because it is intradomain. Under boundaries `1-75,80-148`, it is classified as
`C_domain_compact_core_candidate`. If it has DCA support and replica/tail
support, it can be reported as C-domain core evidence. This does not prove full
interdomain hinge transfer; it reports a domain-core signal while keeping
`claim_allowed=false`.

No protein-name hardcode is required: the classification is based on declared
domain boundaries and sequence separation.
