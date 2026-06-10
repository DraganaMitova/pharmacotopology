# V14 Unified Protein Esperanto Grammar Readout

V14 is a postprocess-only readout. It does not rerun OpenMM, does not tune thresholds, and does not enable biological transfer claims.

It reads already-generated artifacts for:

- 4AKE: domain/hinge object, role-aware steering reproduction artifacts when present.
- 1UBQ: single-domain compact object, adaptive chemical-policy readout.
- 1CLL: multi-domain composite object, hierarchical domain-core/interdomain topology readout.

The shared grammar axes are:

1. External DCA/coupling signal.
2. Geometry reachability from trajectories.
3. Replica persistence.
4. Purpose/topology role assignment.
5. Role-aware chemical policy.
6. Noise and long-range pollution guards.
7. Claim lock or clean abstain.

V14 should not print `protein folding solved`. The intended milestone is:

```text
unified role-aware evidence grammar is coherent across multiple protein object types;
claim_allowed = false
```

If a source artifact is missing, V14 reports that artifact as missing instead of inventing evidence.
MD trajectories are reused. New simulations are unnecessary for V14.
