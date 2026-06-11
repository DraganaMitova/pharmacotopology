# V15 Three-Protein Dynamic Grammar Panel Locked

This is a claim-disabled milestone.

V15 demonstrates a unified, dynamic, role-aware evidence grammar across three protein object types:

- 4AKE: domain-hinge closure object
- 1UBQ: single-domain compact object
- 1CLL: multi-domain composite object

The panel does not claim universal protein folding. It does not claim full hinge/closure recovery for 4AKE or 1CLL. It does not use a fixed residue-distance cutoff as a decision rule. `claim_allowed` remains false.

## Locked evidence interpretation

- 1UBQ: compact balanced-core evidence selected at 23-48 under adaptive chemical policy.
- 1CLL: C-domain compact-core evidence selected at 97-133 under hierarchical topology grammar.
- 4AKE: D4 domain compact-core evidence selected at 124-135 under balanced-candidate dynamic frequency readout; interdomain closure candidates remain monitor-only unless selected by the grammar.

## Important guardrails

- Sequence separation is context only, not a gate.
- Chemical score is role-aware evidence, not a global hard kill.
- No native precision is used to select evidence.
- Visual/PDB-only material is not promoted into evidence.
- Claim locks remain closed.
