# Real Coordinate Visual Contact Boundary

This layer upgrades the visual proof target from toy locked contact pairs to
native contact maps derived from real coordinate structures.

It does not solve protein folding.

## What Changed

The `real_coordinate_visual_8` benchmark uses eight coordinate-backed rows from
the locked real-10 reference set:

```text
1MBN:A
2LZM:A
1TEN:A
1CSP:A
1TIM:A
1PGA:A
4AKE:A
1CLL:A
```

For each row, the benchmark stores a locked C-alpha coordinate trace derived
from RCSB PDB coordinate files. Native contact maps are then extracted from
those coordinates at runtime:

```text
C-alpha distance <= 8.0 angstrom
minimum sequence separation >= 3
```

The benchmark does not store precomputed native contact pairs in the lock file.
The pairs are derived from the locked coordinates and checked against a native
contact-map hash.

## Prediction Boundary

Prediction remains sequence-only:

```text
sequence
-> sequence-only contact candidates
-> only then coordinate-derived native scoring
```

Coordinates, native contacts, fold labels, and truth axes are not predictor
inputs.

## Claim Boundary

Allowed claim:

```text
The project can score sequence-only contact hypotheses against coarse native
C-alpha contact maps extracted from real coordinate structures.
```

Forbidden claim:

```text
The project has discovered the folding mechanism.
The project can predict atomistic protein folds.
The project has solved protein folding.
```

Required flags:

```text
real_coordinate_native_contacts_extracted = true
toy_locked_contact_targets_used = false
coarse_ca_only = true
full_atomic_folding_available = false
repair_heuristic_applied = false
mechanism_discovery_claim_allowed = false
global_folding_claim_allowed = false
folding_problem_solved = false
```

## Why This Matters

The previous visual workbench was useful because failures became visible, but
the native contacts were a toy/coarse internal target. This benchmark makes the
target stricter: native contacts now come from coordinate geometry.

The honest outcome may look worse numerically. That is acceptable. The goal is
not to make the model look good; it is to expose where the sequence-only contact
hypothesis fails against a more realistic target.

## Still Not Enough

This layer is still coarse:

```text
C-alpha only
single-chain coordinate traces
no all-atom side-chain contacts
no solvent, dynamics, ligand, membrane, or ensemble modeling
no learned structural model
no blind external repair validation
```

The next biology work should be driven by visible failure cohorts on this
coordinate benchmark, not by the earlier toy benchmark.
