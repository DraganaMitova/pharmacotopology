# GEODESIC_CONSENSUS_FILTER_MD_V0

This batch tests whether the dense-anchor collapse came from naive geodesic filling.
It keeps the structure-first C-alpha relaxation path, but adds native-free filters before
geodesic contacts become restraints:

- contact-map geodesic distance filter;
- sequence-only bending-energy filter;
- deterministic multi-path consensus filter.

Native/coordinate truth remains audit-only after final structure-to-contact conversion.
No ESMFold/AlphaFold/template/learned geometry prior is used before selection.

Claim gate remains strict: precision, recall, and long-range recall must all be >= 0.70.
MD implementation remains dependency-free and bounded; no full suite is required for this layer.

## Observed default 8-row result

- mean precision: 0.089742
- mean recall: 0.605346
- mean long-range recall: 0.479945
- mean F1: 0.155941
- solved: false
- rejection: precision below 0.70

## 4AKE default result

- precision: 0.092581
- recall: 0.513417
- long-range recall: 0.181347
- F1: 0.156873
- solved: false

## Strict precision probe

With only a tiny accepted geodesic layer, precision rises relative to the dense default but long-range order collapses:

- mean precision: 0.126171
- mean recall: 0.339065
- mean long-range recall: 0.182184
- 4AKE precision: 0.187175
- 4AKE recall: 0.386404
- 4AKE long-range recall: 0.056995

## Interpretation

Consensus filtering rejects many naive geodesic fills, but it does not solve the central tradeoff. The system can be tuned toward precision or toward recall, but not both at the ESMFold-like level. This is a negative result for simple geodesic filtering, bending-cost filtering, and small deterministic path consensus as replacements for a learned dense global geometry prior.
