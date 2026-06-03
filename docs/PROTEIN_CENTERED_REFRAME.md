# Protein-Centered Reframe

Pharmacotopology is being reframed from compound-first pharmacology toward
protein-centered mechanism topology.

The compound is not the main modeled object. A compound is treated only as an
abstract perturbation source. The modeled layer is:

```text
abstract_compound_class
-> protein_family
-> protein_state_shift
-> pathway/network perturbation
-> topology_delta
-> collapse_cost
-> uncertainty_radius
-> evidence_readiness
```

This keeps the workbench closer to mechanism literacy. Receptors, enzymes, ion
channels, transporters, antibodies, cytokines, and signaling proteins are the
places where abstract perturbations become biological state changes. In this
repository, those state changes remain simulated and hypothesis-only.

## Boundary

This is not a drug-design tool. It does not propose molecules, optimize
compounds, recommend treatments, infer patient state, map brand names, or make
clinical claims.

The useful question is:

```text
What kind of protein mechanism is being touched,
what abstract state shift is being assumed,
what topology dimensions move,
and what evidence would be needed before any claim could exist?
```

The unsafe question is:

```text
Which real substance should someone use?
```

That question is outside the project boundary.

## Current Metadata

Each mechanism vector now carries:

```text
abstract_compound_class
protein_family
protein_mechanism_class
protein_state_shift
pathway_network_perturbation
```

These fields are labels for hypothesis review. They are not validated target
claims.
