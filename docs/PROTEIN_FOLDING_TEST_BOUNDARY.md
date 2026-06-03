# Protein Folding Test Boundary

This branch adds a small protein-folding hypothesis lab to the project. The
goal is not to solve folding. The goal is to create a safe benchmark shell
where topology signatures can be compared against folding-relevant reference
summaries.

The test claim is deliberately narrow:

```text
topology signatures extracted from protein sequences and structures
can be compared against folding-relevant benchmark summaries
under explicit uncertainty and failure labels
```

## Modeled Chain

```text
protein sequence
-> predicted topology signature
-> reference topology signature
-> contact-map proxy similarity
-> fold-class match
-> uncertainty radius
-> evidence readiness
```

The signature dimensions are:

```text
sequence_complexity
secondary_structure_balance
contact_map_closure
hydrophobic_core_closure
loop_disorder_pressure
domain_boundary_stability
long_range_contact_order
conformational_flexibility
knot_or_entanglement_signature
uncertainty_radius
```

## Current Benchmark Status

The default benchmark rows are placeholders. They are useful for exercising the
code path, output schema, tests, and safety boundary, but they are not external
validation.

Real calibration would require reference rows derived from external structure
resources such as contact maps, fold-class annotations, structure ensembles, or
curated folding benchmark sets. Those sources are intentionally not bundled in
this dependency-free prototype.

## Output Fields

Each benchmark comparison writes:

```text
protein_id
sequence_length
reference_structure_source
predicted_topology_signature
reference_topology_signature
contact_map_similarity
fold_class_match
uncertainty_radius
evidence_readiness
failure_reason
```

The default `failure_reason` is expected to say that an external structure
benchmark is not attached. That is a useful failure: it prevents placeholder
numbers from being mistaken for evidence.

## Boundary

This module does not:

```text
solve protein folding
predict atomic structure
generate protein sequences
design molecules
recommend compounds
make clinical claims
```

It is a benchmark adapter and hypothesis-review shell. Its best use is to make
the next evidence step obvious.
