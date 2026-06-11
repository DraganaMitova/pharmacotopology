# E61 Failure-Driven Engine Revision

E61 is the first engine revision extracted from the V61 frozen RCSB
nonredundant batch.

V61 stayed frozen on `E60` and reported systematic accepted failures:

- cofactor / ligand / metal context missing
- oligomer / assembly context missing
- membrane context too brittle
- weak-sequence cases correctly preserved as abstention candidates

The E61 engine revision adds explicit Protein Esperanto context marks for:

- `cofactor_context`
- `ligand_context`
- `metal_context`
- `heme_context`
- `nucleotide_context`
- `oligomer_context`
- `assembly_context`
- `partner_copy_context`
- `heteromeric_context`
- `homomeric_context`
- strong membrane context marks such as `transmembrane_context`

It also adds real operator and trajectory support for:

- `cofactor_ligand_assisted_stabilization`
- `oligomerization_controlled_folding`

This revision does not rewrite V61 results. V61 remains the frozen discovery
batch. E61 is the next engine version that should be tested on a fresh frozen
replication batch.
