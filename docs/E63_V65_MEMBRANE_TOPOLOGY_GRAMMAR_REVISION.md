# E63 V65 Membrane Topology Grammar Revision

E63 is derived from the V65 membrane topology panel.

V65/E62 showed:

- true transmembrane missed: `0`
- soluble hydrophobic false membrane calls: `14`
- peripheral false transmembrane calls: `7`
- failed accepted: `21`

E63 changes the selector grammar without adding mechanism classes or operators:

- negative membrane-topology evidence and hydrophobicity-only context now block
  generic membrane promotion
- monotopic/peripheral membrane association now abstains unless an explicit
  cofactor or oligomer explanation is present
- cofactor and oligomer context can explain hydrophobic pockets or interfaces
  before topology-conflict abstention
- strong explicit transmembrane/topology context still selects membrane grammar

The E63 certificate is:

- `data/protein_esperanto_engine/E63/e63_v65_membrane_topology_grammar_revision_certificate.json`

E63 is not a broad claim gate. It is a narrow grammar repair that must be
replayed before any V66 1000-target expansion.
