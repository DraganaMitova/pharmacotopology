# E62 V63 Failure-Driven Engine Revision

E62 is derived from the V63 RCSB 500 discovery batch.

V63 stayed on `E61` and reported:

- `targets_total = 500`
- `supported_count = 238`
- `failed_accepted_count = 262`
- top failure mode: `membrane_misread = 216`

The V63 failure distribution showed that E61 often selected:

- `oligomerization_controlled_folding` for membrane targets with incidental
  copy/assembly metadata
- `cofactor_ligand_assisted_stabilization` for membrane targets with incidental
  ion/ligand metadata

E62 changes the engine grammar priority:

```text
strong membrane/topology context outranks incidental ligand/cofactor or assembly context
```

This is a grammar-priority correction, not a claim gate. V63 remains the E61
discovery result. E62 must be tested in the next repair/expansion batch before
any broader claim is allowed.
