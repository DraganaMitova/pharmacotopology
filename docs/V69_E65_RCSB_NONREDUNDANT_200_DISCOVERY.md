# V69 E65 RCSB Nonredundant 200 Discovery

V69 ran E65 as a fresh four-shard discovery batch. It excluded proteins already used in V61, V63, V64, V65, V66, V67, and V68, with V62 excluded as extra safety because it replays the V61 surface.

## Composition

| shard | count |
| --- | ---: |
| V69A_BROAD_RCSB_NONREDUNDANT | 50 |
| V69B_COFACTOR_LIGAND_METAL_ENRICHED | 50 |
| V69C_ASSEMBLY_COMPLEX_OLIGOMER_ENRICHED | 50 |
| V69D_HARD_TOPOLOGY_ENRICHED | 50 |

## Result

| metric | value |
| --- | ---: |
| targets_total | 200 |
| accepted_count | 173 |
| accepted_supported | 81 |
| clean_abstain_supported | 27 |
| failed_accepted | 92 |
| accepted_accuracy | 0.4682080924855491 |
| coverage | 0.865 |
| controls_passed | true |
| sentinel_regressions | 0 |

## Failure Mode

| failure_mode | count |
| --- | ---: |
| metal_cluster_geometry | 49 |
| ligand_locked_basin | 14 |
| disorder_misread | 7 |
| membrane_topology_missed_or_misread | 7 |
| wrong_regime | 7 |
| coiled_coil_register | 4 |
| soluble_beta_barrel_vs_membrane_barrel | 4 |

V69 selected `metal_cluster_geometry` as the dominant next missing word. The next engine revision is `E66_METAL_CLUSTER_AND_LIGAND_LOCKED_BASIN_GRAMMAR`.

## Boundary

V69 is a discovery batch. It does not claim broad protein folding is solved.
