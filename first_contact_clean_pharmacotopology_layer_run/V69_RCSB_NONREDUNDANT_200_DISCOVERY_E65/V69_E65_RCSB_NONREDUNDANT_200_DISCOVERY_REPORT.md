# V69 E65 RCSB Nonredundant 200 Discovery

Status: `V69_E65_RCSB_NONREDUNDANT_200_DISCOVERY_FAILURES_MINED_REVIEW_REQUIRED`
Targets total: `200`
Accepted count: `173`
Accepted supported: `81`
Clean abstain supported: `27`
Failed accepted: `92`
Accepted accuracy: `0.4682080924855491`
Coverage: `0.865`
Controls: `18/18`
Sentinel regressions: `0`
Top missing Esperanto word: `metal_cluster_geometry`
Recommended next engine revision: `E66_METAL_CLUSTER_AND_LIGAND_LOCKED_BASIN_GRAMMAR`
Next required batch: `E66_AND_V70_REPAIR_PANEL`

## Shard Dashboard

| shard | targets_total | accepted_count | accepted_supported | failed_accepted | clean_abstain | accepted_accuracy | coverage | top_missing_esperanto_word |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `V69A_BROAD_RCSB_NONREDUNDANT` | `50` | `40` | `33` | `7` | `10` | `0.825` | `0.8` | `ligand_locked_basin` |
| `V69B_COFACTOR_LIGAND_METAL_ENRICHED` | `50` | `50` | `8` | `42` | `0` | `0.16` | `1.0` | `metal_cluster_geometry` |
| `V69C_ASSEMBLY_COMPLEX_OLIGOMER_ENRICHED` | `50` | `39` | `13` | `26` | `11` | `0.3333333333333333` | `0.78` | `ligand_locked_basin` |
| `V69D_HARD_TOPOLOGY_ENRICHED` | `50` | `44` | `27` | `17` | `6` | `0.6136363636363636` | `0.88` | `ligand_locked_basin` |
| `TOTAL` | `200` | `173` | `81` | `92` | `27` | `0.4682080924855491` | `0.865` | `metal_cluster_geometry` |

## Failed Accepted By Failure Mode

| failure_mode | count |
| --- | ---: |
| `metal_cluster_geometry` | `49` |
| `ligand_locked_basin` | `14` |
| `disorder_misread` | `7` |
| `membrane_topology_missed_or_misread` | `7` |
| `wrong_regime` | `7` |
| `coiled_coil_register` | `4` |
| `soluble_beta_barrel_vs_membrane_barrel` | `4` |

## Boundary
V69 is a fresh discovery/mining batch. E65 is fixed during the run; any engine revision belongs to E66 after this certificate.
