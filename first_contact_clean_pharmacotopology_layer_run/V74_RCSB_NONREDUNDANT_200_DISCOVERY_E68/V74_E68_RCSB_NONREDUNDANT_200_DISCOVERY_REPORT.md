# V74 E68 RCSB Nonredundant 200 Discovery

Status: `V74_E68_RCSB_NONREDUNDANT_200_DISCOVERY_FAILURES_MINED_REVIEW_REQUIRED`
Targets total: `200`
Accepted count: `107`
Accepted supported: `11`
Clean abstain supported: `93`
Failed accepted: `96`
Accepted accuracy: `0.102803738317757`
Coverage: `0.535`
Controls: `18/18`
Sentinel regressions: `0`
Withheld context leakage detected: `False`
Top failure mode: `multidomain_allostery`
Top missing Esperanto word: `multidomain_allostery`
Recommended next engine revision: `E69_MULTIDOMAIN_ALLOSTERIC_ARCHITECTURE_GRAMMAR`
Next required batch: `E69_AND_V75_REPAIR_PANEL`

## Shard Dashboard

| shard | targets_total | accepted_count | accepted_supported | failed_accepted | clean_abstain_supported | accepted_accuracy | coverage | top_failure_mode | top_missing_esperanto_word |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| `V74A_BROAD_RCSB_NONREDUNDANT` | `50` | `30` | `11` | `19` | `20` | `0.36666666666666664` | `0.6` | `multidomain_allostery` | `multidomain_allostery` |
| `V74B_COILED_COIL_REPEAT_SOLENOID_ENRICHED` | `50` | `32` | `0` | `32` | `18` | `0.0` | `0.64` | `coiled_coil_register` | `coiled_coil_register` |
| `V74C_DISULFIDE_SECRETORY_EXTRACELLULAR_ENRICHED` | `50` | `25` | `0` | `25` | `25` | `0.0` | `0.5` | `disulfide_secretory_redox_context` | `disulfide_secretory_redox_context` |
| `V74D_MULTIDOMAIN_DOMAIN_SWAP_ALLOSTERY_UNUSUAL_ENRICHED` | `50` | `20` | `0` | `20` | `30` | `0.0` | `0.4` | `multidomain_allostery` | `multidomain_allostery` |
| `TOTAL` | `200` | `107` | `11` | `96` | `93` | `0.102803738317757` | `0.535` | `multidomain_allostery` | `multidomain_allostery` |

## Failed Accepted By Failure Mode

| failure_mode | count |
| --- | ---: |
| `coiled_coil_register` | `17` |
| `repeat_solenoid_topology` | `15` |
| `disulfide_secretory_redox_context` | `28` |
| `domain_swapping` | `1` |
| `multidomain_allostery` | `33` |
| `signal_peptide_vs_true_TM` | `0` |
| `knotted_topology` | `0` |
| `membrane_topology_missed_or_misread` | `0` |
| `assembly_required_missed` | `0` |
| `closed_beta_topology` | `2` |
| `metal_ligand_basin` | `0` |
| `disorder_misread` | `0` |
| `other` | `0` |

## Failed Accepted By Missing Word

| missing_esperanto_word | count |
| --- | ---: |
| `multidomain_allostery` | `33` |
| `disulfide_secretory_redox_context` | `28` |
| `coiled_coil_register` | `17` |
| `repeat_solenoid_topology` | `15` |
| `closed_beta_topology` | `2` |
| `domain_swapping` | `1` |

## Boundary
V74 is a fresh discovery/mining batch. E68 is fixed during the run; withheld context isolation is a required control.
