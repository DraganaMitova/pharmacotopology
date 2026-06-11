# V71 E66 RCSB Nonredundant 200 Discovery

Status: `V71_E66_RCSB_NONREDUNDANT_200_DISCOVERY_FAILURES_MINED_REVIEW_REQUIRED`
Targets total: `200`
Accepted count: `181`
Accepted supported: `111`
Clean abstain supported: `19`
Failed accepted: `70`
Accepted accuracy: `0.6132596685082873`
Coverage: `0.905`
Controls: `18/18`
Sentinel regressions: `0`
Top failure mode: `disorder_misread`
Top missing Esperanto word: `disorder_misread`
Recommended next engine revision: `E67_DISORDER_BOUNDARY_AND_FOLD_UPON_BINDING_GRAMMAR`
Next required batch: `E67_AND_V72_REPAIR_PANEL`

## Shard Dashboard

| shard | targets_total | accepted_count | accepted_supported | failed_accepted | clean_abstain_supported | accepted_accuracy | coverage | top_failure_mode | top_missing_esperanto_word |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| `V71A_BROAD_RCSB_NONREDUNDANT` | `50` | `45` | `42` | `3` | `5` | `0.9333333333333333` | `0.9` | `membrane_topology_missed_or_misread` | `membrane_topology_missed_or_misread` |
| `V71B_DISORDER_LOW_COMPLEXITY_FLEXIBLE_REGION_ENRICHED` | `50` | `45` | `13` | `32` | `5` | `0.28888888888888886` | `0.9` | `disorder_misread` | `disorder_misread` |
| `V71C_BETA_BARREL_PROPELLER_REPEAT_SOLENOID_ENRICHED` | `50` | `44` | `12` | `32` | `6` | `0.2727272727272727` | `0.88` | `closed_beta_topology` | `closed_beta_topology` |
| `V71D_COILED_COIL_HELIX_BUNDLE_MULTIDOMAIN_ENRICHED` | `50` | `47` | `44` | `3` | `3` | `0.9361702127659575` | `0.94` | `other` | `other` |
| `TOTAL` | `200` | `181` | `111` | `70` | `19` | `0.6132596685082873` | `0.905` | `disorder_misread` | `disorder_misread` |

## Failed Accepted By Failure Mode

| failure_mode | count |
| --- | ---: |
| `disorder_misread` | `31` |
| `coiled_coil_register` | `0` |
| `soluble_beta_barrel_vs_membrane_barrel` | `0` |
| `repeat_solenoid_topology` | `0` |
| `beta_propeller_closure` | `0` |
| `closed_beta_topology` | `30` |
| `wrong_regime` | `0` |
| `membrane_topology_missed_or_misread` | `5` |
| `multidomain_allostery` | `1` |
| `other` | `3` |

## Failed Accepted By Missing Word

| missing_esperanto_word | count |
| --- | ---: |
| `disorder_misread` | `31` |
| `closed_beta_topology` | `30` |
| `membrane_topology_missed_or_misread` | `5` |
| `other` | `3` |
| `multidomain_allostery` | `1` |

## Boundary
V71 is a fresh discovery/mining batch. E66 is fixed during the run; any engine revision belongs to E67 after this certificate.
