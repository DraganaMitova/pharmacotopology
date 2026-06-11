# V70 E66 Metal Ligand Repair Panel 200

V70 is the repair panel for V69's dominant missing word: `metal_cluster_geometry`, with `ligand_locked_basin` repaired in the same E66 grammar.

## Composition

| group | count |
| --- | ---: |
| V69_METAL_CLUSTER_FAILURE_REPLAY | 49 |
| V69_LIGAND_LOCKED_FAILURE_REPLAY | 14 |
| METAL_CLUSTER_POSITIVE_EXPANSION | 50 |
| LIGAND_LOCKED_POSITIVE_EXPANSION | 30 |
| TRUE_TM_METAL_CONFLICT_SENTINEL | 25 |
| ASSEMBLY_REQUIRED_METAL_CONFLICT_SENTINEL | 20 |
| V69_NON_METAL_FAILURE_TRACKING | 12 |

## Result

| metric | value |
| --- | ---: |
| status | V70_E66_METAL_CLUSTER_LIGAND_REPAIR_PANEL_PASSED_REVIEW_REQUIRED |
| targets_total | 200 |
| accepted_supported | 193 |
| clean_abstain_supported | 0 |
| failed_accepted | 2 |
| targeted_failed_accepted | 0 |
| accepted_accuracy | 0.9897435897435898 |
| coverage | 0.975 |
| controls_passed | true |
| sentinel_regressions | 0 |

## Repair Checks

| check | value |
| --- | ---: |
| V69 metal failures repaired | 49 / 49 |
| V69 ligand failures repaired | 14 / 14 |
| true TM preserved | 25 / 25 |
| assembly required preserved | 20 / 20 |
| non-metal tracking remaining | 7 |

The remaining failures are carried-forward non-metal tracking rows. They are not E66 metal/ligand repair regressions.

## Next

`V71_RCSB_NONREDUNDANT_200_DISCOVERY_E66`
