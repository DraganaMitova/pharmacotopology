# V74 E68 RCSB Nonredundant 200 Discovery

V74 ran E68 unchanged on a fresh four-shard RCSB discovery blade.

## Result

| metric | value |
| --- | ---: |
| targets_total | 200 |
| accepted_count | 107 |
| accepted_supported | 11 |
| clean_abstain_supported | 93 |
| failed_accepted | 96 |
| accepted_accuracy | 0.102803738317757 |
| coverage | 0.535 |
| sentinel_regressions | 0 |
| withheld_context_leakage_detected | false |
| controls_passed | true |

## Failure Taxonomy

| failure_mode | count |
| --- | ---: |
| multidomain_allostery | 33 |
| disulfide_secretory_redox_context | 28 |
| coiled_coil_register | 17 |
| repeat_solenoid_topology | 15 |
| closed_beta_topology | 2 |
| domain_swapping | 1 |
| assembly_required_missed | 0 |
| disorder_misread | 0 |
| knotted_topology | 0 |
| membrane_topology_missed_or_misread | 0 |
| metal_ligand_basin | 0 |
| signal_peptide_vs_true_TM | 0 |

## Shards

| shard | accepted | accepted_supported | failed_accepted | clean_abstain_supported | top_missing_word |
| --- | ---: | ---: | ---: | ---: | --- |
| V74A broad | 30 | 11 | 19 | 20 | multidomain_allostery |
| V74B coiled/repeat/solenoid | 32 | 0 | 32 | 18 | coiled_coil_register |
| V74C disulfide/secretory/extracellular | 25 | 0 | 25 | 25 | disulfide_secretory_redox_context |
| V74D multidomain/domain-swap/allostery | 20 | 0 | 20 | 30 | multidomain_allostery |

## Interpretation

V74 is a useful failure, not a pass. E68 preserved the withheld-context boundary and prior repaired classes, but the accepted predictions over-fired on fresh multidomain/allostery targets.

The next revision is `E69_MULTIDOMAIN_ALLOSTERIC_ARCHITECTURE_GRAMMAR`; the next required batch is `V75_MULTIDOMAIN_ALLOSTERY_REPAIR_PANEL_200`.
