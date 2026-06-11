# V62 E61 Repair and Saturation Rerun

Status: `V62_E61_REPAIR_DIRECTIONAL_IMPROVEMENT_REVIEW_REQUIRED`
Batch mode: `repair`
Engine lineage: `E60 -> E61`
Same targets as V61: `True`
Known failure repair probe: `True`

## E60 vs E61

| Metric | E60 V61 | E61 V62 | Delta |
| --- | ---: | ---: | ---: |
| Accepted | `81` | `96` | `15` |
| Supported | `8` | `90` | `82` |
| Failed accepted | `73` | `6` | `-67` |
| Abstain | `19` | `4` | `-15` |
| Accepted accuracy | `0.09876543209876543` | `0.9375` | `0.8387345679012346` |
| Raw accuracy | `0.08` | `0.9` | `0.8200000000000001` |
| Coverage | `0.81` | `0.96` | `0.1499999999999999` |

## Repair Categories
- `repaired`: `82`
- `stable_supported`: `8`
- `unchanged_failure`: `10`

## Failure Mode Distribution Change
- `cofactor_ligand_missing`: E60 `43`, E61 `0`, delta `-43`
- `disorder_misread`: E60 `4`, E61 `4`, delta `0`
- `membrane_misread`: E60 `11`, E61 `0`, delta `-11`
- `oligomer_state_misread`: E60 `13`, E61 `0`, delta `-13`
- `weak_sequence_signal`: E60 `19`, E61 `4`, delta `-15`
- `wrong_regime`: E60 `2`, E61 `2`, delta `0`

## Remaining Target Scores
- `V62_102L_1` predicted `globular_closure` expected `short_region_host_interface_hijacking` label `contradicted`
- `V62_10OP_3` predicted `insufficient_evidence_clean_abstain` expected `short_region_host_interface_hijacking` label `abstained`
- `V62_10OP_4` predicted `insufficient_evidence_clean_abstain` expected `short_region_host_interface_hijacking` label `abstained`
- `V62_10PX_11` predicted `insufficient_evidence_clean_abstain` expected `intrinsic_disorder_phase_separation` label `abstained`
- `V62_10QM_29` predicted `globular_closure` expected `intrinsic_disorder_phase_separation` label `contradicted`
- `V62_10XZ_1` predicted `globular_closure` expected `intrinsic_disorder_phase_separation` label `contradicted`
- `V62_10ZO_1` predicted `globular_closure` expected `short_region_host_interface_hijacking` label `contradicted`
- `V62_11BE_1` predicted `insufficient_evidence_clean_abstain` expected `intrinsic_disorder_phase_separation` label `abstained`
- `V62_11BQ_27` predicted `globular_closure` expected `intrinsic_disorder_phase_separation` label `contradicted`
- `V62_11DG_26` predicted `globular_closure` expected `intrinsic_disorder_phase_separation` label `contradicted`

## Boundary
V62 is a same-target E61 repair probe. It may show whether E61 repaired V61 failure classes, but it does not license a broad protein-space claim. Coordinates, contacts, ligand geometry, AlphaFold-style models, holdout annotations, and internal runtime artifacts remain blocked before sealing.
