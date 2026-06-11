# V64 E62 RCSB 500 Membrane Repair Replay

Status: `V64_E62_MEMBRANE_REPAIR_DIRECTIONAL_IMPROVEMENT_REVIEW_REQUIRED`
Targets total: `500`
Accepted: `500`
Supported: `328`
Failed accepted: `172`
Abstain: `0`
Accepted accuracy: `0.656`
Raw accuracy: `0.656`
Coverage: `1.0`
Controls: `20/20`
Engine modified during batch: `False`

## E61 Baseline
- supported: `238`
- failed accepted: `262`
- membrane_misread: `216`

## E62 Replay
- supported: `328`
- failed accepted: `172`
- membrane_misread: `123`
- new failure modes from regressions: `{'oligomer_state_misread': 3}`
- abstain: `0`
- accepted accuracy: `0.656`

## Repair Map
- `stable_supported`: `235`
- `repaired`: `93`
- `persistent_failure`: `169`
- `new_failure`: `3`

## Failure-Mode Movement
- `disorder_misread`: `15 -> 15` (`0`)
- `membrane_misread`: `216 -> 123` (`-93`)
- `oligomer_state_misread`: `2 -> 5` (`3`)
- `wrong_regime`: `29 -> 29` (`0`)

## Soluble False Membrane Check
- soluble false membrane calls: `3`
- cofactor/oligomer regressions: `3`

## Mechanism Distribution
- `cofactor_ligand_assisted_stabilization`: predicted `175`, expected `112`
- `fold_upon_binding_disorder`: predicted `1`, expected `0`
- `globular_closure`: predicted `2`, expected `1`
- `intrinsic_disorder_phase_separation`: predicted `0`, expected `15`
- `membrane_multidomain_folding_proteostasis`: predicted `99`, expected `219`
- `oligomerization_controlled_folding`: predicted `222`, expected `124`
- `short_region_host_interface_hijacking`: predicted `1`, expected `29`

## Boundary
V64 is a paired E62 repair replay on the exact V63 target set. It measures repair direction and regression shape; it does not make a broad saturation or solved claim.
