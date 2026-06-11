# V63 RCSB 500 Discovery Batch

Status: `V63_RCSB_500_DISCOVERY_FAILURES_REVIEW_REQUIRED`
Targets total: `500`
Accepted: `500`
Supported: `238`
Failed accepted: `262`
Abstain: `0`
Accepted accuracy: `0.476`
Raw accuracy: `0.476`
Coverage: `1.0`
Controls: `19/19`
Engine modified during batch: `False`
README check skipped by user instruction: `True`

## Top Failure Modes
- `membrane_misread`: `216`
- `wrong_regime`: `29`
- `disorder_misread`: `15`
- `oligomer_state_misread`: `2`

## Mechanism Distribution
- `cofactor_ligand_assisted_stabilization`: predicted `200`, expected `112`
- `fold_upon_binding_disorder`: predicted `1`, expected `0`
- `globular_closure`: predicted `2`, expected `1`
- `intrinsic_disorder_phase_separation`: predicted `0`, expected `15`
- `membrane_multidomain_folding_proteostasis`: predicted `3`, expected `219`
- `oligomerization_controlled_folding`: predicted `293`, expected `124`
- `short_region_host_interface_hijacking`: predicted `1`, expected `29`

## Boundary
V63 is a broad discovery/mining batch. It records E61 behavior on 500 nonredundant RCSB protein entities and preserves failures for E62 grammar mining. It does not make a broad saturation or solved claim.
