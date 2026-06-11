# V66 E63 Fast Membrane Repair Replay 200

Status: `V66_E63_FAST_MEMBRANE_REPAIR_PASSED_REVIEW_REQUIRED`
Targets total: `200`
Supported: `130`
Failed accepted: `70`
Abstain: `28`
Accepted accuracy: `0.5930232558139535`
Raw accuracy: `0.65`
Coverage: `0.86`
Controls: `15/15`

## E62 Baseline On Same 200
- supported: `96`
- failed accepted: `104`
- abstain: `7`
- accepted accuracy: `0.46113989637305697`

## E63 Repair Movement
- `supported_delta`: `34`
- `failed_accepted_delta`: `-34`
- `abstain_delta`: `21`
- `accepted_accuracy_delta`: `0.13188335944089657`
- `coverage_delta`: `-0.10499999999999998`

## False Membrane Repair
- `e62_false_membrane_calls`: `24`
- `e63_false_membrane_calls`: `0`
- `e62_peripheral_as_tm`: `7`
- `e63_peripheral_as_tm`: `0`
- `true_tm_missed_under_e63`: `0`
- `oligomer_regressions_remaining_under_e63`: `0`
- `sentinel_regressions_under_e63`: `0`

## Failure Modes
- `membrane_misread`: `70`

## Boundary
V66 is an adaptive 200-target repair replay. It validates E63's membrane/topology repair direction and leaves broad discovery to V67.
