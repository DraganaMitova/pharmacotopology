# V40 Mechanism Perturbation Pressure Tests Protocol

## Purpose

V40 tests whether non-coordinate mechanism grammar can identify perturbation pressure: which mechanism operators should break, weaken, shift, or block a target mechanism.

This is not atomistic structure prediction, does not run MD, and does not mark the folding problem solved.

## Targets

- `KcsA`
- `XCL1_lymphotactin`
- `alpha_synuclein_SNCA`

## Allowed Inputs

- V36 evidence dossiers
- V37 mechanism maps
- V39 prediction packets and validation results
- UniProt feature, variant, and function annotations
- DisProt annotations
- InterPro/Pfam family signatures
- Literature-derived mutation, state, and function annotations
- Non-coordinate experimental perturbation evidence

## Blocked Inputs

- PDB coordinate contacts
- Coordinate-derived contact CSVs
- Native coordinate metrics
- AlphaFold, ESMFold, or RoseTTAFold predicted coordinates
- V33/V34 KcsA coordinate-derived CSVs
- `first_contact_clean_pharmacotopology_layer_run` files as evidence
- V38 answer key as assignment evidence

## Artifacts

V40 writes target-level inputs under:

- `data/mechanism_perturbations/V40/sources/<target>/perturbation_source_manifest.json`
- `data/mechanism_perturbations/V40/predictions/<target>/perturbation_prediction_packet.json`
- `data/mechanism_perturbations/V40/predictions/<target>/perturbation_table.csv`
- `data/mechanism_perturbations/V40/validation/<target>/perturbation_validation_result.json`
- `data/mechanism_perturbations/V40/validation/<target>/scientist_question_answer.md`

The run-level certificate and report are written under:

- `first_contact_clean_pharmacotopology_layer_run/V40_MECHANISM_PERTURBATION_PRESSURE_TESTS/`

## Validation Levels

- `perturbation_supported_by_holdout`
- `perturbation_partially_supported_clean_abstain`
- `perturbation_contradicted_or_failed`
- `blocked_for_leakage`

## Controls

V40 requires all controls to pass:

1. Prediction perturbations use no coordinate-derived evidence.
2. Validation holdouts use no coordinate-derived evidence.
3. Internal runtime sources are blocked.
4. V38 answer key leakage attempt is blocked.
5. KcsA generic-channel-only perturbation does not preserve hard KcsA grammar.
6. KcsA filter/ion-selectivity removed weakens or invalidates hard KcsA grammar.
7. XCL1 state A removed becomes partial/invalid.
8. XCL1 state B removed becomes partial/invalid.
9. XCL1 mixed-state pooling forced is blocked.
10. SNCA disorder evidence removed becomes partial/invalid.
11. SNCA compact single-fold forcing is blocked.
12. SNCA aggregation context overpromoted to native fold is blocked.
13. Swapped perturbation holdouts are detected or fail validation.
14. Placeholder citations/source rows are blocked.

## Status Rule

V40 passes only if all three targets generate perturbation predictions, each target has at least three perturbation buckets, at least two targets are supported by holdouts, no target is contradicted, all required controls pass, no coordinate/internal/answer-key leakage occurs, and claim flags remain disabled.

The pass status is:

`V40_MECHANISM_PERTURBATION_PRESSURE_PASSED_CLAIM_DISABLED`

The pass does not permit a solved-folding claim.
