# V57 Blind Protein Esperanto Generalization Gate

V57 freezes the V50-V56 engine and tests fresh targets that were not used to extract the grammar.

## Frozen Engine Rule

Allowed changes:

- new target manifests
- new source manifests
- new holdout manifests
- new reports
- new certificates
- new tests

Not allowed:

- new operators
- new mechanism classes
- new scoring rules
- new control definitions
- README updates
- post-result tuning

If the frozen engine needs tuning to pass V57, the correct status is `V57_GENERALIZATION_BLOCKED_ENGINE_NEEDS_REVISION`.

## Fresh Target Battery

V57 uses five fresh regime slots:

| Slot | Target |
| --- | --- |
| globular closure | Streptococcal protein G B1 domain |
| IDP / phase or disorder | hnRNPA1-like prion LCD |
| membrane / proteostasis | human rhodopsin P23H gate |
| metamorphic / switching | XCL1 lymphotactin fold-switch gate |
| short interface / hijacking | HIV-1 Tat host-interface gate |

## Pre-Registered Outputs

Before holdout validation, each target must produce:

- selected grammar
- operator field
- trajectory type
- dominant regions
- predicted perturbation directions
- falsification criteria
- abstention criteria

## Required Controls

Every target must reject one forced wrong grammar. Every target must separate correct perturbations from wrong-region or neutral controls. Global controls must preserve coordinate leakage blocking, holdout-before-seal blocking, random-sequence abstention, and `folding_problem_solved = false`.

## Status Values

- `V57_BLIND_GENERALIZATION_PASSED_REVIEW_REQUIRED`
- `V57_PARTIAL_GENERALIZATION_WITH_ABSTENTIONS_REVIEW_REQUIRED`
- `V57_GENERALIZATION_BLOCKED_ENGINE_NEEDS_REVISION`
- `V57_BLOCKED_FOR_LEAKAGE`

## Claim Boundary

Allowed if passed:

```text
The frozen Protein Esperanto engine generalized to fresh targets across multiple regimes under leakage-controlled blind validation.
```

Still forbidden:

```text
Universal protein folding is solved.
Coordinates were predicted de novo.
Atomistic folding was solved.
External review is unnecessary.
```
