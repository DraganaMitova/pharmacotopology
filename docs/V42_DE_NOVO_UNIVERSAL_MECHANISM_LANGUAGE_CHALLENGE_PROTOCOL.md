# V42 De Novo Universal Mechanism Language Challenge Protocol

## Purpose

V42 is the first sealed de novo mechanism-language benchmark. It forces the system to predict mechanism class, operator regions, perturbation pressure, low-resolution structural or ensemble consequence, and falsification criteria before holdout validation is opened.

This is not a solved protein-folding claim. It does not run MD and does not use coordinates before prediction sealing.

## Panel

The panel contains 24 targets:

- 4 known anchors: KcsA, XCL1, alpha-synuclein, and 4AKE adenylate kinase.
- 4 membrane/channel/transporter proteins.
- 4 soluble compact single-domain proteins.
- 4 weak or shallow evolutionary-information targets.
- 4 disordered or partially disordered proteins.
- 4 multistate, allosteric, metamorphic, folding-upon-binding, or conformational-switch proteins.

## Seal Rule

Prediction must be sealed before validation:

- `prediction_timestamp`
- `prediction_hash`
- `sealed_prediction_packet.json`
- `prediction_inputs_manifest.json`
- `blocked_inputs_manifest.json`

Holdout manifests are written only after sealed prediction packets exist.

## Prediction Inputs

Allowed:

- sequence-level and annotation-level evidence
- UniProt feature/function text
- Pfam/InterPro signatures
- DisProt disorder annotations
- MSA/evolutionary metadata if available
- literature-derived state/function annotations
- non-coordinate experimental state evidence

Blocked:

- PDB coordinates before sealing
- coordinate-derived contacts
- native contact maps
- AlphaFold/ESMFold/RoseTTAFold coordinates
- prior V33/V34 coordinate-derived KcsA CSVs
- runtime reports as biological evidence
- validation holdouts before prediction sealing
- answer key or class labels during assignment

## Metrics

V42 reports:

- mechanism-class accuracy
- hard-class precision and recall
- operator-region support rate
- perturbation-prediction support rate
- low-resolution structure/ensemble support rate
- false hard grammar promotion count
- false single-fold promotion count
- clean abstain count
- leakage count

## Baselines

Required baselines:

- `random_class_baseline`
- `annotation_keyword_baseline`
- `majority_class_baseline`

Baselines are scored separately and cannot be counted as mechanism grammar.

## Required Controls

1. Prediction packets sealed before holdout validation.
2. Holdout files unavailable to prediction function.
3. PDB coordinates supplied before sealing are blocked.
4. Internal runtime reports supplied as biological evidence are blocked.
5. V33/V34 coordinate KcsA CSVs supplied as prediction evidence are blocked.
6. Answer-key leakage attempt is blocked.
7. Target-name-only assignment is blocked.
8. Generic annotation-only evidence cannot produce high-confidence hard grammar.
9. IDP forced into compact single-fold grammar is blocked.
10. Membrane generic-only evidence cannot become KcsA-like pore grammar without selectivity/filter/signature evidence.
11. Metamorphic/multistate claim requires state-separated evidence.
12. Random baseline does not count as mechanism.
13. Annotation-keyword baseline is reported separately.
14. Contradicting holdout marks the target failed, not silently rescued.
15. Fewer than 12 real targets clean-abstains as insufficient panel size.

## Status

Safe pass:

`V42_DE_NOVO_MECHANISM_LANGUAGE_CHALLENGE_PASSED_CLAIM_DISABLED`

The pass remains claim-disabled and does not set `folding_problem_solved`.
