# V41 Protein Folding Claim Gate

Status: `V41_MECHANISM_CLAIM_ALLOWED_C5_BLOCKED`
Max allowed claim level: `C3_FALSIFIABLE_OPERATOR_MECHANISM`
Protein-folding-solved claim allowed: `False`
C5 claim allowed: `False`
Controls: `21` / `21`

## Allowed Claim
We have a claim-disabled, coordinate-free mechanism-grammar prototype that classifies hard protein-folding regimes and predicts perturbation-sensitive operator constraints for KcsA, XCL1, and alpha-synuclein, with masked-target and decoy controls, no MD, and no coordinate leakage. This is not solved protein folding, is not a de novo protein-structure predictor, and does not solve protein folding.

## Forbidden Claims
- we solved protein folding
- we can predict every protein
- we can predict de novo 3D structure from sequence
- we outperform AlphaFold
- KcsA/XCL1/SNCA prove universal folding
- coordinates were not needed to solve structure
- this is a validated drug-discovery engine

## Missing Evidence For C5
- broad prospective blind panel
- de novo coordinate or ensemble predictions
- independent experimental structure/function validation
- quantitative comparison against strong baselines
- failure modes documented across broad protein classes
- reproducibility certificate
- external review-ready report

## Scientist Attack Surface
- Only three focal mechanism exemplars are supported at C3: KcsA, XCL1, and alpha-synuclein.
- V38 is masked and decoy-controlled, but V41 does not add a new prospective unseen panel.
- The claim is coordinate-free mechanism grammar, not de novo 3D coordinate prediction.
- Static structure, folding pathway, dynamics, membrane context, and ensemble prediction remain distinct claims.
- External benchmark comparison against AlphaFold/ESMFold/RoseTTAFold is not present in V41.
- Perturbation pressure is literature/annotation-supported, not newly experimentally validated here.
- Drug-discovery utility, therapeutic validation, and clinical relevance are not established.
- C5 requires external blind validation and broad reproducibility evidence that is not present.

## Plain English Interpretation
V41 allows a bounded C3 mechanism-grammar claim and blocks C5. The project can say it has a coordinate-free, claim-disabled mechanism grammar with falsifiable perturbation-sensitive operators for the three audited regimes; it cannot say protein folding is solved, cannot claim de novo 3D structure prediction, and cannot claim superiority to AlphaFold.
