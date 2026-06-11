# V18 p53 Context Contrast Lock

This layer closes the first p53/MDM2 pressure-evidence milestone without running MD and without upgrading any folding claim.

## Locked pieces

1. `V18_P53_PARTNER_INDUCED_EVIDENCE_LOCK`
   - freezes the V18 partner-bound evidence result
   - keeps `positive_pressure_evidence_found=true`
   - keeps `positive_folding_evidence_found=false`
   - keeps `claim_allowed=false`

2. `V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT`
   - preserves the bound-context partner-induced interface/helix evidence
   - keeps isolated p53 TAD autonomous compact-core/fold claims forbidden
   - treats missing isolated TAD evidence as a clean abstain context, not as a false positive

## Boundary

This is a pressure-context role test, not a universal protein folding claim.

The locked interpretation is:

```text
isolated p53 TAD -> clean abstain / no autonomous compact fold claim
p53 TAD + MDM2 -> partner-induced interface/helix role evidence
claim_allowed = false
```
