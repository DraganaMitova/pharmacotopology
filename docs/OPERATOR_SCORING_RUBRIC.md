# Operator Scoring Rubric

V49 freezes operator-level scoring as pre-registered rows. A run does not earn "validated" wording unless its explicit predictions map to post-seal evidence under this schema.

## Required Prediction Row Fields

```text
prediction_id
target
mechanism_class
operator_bucket
region_or_state
predicted_effect
perturbation
expected_direction
confidence
prediction_source_ids
falsification_criteria
holdout_evidence_ids
score_label
score_reason
scoring_pre_registered
```

`scoring_pre_registered` must be `true` for rows that are eligible for support claims.

## Score Labels

Allowed `score_label` values:

```text
supported
partially_supported
contradicted
not_testable
blocked_for_leakage
```

Meaning:

- `supported`: the post-seal holdout supports the predicted direction and mechanism bucket.
- `partially_supported`: the holdout supports part of the prediction but leaves a material condition, region, or causal step unresolved.
- `contradicted`: the holdout opposes the predicted direction or mechanism bucket.
- `not_testable`: the available holdout does not test the prediction.
- `blocked_for_leakage`: the prediction cannot be scored because the evidence boundary was violated.

## Frozen CFTR F508del Mapping Requirement

For V46-style CFTR F508del claims, scoring must map each claim to a row:

- F508del weakens NBD1 stability.
- NBD1-only rescue is partial.
- Interdomain/interface correction is needed.
- Trafficking/proteostasis is maturation context, not atomic proof.
- Generic channel annotation is insufficient.
- Coordinate evidence before sealing blocks the claim.

The row-level result must be visible in the certificate or report, not only summarized as "validated."
