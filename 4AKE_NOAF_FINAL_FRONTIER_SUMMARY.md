# 4AKE no-AlphaFold final frontier summary

## What was added in this final pass

A new fallback branch was implemented after the iterative distance-geometry diffusion result:

```text
leave-one-out empirical sequence-contact prior
```

This branch trains only on the other locked coordinate rows:

```text
1MBN:A, 2LZM:A, 1TEN:A, 1CSP:A, 1TIM:A, 1PGA:A, 1CLL:A
```

It explicitly excludes 4AKE from training, predicts from 4AKE sequence-derived features, and attaches 4AKE native contacts only after selection for audit metrics.

No AlphaFold, no MSA, no GIF, no raw sequence persistence, no hanging full suite.

## Final audited branches

| branch | precision | recall | F1 | long-range recall | verdict |
|---|---:|---:|---:|---:|---|
| iterative DG diffusion | 0.181818 | 0.143113 | 0.160160 | 0.088083 | not solved |
| leave-one-out empirical sequence prior | 0.487528 | 0.384615 | 0.430000 | 0.000000 | not solved |
| empirical local + DG long-range fusion | 0.417266 | 0.415027 | 0.416143 | 0.088083 | not solved |

Solve boundary remains:

```text
precision >= 0.70
recall >= 0.70
long_range_recall >= 0.70
```

Final verdict:

```text
folding_problem_solved_after_native_audit = false
claim_rejection_reason = no_noaf_branch_reached_0_70_precision_0_70_recall_0_70_long_range_recall
```

## Scientific interpretation

The new empirical prior is the strongest no-AlphaFold branch by F1. It raises 4AKE from the earlier no-AF zone into a much better local-contact regime:

```text
DG diffusion F1: 0.160160
empirical prior F1: 0.430000
```

But it does this by capturing mostly local/short-range structural regularities. It does not recover the long-range tertiary contacts needed for a real fold:

```text
empirical prior long_range_recall = 0.000000
```

The fused branch adds DG long-range contacts back in, raising total recall to 0.415027, but long-range recall remains only 0.088083.

## Final conclusion

4AKE is still not cracked without AlphaFold/MSA/a real MSA-free learned structure source.

What is now proven by code and tests:

```text
hand physics alone: not enough
iterative distance-geometry diffusion: not enough
small leave-one-out empirical prior: improves strongly, but local-only
empirical + DG fusion: better recall, still not enough
```

The missing piece is now sharper than before:

```text
A global MSA-free learned distance-geometry prior is required.
```

That means ESMFold/OmegaFold/SPIRED-style output, or a new trained generative model with global tertiary geometry. More thresholds or hand rules are not enough for 4AKE.

## Validation

```text
7 passed, 2 skipped in 16.30s
```

The skipped tests are opt-in visual / heavier self-consistency paths. GIF generation remains off by default.
