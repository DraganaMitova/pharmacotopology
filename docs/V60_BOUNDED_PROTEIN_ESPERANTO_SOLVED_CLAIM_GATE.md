# V60 Bounded Protein Esperanto Solved-Claim Gate

V60 is a claim reader, not a predictor.

It reads:

- V50-V56 engine suite
- V57 blind generalization
- V58 real-sequence gate
- V58 calibration/autopsy
- V59 real folding-process panel

Then it computes the allowed claim.

## Computed Outcome

```text
V60_BOUNDED_PROCESS_REPLICATION_CLAIM_PASSED_REVIEW_REQUIRED
```

Allowed claim:

```text
The protein folding problem is not universally solved here; however, after the V59-required bounded engine hardening, the current Protein Esperanto engine supports a bounded solution claim on accepted targets: it predicts folding regime, important regions, topology/observable behavior, and selected folding-process observables under source-separated validation.
```

Still forbidden:

- universal protein folding is solved
- all sequences are solved
- atomistic folding is solved
- AlphaFold is replaced
- unsupported targets can be answered instead of abstained
- the V59 engine revision can be hidden under a frozen-engine claim

## Claim Flags

```text
protein_esperanto_engine_claim_allowed = true
bounded_accepted_target_process_replication_supported = true
universal_folding_problem_solved = false
atomistic_folding_solved = false
all_sequence_prediction_solved = false
engine_revision_after_v58_required = true
```
