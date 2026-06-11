# Negative Control Inventory

V49 makes negative and null controls README-visible. Future live packets should expose the applicable controls in their protocol, certificate, or report.

## Required Controls

- `random_sequence_control`: random sequence input must not create a target-specific packet.
- `shuffled_sequence_control`: shuffled sequence input must weaken or block operator assignment.
- `swapped_evidence_control`: evidence from another target must not validate target-specific predictions.
- `wrong_target_control`: wrong target labels must not pass as the live target.
- `generic_annotation_only_control`: generic function annotation alone must not create a full packet.
- `coordinate_leakage_control`: coordinate-derived sources before sealing must block prediction use.
- `internal_runtime_leakage_control`: internal reports and certificates must not become biological evidence.
- `forced_wrong_grammar_control`: forcing the wrong mechanism grammar must fail or cleanly abstain.
- `failed_prediction_not_repaired_after_holdout`: failed predictions must remain failed after holdout opening.
- `holdout_opened_before_seal_control`: holdout evidence opened before sealing must block the run.

## Control Interpretation

Passing these controls does not solve protein folding and does not prove an atomic structure. It only says the packet stayed source-separated, resisted obvious shortcuts, and preserved the claim boundary while making falsifiable operator-level predictions.
