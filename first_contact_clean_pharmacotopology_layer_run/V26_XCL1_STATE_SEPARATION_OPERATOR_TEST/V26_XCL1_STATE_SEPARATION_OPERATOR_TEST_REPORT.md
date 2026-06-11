# V26 XCL1 State-Separation Operator Test

Status: `V26_XCL1_STATE_SEPARATION_OPERATOR_TEST_PASSED_CLAIM_DISABLED`
Mechanism operator: `state_separation`
State-separation operator passed: `True`
Positive pressure evidence: `True`
Positive folding evidence: `False`
Claim allowed: `False`
New MD executed: `False`

## Locked interpretation
V26 tests the repeated mechanism operator `state_separation` on XCL1. State A and state B remain separate; mixed-state fake core, single-fold forcing, fold-switch claim, native-metric selection, and MD remain forbidden. This is mechanism-operator evidence, not positive folding evidence.

## Next decision
{
  "claim_allowed": false,
  "decision_status": "V26_OPERATOR_TEST_LOCKED_NEXT_EVIDENCE_ACQUISITION",
  "kind": "V26_NEXT_MECHANISM_DECISION_v0",
  "new_MD_allowed": false,
  "new_MD_allowed_policy": "only_after_condition_labels_and_state_specific_external_constraints_are_locked_and_false_win_risk_is_low",
  "new_MD_recommended": false,
  "parallel_MD_paths_allowed": false,
  "positive_folding_evidence_targets": [],
  "reason": "XCL1 state-separation operator is clean and false-win risk remains controlled; next step is condition/coupling acquisition, not MD",
  "required_before_any_MD": [
    "condition_labels_or_state_context_annotations",
    "state_specific_external_couplings_or_constraints_if_available",
    "explicit_mixed_state_leakage_guard_preserved"
  ],
  "selected_next_panel": "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION"
}
