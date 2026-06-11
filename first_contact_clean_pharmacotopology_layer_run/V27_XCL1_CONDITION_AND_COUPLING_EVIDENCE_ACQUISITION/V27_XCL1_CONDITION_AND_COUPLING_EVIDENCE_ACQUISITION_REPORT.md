# V27 XCL1 Condition and Coupling Evidence Acquisition

Status: `V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION_PASSED_CLAIM_DISABLED`
Condition label context locked: `True`
State-specific couplings present: `False`
Positive pressure evidence: `True`
Positive folding evidence: `False`
Claim allowed: `False`
New MD executed: `False`

## Locked interpretation
V27 acquires the XCL1 condition/state-label layer and scans for state-specific couplings without inventing them. State A and state B remain separate; mixed-state fake core, single-fold forcing, fold-switch claim, MD, threshold tuning, and claim upgrade remain forbidden.

## Next decision
{
  "claim_allowed": false,
  "decision_status": "V28_TARGET_SELECTED",
  "kind": "V27_XCL1_NEXT_DECISION_v0",
  "new_MD_allowed": false,
  "new_MD_allowed_policy": "only_after_state_specific_external_constraints_are_locked_and_false_win_risk_is_low",
  "new_MD_recommended": false,
  "parallel_MD_paths_allowed": false,
  "positive_folding_evidence_targets": [],
  "reason": "XCL1 state-condition labels and mixed-state leakage guard are locked; state-specific couplings remain optional for V28 and required before any MD/contact claim",
  "required_before_any_MD": [
    "state_specific_external_couplings_or_constraints_if_available",
    "condition_labels_or_state_context_annotations_preserved",
    "mixed_state_leakage_guard_preserved"
  ],
  "selected_V28_target": "XCL1_lymphotactin",
  "selected_V28_test": "XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST",
  "selected_next_panel": "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST"
}
