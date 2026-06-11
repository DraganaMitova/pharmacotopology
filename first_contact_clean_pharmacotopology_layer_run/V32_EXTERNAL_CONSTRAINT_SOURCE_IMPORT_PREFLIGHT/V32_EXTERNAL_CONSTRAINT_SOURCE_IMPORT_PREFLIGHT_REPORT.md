# V32 External Constraint Source Import Preflight

Status: `V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED`
Selected V31 targets: `['XCL1_lymphotactin', 'KcsA']`
Selected V33 target: `KcsA`
Claim allowed: `False`
New MD allowed: `False`

## Locked interpretation
V32 is an acquisition/import preflight, not a folding win. It requires real external constraints for XCL1 or KcsA before any constraint-backed operator readout. Generated internal reports remain audit-only; annotations remain context-only; no MD and no claim are allowed.

## Target rows
### XCL1_lymphotactin
Status: `no_import_rows_for_target`
Ready for V33: `False`
Missing: `['real_external_constraint_import_rows']`
Valid real constraints: `0`
Invalid/excluded rows: `0`

### KcsA
Status: `target_ready_for_V33_constraint_backed_readout`
Ready for V33: `True`
Missing: `[]`
Valid real constraints: `2`
Invalid/excluded rows: `0`

## Next decision
{
  "claim_allowed": false,
  "decision_status": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED",
  "kind": "V32_NEXT_CONSTRAINT_BACKED_OPERATOR_READOUT_DECISION_v0",
  "new_MD_allowed": false,
  "new_MD_recommended": false,
  "next_action": "run_V33_constraint_backed_operator_readout_for_auto_selected_target_no_MD",
  "reason": "V32 imports only user-provided or locally present real external constraint sources. It does not treat V30/V31 runtime reports as evidence and does not run MD.",
  "selected_V33_panel": "V33_CONSTRAINT_BACKED_OPERATOR_READOUT",
  "selected_V33_target": "KcsA"
}
