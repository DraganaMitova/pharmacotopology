# V32 External Constraint Source Import Preflight

Status: `V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED`
Selected V31 targets: `['XCL1_lymphotactin', 'KcsA']`
Selected V33 target: `None`
Claim allowed: `False`
New MD allowed: `False`

## Locked interpretation
V32 is an acquisition/import preflight, not a folding win. It requires real external constraints for XCL1 or KcsA before any constraint-backed operator readout. Generated internal reports remain audit-only; annotations remain context-only; no MD and no claim are allowed.

## Target rows
### XCL1_lymphotactin
Status: `partial_or_missing_state_specific_external_constraints`
Ready for V33: `False`
Missing: `['state_A_real_external_constraint', 'state_B_real_external_constraint']`
Valid real constraints: `0`
Invalid/excluded rows: `2`

### KcsA
Status: `partial_or_missing_kcsa_coupling_interface_constraints`
Ready for V33: `False`
Missing: `['pore_filter_or_external_coupling_constraint', 'assembly_or_chain_interface_constraint']`
Valid real constraints: `0`
Invalid/excluded rows: `2`

## Next decision
{
  "claim_allowed": false,
  "decision_status": "V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED",
  "kind": "V32_NEXT_CONSTRAINT_BACKED_OPERATOR_READOUT_DECISION_v0",
  "new_MD_allowed": false,
  "new_MD_recommended": false,
  "next_action": "place_real_external_constraint_files_under_data/external_constraints_and_fill_v32_import_manifest_then_rerun_V32",
  "reason": "V32 imports only user-provided or locally present real external constraint sources. It does not treat V30/V31 runtime reports as evidence and does not run MD.",
  "selected_V33_panel": null,
  "selected_V33_target": null
}
