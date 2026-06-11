# V16 Target Manifest and Role-Expectation Lock

This is not a new tuning panel and not a folding-solved claim.

Lock status: `V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCKED`
Claim allowed: `False`
New MD executed: `False`
Data preflight executed: `False`

## Locked interpretation

V16 tests whether the locked V15 dynamic role grammar can transfer from protein object types to new pressure regimes: disorder/binding, membrane/pore/oligomer, and fold-switching. It is not a new tuning panel, does not run MD at lock time, forbids fixed target-specific thresholds, does not use native metrics for selection, and keeps claim_allowed=false.

## Targets
### p53_TAD_MDM2
- Pressure class: `disorder_partner_induced_binding`
- Expected role class: `disorder_partner_induced_object`
- Clean abstain allowed: `True`
- Claim allowed: `False`

Allowed evidence roles:
- `clean_abstain_no_autonomous_core`
- `partner_induced_interface_or_helix_signal_found`
- `binding_stabilized_interface_core`
- `disorder_signal_report_only`

Forbidden misclassification:
- `isolated_TAD_as_soluble_compact_core`
- `partner_interface_as_autonomous_fold_claim`
- `binding_helix_as_universal_p53_fold`

### KcsA
- Pressure class: `membrane_pore_oligomer_environment`
- Expected role class: `membrane_pore_oligomer_object`
- Clean abstain allowed: `True`
- Claim allowed: `False`

Allowed evidence roles:
- `transmembrane_helix_scaffold`
- `pore_selectivity_filter_core`
- `tetramer_interface_support`
- `membrane_environment_pressure_context`
- `ion_filter_diagnostic_shell`
- `clean_abstain_missing_membrane_or_oligomer_context`

Forbidden misclassification:
- `membrane_pore_as_soluble_single_domain_compact_core`
- `ion_filter_geometry_as_whole_fold_claim`
- `tetramer_interface_as_monomer_autonomous_core`

### XCL1_lymphotactin
- Pressure class: `metamorphic_fold_switching`
- Expected role class: `metamorphic_switch_object`
- Clean abstain allowed: `True`
- Claim allowed: `False`

Allowed evidence roles:
- `state_A_chemokine_monomer_support`
- `state_B_beta_sandwich_dimer_support`
- `shared_local_support`
- `pressure_condition_switch_context`
- `clean_state_specific_abstain`

Forbidden misclassification:
- `forcing_single_canonical_fold`
- `mixing_two_states_into_one_false_core`
- `claiming_fold_switch_without_state_specific_support`
