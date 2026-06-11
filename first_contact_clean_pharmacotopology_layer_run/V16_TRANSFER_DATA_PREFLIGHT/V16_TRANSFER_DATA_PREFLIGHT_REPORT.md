# V16 Transfer Data Preflight

This is a data/material preflight only. It does not run MD, tune thresholds, or claim solved folding.

Status: `V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT`
Claim allowed: `False`
New MD executed: `False`
Ready targets: `['p53_TAD_MDM2', 'KcsA', 'XCL1_lymphotactin']`
Blocked targets: `[]`

## Targets
### p53_TAD_MDM2
- Pressure class: `disorder_partner_induced_binding`
- Expected role class: `disorder_partner_induced_object`
- Target material status: `ready_for_zero_md_role_readout`
- Role preflight status: `partner_complex_context_present_isolated_disorder_material_later_no_claim`
- Missing required material: `[]`

Required material:
- `p53_TAD_MDM2_complex_structure` / `1YCR`: `present_and_provenanced`
  - Failed checks: `[]`
  - CA atoms: `98`
  - Chains with CA: `2`

### KcsA
- Pressure class: `membrane_pore_oligomer_environment`
- Expected role class: `membrane_pore_oligomer_object`
- Target material status: `ready_for_zero_md_role_readout`
- Role preflight status: `membrane_pore_oligomer_context_present_membrane_annotation_later_no_claim`
- Missing required material: `[]`

Required material:
- `KcsA_Fab_complex_structure` / `1K4C`: `present_and_provenanced`
  - Failed checks: `[]`
  - CA atoms: `534`
  - Chains with CA: `3`

### XCL1_lymphotactin
- Pressure class: `metamorphic_fold_switching`
- Expected role class: `metamorphic_switch_object`
- Target material status: `ready_for_zero_md_role_readout`
- Role preflight status: `two_state_metamorphic_context_present_state_specific_readout_later_no_claim`
- Missing required material: `[]`

Required material:
- `XCL1_state_A_chemokine_like_structure` / `2HDM`: `present_and_provenanced`
  - Failed checks: `[]`
  - CA atoms: `1480`
  - Chains with CA: `1`
- `XCL1_state_B_alternative_conformation_structure` / `2JP1`: `present_and_provenanced`
  - Failed checks: `[]`
  - CA atoms: `2400`
  - Chains with CA: `2`
