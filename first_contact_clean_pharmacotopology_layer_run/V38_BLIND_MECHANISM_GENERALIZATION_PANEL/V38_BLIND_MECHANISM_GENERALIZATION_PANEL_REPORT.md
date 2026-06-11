# V38 Blind Mechanism Generalization Panel

Status: `V38_BLIND_MECHANISM_GENERALIZATION_PASSED_CLAIM_DISABLED`
Known positive accuracy: `3` / `3`
False positive hard-grammar count: `0`
False negative known-positive count: `0`
Controls: `12` / `12`
Answer key used for assignment: `False`
Target names masked: `True`

## Masked Target Table
| Masked ID | Assigned class | Strength | Abstain reason |
|---|---|---:|---|
| TARGET_001 | membrane_pore_filter_oligomeric_ion_selectivity | 4 |  |
| TARGET_002 | metamorphic_two_state_fold_switch | 4 |  |
| TARGET_003 | intrinsic_disorder_contextual_ensemble | 4 |  |
| TARGET_004 | other_membrane_or_transport_context | 3 |  |
| TARGET_005 | other_membrane_or_transport_context | 3 |  |
| TARGET_006 | soluble_single_or_contextual_fold_not_metamorphic | 3 |  |
| TARGET_007 | soluble_single_or_contextual_fold_not_metamorphic | 3 |  |
| TARGET_008 | soluble_single_or_contextual_fold_not_metamorphic | 3 |  |
| TARGET_009 | soluble_single_or_contextual_fold_not_metamorphic | 3 |  |

## Answer-Key Evaluation
| Masked ID | Target | Expected | Assigned |
|---|---|---|---|
| TARGET_001 | KcsA | membrane_pore_filter_oligomeric_ion_selectivity | membrane_pore_filter_oligomeric_ion_selectivity |
| TARGET_002 | XCL1_lymphotactin | metamorphic_two_state_fold_switch | metamorphic_two_state_fold_switch |
| TARGET_003 | alpha_synuclein_SNCA | intrinsic_disorder_contextual_ensemble | intrinsic_disorder_contextual_ensemble |
| TARGET_004 | bacteriorhodopsin | other_membrane_or_transport_context | other_membrane_or_transport_context |
| TARGET_005 | AQP1 | other_membrane_or_transport_context | other_membrane_or_transport_context |
| TARGET_006 | CXCL8 | soluble_single_or_contextual_fold_not_metamorphic | soluble_single_or_contextual_fold_not_metamorphic |
| TARGET_007 | ubiquitin | soluble_single_or_contextual_fold_not_metamorphic | soluble_single_or_contextual_fold_not_metamorphic |
| TARGET_008 | lysozyme | soluble_single_or_contextual_fold_not_metamorphic | soluble_single_or_contextual_fold_not_metamorphic |
| TARGET_009 | myoglobin | soluble_single_or_contextual_fold_not_metamorphic | soluble_single_or_contextual_fold_not_metamorphic |

## Controls
- `baseline_masked_panel_assignment`: `PASS`
- `unmasked_vs_masked_consistency_check`: `PASS`
- `swapped_known_positive_dossiers_detected`: `PASS`
- `membrane_decoy_not_kcsa_without_filter_ion`: `PASS`
- `chemokine_decoy_not_xcl1_without_two_state`: `PASS`
- `folded_decoy_not_snca_without_disorder`: `PASS`
- `generic_only_annotations_clean_abstain`: `PASS`
- `coordinate_derived_source_blocked`: `PASS`
- `internal_runtime_source_blocked`: `PASS`
- `placeholder_source_blocked`: `PASS`
- `target_name_only_assignment_attempt_blocked`: `PASS`
- `answer_key_leakage_attempt_blocked`: `PASS`

## Plain English Interpretation
V38 passing means the mechanism grammar survived masked target names and near-decoy pressure: the three known positives were recovered, and membrane/soluble/folded decoys were not promoted into the hard positive grammars. This is mechanism-class discrimination, not atom prediction or a folding-solved claim.
