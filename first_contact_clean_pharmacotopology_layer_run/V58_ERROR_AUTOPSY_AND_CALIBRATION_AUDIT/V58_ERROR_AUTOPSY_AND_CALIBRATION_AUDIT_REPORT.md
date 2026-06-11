# V58 Error Autopsy And Calibration Audit

Status: `V58_ERROR_AUTOPSY_AND_CALIBRATION_AUDIT_COMPLETED_REVIEW_REQUIRED`
Targets: `20`
Raw supported: `14`
Raw failures preserved: `6`
Accepted count: `14`
Strict accepted count: `9`
Abstain recommended: `6`
Accepted accuracy: `1.0`
Strict accepted accuracy: `1.0`
Overall accuracy: `0.7`
Engine biology modified: `False`

## Failure Buckets
- `globular_vs_interface_or_oligomer_ambiguity`: `3`
- `globular_vs_switch_ambiguity`: `1`
- `insufficient_evidence_or_missing_context`: `2`

## Target Decisions
- `V58_10QM_4` selected `globular_closure` expected `globular_closure` decision `accepted` bucket `none`
- `V58_28ZU_37` selected `globular_closure` expected `globular_closure` decision `accepted` bucket `none`
- `V58_10QM_8` selected `membrane_multidomain_folding_proteostasis` expected `membrane_multidomain_folding_proteostasis` decision `accepted_with_caution` bucket `none`
- `V58_10PJ_1` selected `oligomerization_controlled_folding` expected `globular_closure` decision `abstain_recommended` bucket `globular_vs_interface_or_oligomer_ambiguity`
- `V58_10BT_2` selected `intrinsic_disorder_phase_separation` expected `intrinsic_disorder_phase_separation` decision `accepted` bucket `none`
- `V58_10QM_6` selected `globular_closure` expected `globular_closure` decision `accepted` bucket `none`
- `V58_28ZX_47` selected `globular_closure` expected `globular_closure` decision `accepted` bucket `none`
- `V58_10QM_17` selected `membrane_multidomain_folding_proteostasis` expected `membrane_multidomain_folding_proteostasis` decision `accepted_with_caution` bucket `none`
- `V58_12ZJ_1` selected `oligomerization_controlled_folding` expected `globular_closure` decision `abstain_recommended` bucket `globular_vs_interface_or_oligomer_ambiguity`
- `V58_10BT_1` selected `intrinsic_disorder_phase_separation` expected `intrinsic_disorder_phase_separation` decision `accepted` bucket `none`
- `V58_10QM_7` selected `globular_closure` expected `globular_closure` decision `accepted` bucket `none`
- `V58_10QM_30` selected `membrane_multidomain_folding_proteostasis` expected `membrane_multidomain_folding_proteostasis` decision `accepted_with_caution` bucket `none`
- `V58_24JD_1` selected `oligomerization_controlled_folding` expected `globular_closure` decision `abstain_recommended` bucket `globular_vs_interface_or_oligomer_ambiguity`
- `V58_10QM_3` selected `intrinsic_disorder_phase_separation` expected `intrinsic_disorder_phase_separation` decision `accepted_with_caution` bucket `none`
- `V58_10QM_10` selected `globular_closure` expected `globular_closure` decision `accepted` bucket `none`
- `V58_10QM_5` selected `membrane_multidomain_folding_proteostasis` expected `membrane_multidomain_folding_proteostasis` decision `accepted_with_caution` bucket `none`
- `V58_24SW_2` selected `metamorphic_fold_switching` expected `globular_closure` decision `abstain_recommended` bucket `globular_vs_switch_ambiguity`
- `V58_10QM_9` selected `insufficient_evidence_clean_abstain` expected `intrinsic_disorder_phase_separation` decision `abstain_recommended` bucket `insufficient_evidence_or_missing_context`
- `V58_10QM_11` selected `globular_closure` expected `globular_closure` decision `accepted` bucket `none`
- `V58_10YZ_1` selected `insufficient_evidence_clean_abstain` expected `membrane_multidomain_folding_proteostasis` decision `abstain_recommended` bucket `insufficient_evidence_or_missing_context`

## Boundary
No biological operators, mechanism classes, scoring controls, or engine weights were changed. The audit uses self-consistency, ambiguity labels, and failure buckets instead of tuning predictions toward universal acceptance.
