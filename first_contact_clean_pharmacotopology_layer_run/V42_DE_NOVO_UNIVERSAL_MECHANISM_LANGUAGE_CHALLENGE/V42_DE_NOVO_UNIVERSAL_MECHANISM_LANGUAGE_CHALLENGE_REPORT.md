# V42 De Novo Universal Mechanism Language Challenge

Status: `V42_DE_NOVO_MECHANISM_LANGUAGE_CHALLENGE_PASSED_CLAIM_DISABLED`
Panel target count: `24`
Sealed predictions: `24`
Mechanism-class accuracy: `1.000`
Operator-region support rate: `1.000`
Perturbation support rate: `1.000`
Controls: `15` / `15`

## Baselines
- `random_class_baseline`: `0.125`
- `annotation_keyword_baseline`: `0.792`
- `majority_class_baseline`: `0.208`

## Target Results
- `TARGET_001` `KcsA`: `membrane_pore_filter_or_transport` -> `supported`
- `TARGET_002` `XCL1_lymphotactin`: `metamorphic_or_multistate_switch` -> `supported`
- `TARGET_003` `alpha_synuclein_SNCA`: `intrinsic_disorder_contextual_ensemble` -> `supported`
- `TARGET_004` `4AKE_adenylate_kinase`: `compact_single_fold_core_closure` -> `supported`
- `TARGET_005` `bacteriorhodopsin_BR`: `membrane_pore_filter_or_transport` -> `supported`
- `TARGET_006` `AQP1_aquaporin`: `membrane_pore_filter_or_transport` -> `supported`
- `TARGET_007` `LacY_lactose_permease`: `membrane_pore_filter_or_transport` -> `supported`
- `TARGET_008` `CFTR_NBD1`: `metamorphic_or_multistate_switch` -> `supported`
- `TARGET_009` `ubiquitin_1UBQ`: `compact_single_fold_core_closure` -> `supported`
- `TARGET_010` `lysozyme_HEWL`: `compact_single_fold_core_closure` -> `supported`
- `TARGET_011` `myoglobin_Mb`: `compact_single_fold_core_closure` -> `supported`
- `TARGET_012` `barnase`: `compact_single_fold_core_closure` -> `supported`
- `TARGET_013` `SARS_CoV_2_ORF6`: `weak_evolutionary_information_abstain_or_low_confidence` -> `clean_abstain_supported`
- `TARGET_014` `SARS_CoV_2_ORF8`: `weak_evolutionary_information_abstain_or_low_confidence` -> `clean_abstain_supported`
- `TARGET_015` `p53_TAD`: `weak_evolutionary_information_abstain_or_low_confidence` -> `clean_abstain_supported`
- `TARGET_016` `Ebola_VP35_IID`: `weak_evolutionary_information_abstain_or_low_confidence` -> `clean_abstain_supported`
- `TARGET_017` `tau_K18`: `intrinsic_disorder_contextual_ensemble` -> `supported`
- `TARGET_018` `FUS_low_complexity_domain`: `intrinsic_disorder_contextual_ensemble` -> `supported`
- `TARGET_019` `TDP43_C_terminal_LCD`: `intrinsic_disorder_contextual_ensemble` -> `supported`
- `TARGET_020` `hnRNPA1_LCD`: `intrinsic_disorder_contextual_ensemble` -> `supported`
- `TARGET_021` `RfaH_CTD`: `metamorphic_or_multistate_switch` -> `supported`
- `TARGET_022` `Mad2`: `metamorphic_or_multistate_switch` -> `supported`
- `TARGET_023` `KaiB`: `metamorphic_or_multistate_switch` -> `supported`
- `TARGET_024` `calmodulin`: `metamorphic_or_multistate_switch` -> `supported`

## Plain English Interpretation
V42 is a sealed de novo mechanism-language challenge: it predicts mechanism class, operator regions, perturbation pressure, and low-resolution ensemble consequences before opening holdouts. Passing supports an attack path toward a universal mechanism language, but it remains claim-disabled and does not solve protein folding.
