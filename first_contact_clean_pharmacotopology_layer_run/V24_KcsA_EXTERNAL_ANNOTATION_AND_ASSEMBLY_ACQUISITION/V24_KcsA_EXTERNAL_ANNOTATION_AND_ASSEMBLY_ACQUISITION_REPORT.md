# V24 KcsA External Annotation and Assembly Acquisition

Status: `V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION_PASSED_CLAIM_DISABLED`
Target role: `membrane_pore_oligomer_object`
Positive pressure evidence found: `True`
Positive folding evidence found: `False`
Claim allowed: `False`
New MD executed: `False`
Membrane MD executed: `False`

## Locked interpretation
V24 locks KcsA external membrane/pore annotation and assembly-context evidence: pore/filter diagnostic, TVGYG motif, potassium-ion shell, transmembrane helix scaffold, membrane context, and chain/interface context are preserved while soluble-core, tetramer-solved, whole-channel fold, membrane-MD, fixed-threshold, and native-metric claims remain forbidden. External couplings are not required for V24 but are marked as required/conditional for future V25 coupling/contact evidence tests.

## Evidence
Available: `['KcsA_external_annotation_source_references_locked', 'pore_selectivity_filter_annotation_or_diagnostic', 'TVGYG_selectivity_filter_motif', 'K_plus_ion_diagnostic_shell', 'transmembrane_helix_scaffold_context', 'membrane_environment_context', 'chain_or_interface_context', 'biological_assembly_context_reference_locked_for_context_only', 'leakage_guard_against_soluble_core_misread']`
Missing: `['external_couplings_if_available', 'expanded_biological_assembly_coordinate_or_author_defined_assembly_model_for_tetramer_claim']`
Failed checks: `[]`

## Next decision
`{'kind': 'V24_KcsA_NEXT_DECISION_v0', 'decision_status': 'V25_READY_FOR_KcsA_COUPLING_AND_INTERFACE_EVIDENCE_READOUT', 'selected_next_panel': 'V25_KcsA_COUPLING_AND_INTERFACE_EVIDENCE_READOUT', 'reason': 'KcsA external membrane/pore annotation layer is locked; pore/filter, TM scaffold, ion shell, and chain/interface context are preserved; no tetramer/folding claim allowed', 'new_md_recommended': False, 'membrane_md_recommended': False, 'parallel_md_paths_allowed': False, 'claim_allowed': False}`
