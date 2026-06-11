# V19 KcsA Membrane/Pore Evidence Readout

Zero-MD, no tuning, no membrane-MD claim. This reads membrane/pore/chain-context evidence and keeps soluble compact-core and whole-channel fold claims forbidden.

Status: `V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED`
Positive pressure evidence found: `True`
Positive folding evidence found: `False`
Claim allowed: `False`
New MD executed: `False`

## Role evidence
Role buckets assigned: `['chain_or_interface_context_support', 'ion_filter_diagnostic_shell', 'membrane_environment_context', 'pore_selectivity_filter_core', 'transmembrane_helix_scaffold']`
Pore/filter annotation present: `True`
Transmembrane role present: `True`
Oligomer/interface context present: `True`
Soluble-core misclassification avoided: `True`

## Pore/filter readout
Motif probe: `{'filter_motif_detected': True, 'strong_TVGYG_motif_detected': True, 'motif_hits': [{'source': 'residue_sequences_by_chain', 'chain': 'C', 'motif': 'TVGYG', 'index': 53}, {'source': 'residue_sequences_by_chain', 'chain': 'C', 'motif': 'TxxYG_diagnostic', 'index': 53, 'window': 'TVGYG'}, {'source': 'seqres_sequences_by_chain', 'chain': 'C', 'motif': 'TVGYG', 'index': 74}, {'source': 'seqres_sequences_by_chain', 'chain': 'C', 'motif': 'TxxYG_diagnostic', 'index': 74, 'window': 'TVGYG'}], 'motif_policy': 'sequence_or_seqres_report_only_no_single_fixed_success_threshold', 'selection_threshold_used': False}`
Annotation probe: `{'potassium_channel_annotation': True, 'pore_annotation': True, 'selectivity_filter_annotation': False, 'membrane_annotation': True, 'annotation_policy': 'pdb_header_keyword_report_only_no_native_selection', 'selection_threshold_used': False}`
Potassium ion probe: `{'potassium_ion_count': 7, 'potassium_ion_records_sample': [{'atom_name': 'K', 'residue_name': 'K', 'chain': 'C', 'residue_number': '3001', 'xyz': (155.336, 155.342, -30.553)}, {'atom_name': 'K', 'residue_name': 'K', 'chain': 'C', 'residue_number': '3002', 'xyz': (155.331, 155.331, -33.953)}, {'atom_name': 'K', 'residue_name': 'K', 'chain': 'C', 'residue_number': '3003', 'xyz': (155.341, 155.323, -37.162)}, {'atom_name': 'K', 'residue_name': 'K', 'chain': 'C', 'residue_number': '3004', 'xyz': (155.33, 155.324, -40.505)}, {'atom_name': 'K', 'residue_name': 'K', 'chain': 'C', 'residue_number': '3005', 'xyz': (155.327, 155.335, -47.577)}, {'atom_name': 'K', 'residue_name': 'K', 'chain': 'C', 'residue_number': '3006', 'xyz': (155.339, 155.327, -22.975)}, {'atom_name': 'K', 'residue_name': 'K', 'chain': 'C', 'residue_number': '3007', 'xyz': (155.343, 155.33, -26.017)}], 'ion_filter_diagnostic_shell_present': True, 'ion_probe_policy': 'ion_context_report_only_no_whole_fold_claim', 'selection_threshold_used': False}`

## Helix/interface readout
Helix records: `{'helix_record_count': 10, 'helix_records_by_chain': {'A': 4, 'B': 3, 'C': 3}, 'transmembrane_helix_scaffold_context_present': True, 'helix_policy': 'pdb_helix_records_report_only_no_membrane_MD_claim', 'selection_threshold_used': False}`
Chain interface readout: `{'chain_ca_counts': {'A': 219, 'B': 212, 'C': 103}, 'chain_count_with_ca': 3, 'multi_radius_contact_counts_by_chain_pair': {'A-B': {'6.0': 6, '8.0': 69, '10.0': 254}, 'A-C': {'6.0': 4, '8.0': 12, '10.0': 45}, 'B-C': {'6.0': 6, '8.0': 21, '10.0': 44}}, 'min_ca_distance_by_chain_pair': {'A-B': 5.046, 'A-C': 4.966, 'B-C': 4.847}, 'interface_pairs_present': ['A-B', 'A-C', 'B-C'], 'interface_context_present': True, 'contact_probe_radii_angstrom_report_only': [6.0, 8.0, 10.0], 'contact_probe_policy': 'multi_radius_report_only_no_single_fixed_selection_threshold', 'selection_threshold_used': False}`

## Guards
Forbidden misclassification violations: `[]`
Missing evidence: `['biological_tetramer_assembly_annotation_for_tetramer_claim', 'external_couplings_if_available']`
