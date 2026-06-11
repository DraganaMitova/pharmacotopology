# V23 p53 TAD/MDM2 External Evidence Contrast Test

Status: `V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST_PASSED_CLAIM_DISABLED`
Positive pressure evidence found: `True`
Positive folding evidence found: `False`
Claim allowed: `False`
New MD executed: `False`

## Contrast
Isolated TAD clean abstain preserved: `True`
Bound partner-induced evidence preserved: `True`
External disorder reference context found: `True`
Partner-bound interface reference context found: `True`

## Interface readout
`{'chain_ca_counts': {'A': 85, 'B': 13}, 'chain_pair': 'A-B', 'contact_probe_policy': 'multi_radius_report_only_no_single_fixed_selection_threshold', 'contact_probe_radii_angstrom_report_only': [6.0, 8.0, 10.0], 'interface_contact_evidence_present': True, 'large_chain': 'A', 'min_ca_distance': 4.884, 'multi_radius_contact_counts': {'10.0': 96, '6.0': 7, '8.0': 28}, 'nearest_ca_pair': {'distance_angstrom': 4.884, 'large_chain': 'A', 'large_residue_number': 50, 'small_chain': 'B', 'small_residue_number': 29}, 'selection_threshold_used': False, 'small_chain': 'B'}`

## Source references
- `DisProt_DP00086` — `disorder_annotation_reference` — https://disprot.org/DP00086
- `RCSB_1YCR` — `partner_bound_interface_structure_reference` — https://www.rcsb.org/structure/1YCR

## Missing but nonblocking for this test
- `raw_local_external_annotation_files_if_needed`
- `external_couplings_if_available`

## Forbidden misclassification violations
`[]`

## Interpretation
V23 locks the p53/MDM2 external-evidence contrast: external disorder-source context supports the isolated TAD abstain guard, while the 1YCR partner-bound complex preserves interface/helix role evidence. This is positive pressure-context evidence, not positive folding evidence.
