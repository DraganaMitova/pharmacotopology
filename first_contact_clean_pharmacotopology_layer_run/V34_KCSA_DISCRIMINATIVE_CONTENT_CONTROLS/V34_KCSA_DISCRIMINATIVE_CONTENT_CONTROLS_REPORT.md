# V34 KcsA Discriminative Content Controls

Status: `V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_PASSED_CLAIM_DISABLED`
Passed controls: `8` / `8`
Claim allowed: `False`
New MD allowed: `False`
Folding solved: `False`

## Controls
### baseline_real_kcsa_content_signature_valid
Passed: `True`
Observed valid: `True`
Expected valid: `True`
Observed failures: `[]`
Expected failure any: `[]`
Reason: Real imported KcsA source rows must carry pore/filter identity, K+ coordination identity, and tetramer/interface identity.

### pore_residue_shuffle_breaks_canonical_filter_identity
Passed: `True`
Observed valid: `False`
Expected valid: `False`
Observed failures: `['missing_canonical_filter_residues_75_78']`
Expected failure any: `['missing_canonical_filter_residues_75_78']`
Reason: A residue-number randomized pore file must not still count as the KcsA canonical filter.

### pore_motif_relabel_breaks_TVGYG_identity
Passed: `True`
Observed valid: `False`
Expected valid: `False`
Observed failures: `['missing_TVGYG_motif_label']`
Expected failure any: `['missing_TVGYG_motif_label']`
Reason: A pore file without the TVGYG label must not pass as KcsA filter grammar.

### pore_ion_relabel_breaks_potassium_identity
Passed: `True`
Observed valid: `False`
Expected valid: `False`
Observed failures: `['missing_potassium_ion_identity']`
Expected failure any: `['missing_potassium_ion_identity']`
Reason: KcsA filter evidence must not survive if K+ identity is removed or relabeled.

### pore_constraint_class_relabel_breaks_filter_contact_identity
Passed: `True`
Observed valid: `False`
Expected valid: `False`
Observed failures: `['missing_pore_filter_potassium_coordination_class']`
Expected failure any: `['missing_pore_filter_potassium_coordination_class']`
Reason: Generic contacts must not be promoted into pore/filter potassium-coordination evidence.

### interface_same_chain_collapse_breaks_tetramer_context
Passed: `True`
Observed valid: `False`
Expected valid: `False`
Observed failures: `['insufficient_distinct_chain_pairs', 'insufficient_oligomer_chains']`
Expected failure any: `['insufficient_oligomer_chains', 'insufficient_distinct_chain_pairs']`
Reason: A collapsed same-chain interface file must not pass as tetramer/interface context.

### singleton_rows_do_not_define_operator_context
Passed: `True`
Observed valid: `False`
Expected valid: `False`
Observed failures: `['assembly_interface_too_few_rows', 'insufficient_distinct_chain_pairs', 'insufficient_oligomer_chains', 'missing_canonical_filter_residues_75_78', 'pore_filter_not_observed_across_four_chains', 'pore_filter_too_few_rows']`
Expected failure any: `['pore_filter_too_few_rows', 'assembly_interface_too_few_rows']`
Reason: One row is not enough to define a stable operator context for either pore or interface.

### swapped_pore_and_interface_schema_fails_content_validation
Passed: `True`
Observed valid: `False`
Expected valid: `False`
Observed failures: `['insufficient_distinct_chain_pairs', 'missing_TVGYG_motif_label', 'missing_assembly_interface_contact_class', 'missing_pore_filter_potassium_coordination_class', 'missing_potassium_ion_identity']`
Expected failure any: `['missing_TVGYG_motif_label', 'missing_assembly_interface_contact_class', 'insufficient_distinct_chain_pairs']`
Reason: Swapping pore/interface schemas must not pass merely because rows exist.

## Baseline content signature
```json
{
  "interface_chain_pair_count": 6,
  "interface_chain_pairs": [
    "A-B",
    "A-C",
    "A-D",
    "B-C",
    "B-D",
    "C-D"
  ],
  "interface_chains": [
    "A",
    "B",
    "C",
    "D"
  ],
  "interface_constraint_classes": [
    "ASSEMBLY_INTERFACE_HEAVY_ATOM_CONTACT"
  ],
  "interface_row_count": 186,
  "pore_chains": [
    "A",
    "B",
    "C",
    "D"
  ],
  "pore_constraint_classes": [
    "PORE_FILTER_POTASSIUM_COORDINATION_CONTACT"
  ],
  "pore_ion_names": [
    "K"
  ],
  "pore_motif_labels": [
    "TVGYG"
  ],
  "pore_residue_numbers": [
    75,
    76,
    77,
    78
  ],
  "pore_row_count": 20
}
```

## Locked interpretation
Passing V34 means the imported KcsA source files are not merely counted: damaged pore/filter identity, K+ identity, constraint-class identity, tetramer/interface identity, singleton rows, and swapped schemas fail. This is still not de-novo folding prediction and does not solve universal protein folding.
