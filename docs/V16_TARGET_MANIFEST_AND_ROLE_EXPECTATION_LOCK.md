# V16 Target Manifest and Role-Expectation Lock

V16 is a locked transfer manifest, not a tuning panel and not a proof that protein folding is solved.

The V15 panel locked a unified, dynamic, role-aware evidence grammar across three protein object types:

- `1UBQ`: single-domain compact object
- `1CLL`: multi-domain composite object
- `4AKE`: domain-hinge closure object with D4 domain-core evidence

V16 asks a different question:

> Can the same locked grammar read new pressure regimes without changing the grammar?

## V16 pressure classes

The manifest locks three targets before any new data preflight or MD:

| Target | Pressure class | Expected role class |
|---|---|---|
| `p53_TAD_MDM2` | disorder + partner-induced binding | `disorder_partner_induced_object` |
| `KcsA` | membrane + pore + oligomer environment | `membrane_pore_oligomer_object` |
| `XCL1_lymphotactin` | metamorphic fold switching | `metamorphic_switch_object` |

## Locked policy

V16 manifest lock requires:

- no new MD at lock time
- no data download at lock time
- no target-specific threshold tuning
- no grammar changes
- no fixed residue-distance cutoff
- native metrics not used for selection
- `claim_allowed=false`

## Valid pass statuses are role-aware

V16 success is not `folded/not folded`.

Examples:

- p53 isolated TAD: `clean_abstain_no_autonomous_core`
- p53–MDM2: `partner_induced_interface_or_helix_signal_found`
- KcsA: `membrane_pore_roles_detected_without_soluble_core_misclassification`
- XCL1: `multiple_state_roles_supported_or_clean_state_specific_abstain`

Clean abstain is valid when role evidence is insufficient.
MD or target-specific data preflight should happen only after this manifest is locked.
