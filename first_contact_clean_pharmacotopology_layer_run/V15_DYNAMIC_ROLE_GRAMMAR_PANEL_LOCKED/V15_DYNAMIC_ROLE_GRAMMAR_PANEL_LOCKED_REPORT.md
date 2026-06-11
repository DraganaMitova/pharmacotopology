# V15 Dynamic Role Grammar Panel Locked

This milestone is claim-disabled. It freezes the current dynamic grammar panel without claiming universal protein folding.

Lock status: `V15_THREE_PROTEIN_DYNAMIC_GRAMMAR_PANEL_LOCKED`
Source global status: `dynamic_separation_grammar_coherent_across_three_object_types_claim_disabled`
Claim allowed: `False`

## Locked interpretation

V15 demonstrates a unified, dynamic, role-aware evidence grammar across three protein object types: single-domain compact, multi-domain composite, and domain-hinge closure. It does not claim universal protein folding, does not claim full hinge/closure recovery for 4AKE or 1CLL, does not use a fixed residue-distance cutoff, and keeps claim_allowed=false.

## Locked rows
### 4AKE
- Artifact status: `present_machine_readable_4ake_role_artifact`
- Target role: `domain_hinge_closure_object`
- Final status: `domain_hinge_object_machine_readable_role_evidence_bridged;claim_allowed=false`
- Claim allowed: `False`

### 1UBQ
- Artifact status: `present`
- Target role: `single_domain_compact`
- Final status: `single_domain_compact_signal_found_with_dynamic_separation_context;claim_allowed=false`
- Claim allowed: `False`

### 1CLL
- Artifact status: `present`
- Target role: `multi_domain_composite`
- Final status: `multi_domain_composite_domain_core_signal_found_with_dynamic_separation_context;interdomain_hinge_not_yet_proven;claim_allowed=false`
- Claim allowed: `False`
