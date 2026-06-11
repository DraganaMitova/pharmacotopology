# V31 Constraint-Backed Operator Readout Preflight

V31 is a provenance firewall after V30. V30 scans broadly for local candidate files and can list generated internal reports beside possible external constraints. V31 classifies every candidate file by provenance and assigns an allowed use before any V32 readout.

## Source classes

- `real_external_constraint_or_coupling`: allowed for constraint-backed operator readout.
- `real_external_alignment_source`: allowed only for deriving constraints in a later preflight, not as a direct constraint readout.
- `annotation_only_external_context`: allowed for role context only.
- `external_structure_source`: allowed for structure context or validation only.
- `generated_internal_report`: allowed for audit only.
- unverified/unusable classes: excluded.

## Pass/abstain rule

V31 passes into V32 only when a selected V31 target has real external constraint/coupling evidence with clean provenance. If selected targets only have annotations or internal reports, V31 clean-abstains and recommends external constraint import.

## Hard locks

- `claim_allowed=false`
- `new_MD_allowed=false`
- `new_MD_recommended=false`
- generated reports cannot be evidence claims
- annotation-only context cannot be constraint claims
