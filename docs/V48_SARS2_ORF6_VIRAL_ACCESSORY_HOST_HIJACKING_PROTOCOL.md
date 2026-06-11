# V48 SARS-CoV-2 ORF6 Viral Accessory Host-Hijacking Protocol

V48 tests a fourth hard mechanism class: viral short-region host hijacking.

The target is SARS-CoV-2 ORF6 / NS6, a 61-residue accessory protein. The scientific question is not a stable globular fold. The question is whether a sealed packet can identify a C-terminal host-interface grammar involving RAE1/NUP98 binding, nuclear transport disruption, interferon antagonism, and localization context without using coordinates or generic viral-protein shortcuts.

## Scope

- Target: SARS-CoV-2 ORF6 / accessory protein 6 / ns6
- UniProt: P0DTC6 / NS6_SARS2
- Region: full length 1-61, with C-terminal focus 38-61 and motif focus 50-61
- Mechanism class: `viral_accessory_short_region_host_hijacking_disorder_interface`
- Runtime claim: no MD, no coordinate prediction, no universal folding-solved claim

## Allowed Prediction Inputs

- ORF6 sequence
- UniProt non-coordinate function, interaction, localization, and mutagenesis text
- literature-derived host-interaction labels
- short-region motif or MoRF-like annotations
- non-coordinate perturbation/function evidence
- host-interaction labels for RAE1/NUP98, nuclear transport, and IFN antagonism

## Blocked Prediction Inputs

- PDB/mmCIF coordinates before sealing
- PDB-derived contacts
- native contact maps
- AlphaFold/ESMFold/RoseTTAFold coordinates
- ORF6 peptide-complex coordinate templates before sealing
- internal runtime reports as biological evidence
- validation holdouts before sealing
- target-name-only or generic viral accessory assignment
- swapped ORF8/ORF3a/NSP evidence as ORF6 validation

## Sealing Rule

V48 writes `sealed_prediction_packet.json` before any holdout validation.

The sealed packet includes:

- `prediction_hash`
- `prediction_timestamp`
- `prediction_inputs_manifest.json`
- `blocked_inputs_manifest.json`
- `no_holdout_access_before_hash: true`

Only after the hash exists does V48 write post-seal holdouts and validation scores.

## Required Operators

- `C_terminal_host_interaction_operator`
- `RAE1_NUP98_binding_context_operator`
- `nuclear_transport_disruption_operator`
- `interferon_antagonism_context_operator`
- `short_linear_motif_or_MoRF_operator`
- `disorder_to_interface_operator`
- `localization_context_operator`
- `no_globular_single_fold_operator`

## Forbidden Operators

- `compact_single_native_fold_operator`
- `generic_viral_accessory_annotation_only_operator`
- `membrane_channel_operator`
- `IDP_phase_separation_operator`
- `metamorphic_alpha_beta_fold_switch_operator`
- `solved_atomic_structure_operator`
- `coordinate_contact_operator`
- `AlphaFold_confidence_proxy_operator`

## Required Perturbation Logic

V48 must emit at least 10 explicit perturbation predictions. The required classes are:

- C-terminal deletion or disruption weakens host-hijacking grammar
- RAE1/NUP98 interaction evidence localizes the operator to the C-terminal region
- generic viral accessory annotation alone fails
- nuclear transport disruption is downstream of host-interaction grammar
- IFN antagonism is functional context, not atomic fold proof
- short-motif or MoRF-like evidence supports interface grammar, not stable globular fold
- compact single-fold forcing is blocked
- ORF8/ORF3a/NSP evidence does not validate ORF6-specific predictions
- removing RAE1/NUP98 evidence weakens the packet
- removing C-terminal/motif evidence weakens the packet
- C-terminal motif mutation evidence supports or falsifies the predicted operator

## Pass Conditions

V48 passes only if:

- prediction is sealed before holdout validation
- coordinate and internal-runtime leakage counts are zero
- mechanism class is exactly `viral_accessory_short_region_host_hijacking_disorder_interface`
- at least 7 operator buckets are produced
- at least 10 perturbation predictions are produced
- at least 7 perturbations are supported or partially supported by post-seal holdouts
- `host_interaction_support_rate >= 0.6`
- `functional_consequence_support_rate >= 0.6`
- compact globular-fold grammar is rejected
- generic viral annotation-only grammar is rejected
- contradictions are at most 2
- all leakage and shortcut controls pass
- claim boundary remains honest

V48 may set:

```text
live_viral_host_hijacking_solution_packet = true
protein_folding_solved_candidate_strengthened = true
```

only when the pass conditions are met.

V48 must never set:

```text
folding_problem_solved = true
```

## Honest Claim Boundary

Allowed wording:

`We have a sealed live solution packet for SARS-CoV-2 ORF6 as a short-region viral host-hijacking and disorder-interface mechanism. This is not a universal protein-folding solved claim, an atomic ORF6-host complex claim, or a stable globular-fold claim.`

Forbidden wording:

- universal protein folding is solved
- ORF6 is a solved compact globular fold
- ORF6-host atomic coordinates were predicted de novo
- generic viral accessory annotation solves ORF6
- IFN antagonism proves an atomic fold
- external review is unnecessary
