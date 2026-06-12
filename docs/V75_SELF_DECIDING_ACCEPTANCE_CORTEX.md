# V75 Self-Deciding Acceptance Cortex

V75 replaced the static zero-failed-accepted gate idea with a self-decision judge.

The batch used 200 targets:

```text
33 V74 multidomain/allostery failure replays
65 E69 multidomain/allostery positives
90 protected sentinels
12 missing-word candidate clean-abstain controls
```

## Result

| metric | value |
| --- | ---: |
| targets_total | 200 |
| accepted_count | 187 |
| accepted_supported | 187 |
| clean_abstain_supported | 12 |
| blocked_supported | 0 |
| failed_accepted | 0 |
| accepted_accuracy | 1.0 |
| coverage | 0.935 |
| unknown_word_abstentions | 12 |
| sentinel_regressions | 0 |
| controls_passed | true |
| accepted_operator_basis_stable | 187 |
| accepted_cross_view_bound | 187 |
| accepted_temporal_bound | 187 |
| real_physical_calibration_rows | 8 |
| physical_basis_claim_allowed | false |

## Cortex Proof

The V75 packet now records:

```text
mechanism_competition
internal_consensus
dominance_law
cross_view_binding
evidence_masking
wrong_grammar_challenge
counterfactual_sequence_controls
operator_basis_stability_probe
temporal_binding_probe
physical_grounding_gate
contradictions
missing_word_candidate
final_self_decision
```

Accepted predictions were all supported. Candidate words that E69 does not yet know cleanly abstained instead of becoming failed accepted predictions.

One multidomain positive remained a conservative abstain because the endogenous operator-basis probe found the selected mechanism sensitive to assigning the observed operator strengths differently. That lowered coverage without creating a failed accepted prediction.

The acceptance law is:

```text
single_dominant_learned_mechanism_bound_across_views
```

Candidate grammars are intentionally not implemented in V75:

```text
clean_abstain_until_revision_implements_grammar
```

## Physical Calibration

V75 adds real physical calibration inputs from `data/folding_real_coordinate_visual_8.locked.json`.

```text
real_physical_calibration_inputs_used = true
source_coordinate_database = RCSB_PDB
row_count = 8
calibration_input_type = coordinate_derived_native_contact_summary
target_native_contacts_used_before_prediction = false
coordinate_truth_used_as_prediction_input = false
physical_grounding_status = real_physical_calibration_inputs_loaded_openmm_available_but_target_physical_execution_not_run
```

This means the engine now knows about real coordinate-derived calibration context, but it still refuses a physical basis claim until target-specific physical topology, environment, force-field protocol, and independent dynamic holdouts are executed.

## Missing Words

| missing_word_candidate | clean_abstain_supported |
| --- | ---: |
| `disulfide_secretory_redox_context` | 5 |
| `coiled_coil_register` | 4 |
| `repeat_solenoid_topology` | 3 |

## Next

```text
E70_SECRETORY_DISULFIDE_REDOX_TOPOLOGY_GRAMMAR
V76_SECRETORY_DISULFIDE_REPAIR_PANEL_200
```
