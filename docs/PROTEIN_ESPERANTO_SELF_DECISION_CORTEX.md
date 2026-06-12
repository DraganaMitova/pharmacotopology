# Protein Esperanto Self-Decision Cortex

Starting with V75, acceptance is judged by self-consistency, not by a static unresolved-word gate.

The judge asks each target:

```text
Do mechanism readouts agree?
Do masking probes preserve the decision?
Do forced wrong grammars fail?
Do counterfactual sequence controls weaken the signal?
Do endogenous operator-basis probes preserve the qualitative mechanism signature?
Do simulated timepoint readouts bind to the selected trajectory instead of only the final state?
Does the physical grounding gate know whether real calibration inputs exist?
Are evidence, sequence, operators, and trajectory bound rather than independent labels?
Are there unresolved contradictions?
Is a missing Esperanto word competing with learned grammars?
```

Allowed final decisions:

| decision | meaning |
| --- | --- |
| `accepted` | one learned mechanism dominates the internal probes |
| `accepted_with_caution` | reserved for future cautious acceptance |
| `clean_abstain_missing_word` | a candidate missing word competes with learned grammars |
| `clean_abstain_conflict` | internal contradictions remain unresolved |
| `clean_abstain_low_internal_consensus` | no learned mechanism dominates |
| `blocked_for_leakage` | coordinate/internal/holdout leakage crossed the seal |

Batch success still requires:

```text
failed_accepted = 0
accepted_accuracy = 1.0
accepted_operator_basis_stable = accepted_count
accepted_cross_view_bound = accepted_count
accepted_temporal_bound = accepted_count
physical_basis_claim_allowed = false
```

Coverage may drop. That is acceptable when the engine is honestly saying that it does not know the word yet.

The operator-basis probe does not use a fixed coefficient perturbation range. It permutes the observed operator-strength multiset and maps single-operator ablations, then lets the self-decision judge abstain if the selected mechanism depends on a fragile assignment.

The physical gate now accepts real calibration inputs from locked RCSB coordinate rows, records their hash, and still blocks physical claims unless target-specific physical execution and independent dynamic holdouts exist.

## V75 Candidate Missing Words

Candidate grammars are not implemented in V75. Their acceptance role is:

```text
clean_abstain_until_revision_implements_grammar
```

The V75 cortex mined these clean abstentions:

| candidate | count |
| --- | ---: |
| `disulfide_secretory_redox_context` | 5 |
| `coiled_coil_register` | 4 |
| `repeat_solenoid_topology` | 3 |

The next repair target is `disulfide_secretory_redox_context`.
