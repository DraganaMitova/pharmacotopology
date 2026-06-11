# V58 Error Autopsy And Calibration Audit

V58.1/V58.2 reads the sealed V58 real-sequence results and explains the failures without changing the Protein Esperanto engine.

## Rule

No biological logic changes:

- no new operators
- no new mechanism classes
- no engine weight tuning
- no README update
- no failed prediction repair

The audit only adds self-decision labels:

```text
accepted
accepted_with_caution
abstain_recommended
```

## Purpose

V58 showed:

```text
14 / 20 raw support
6 failures preserved
```

The calibration audit asks:

```text
Could the system have avoided failures by saying "not enough signal" or "ambiguous"?
```

It does not try to force 20/20.

## Failure Buckets

The audit classifies failures into:

- `globular_vs_interface_or_oligomer_ambiguity`
- `globular_vs_switch_ambiguity`
- `insufficient_evidence_or_missing_context`
- `right_regime_wrong_topology_or_observable`
- `regime_ambiguity_detectable_pre_holdout`
- `wrong_regime_not_yet_detectable_by_current_confidence`

## Certificate Fields

The certificate reports:

- `accepted_accuracy`
- `strict_accepted_accuracy`
- `abstention_rate`
- `overall_accuracy`
- `failure_preserved`
- `engine_biology_modified = false`

## Boundary

An abstention is not a biological success. It is a reliability behavior. The useful claim is that accepted predictions were correct under the audit, while unsupported cases were preserved and marked for abstention instead of being silently repaired.
