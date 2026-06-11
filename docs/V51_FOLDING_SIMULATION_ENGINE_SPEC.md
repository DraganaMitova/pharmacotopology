# V51 Folding Simulation Engine Spec

V51 freezes the engine contract before implementation drift.

## Input Schema

- `sequence`: amino-acid sequence.
- `prediction_sources`: V49-classified source rows.
- `focus_regions`: optional named residue spans with no coordinates.
- `forced_grammar`: optional control input that must fail or abstain when wrong.

## State Vector Schema

Each residue state includes `position_index`, `residue_identity`, `segment_id`, `charge_mark`, `hydrophobic_mark`, `disorder_mark`, `secondary_propensity_mark`, `interface_mark`, `membrane_mark`, `operator_activations`, `current_state`, and `state_confidence`.

## Operator Schema

Each operator includes `operator`, `acts_on`, `activation_strength`, `activated_by_evidence_ids`, `pushes_toward`, `perturbation_should`, `falsified_by`, and `state_variable`.

## Simulation Cycle

```text
gate evidence
-> build sequence field
-> select or abstain mechanism grammar
-> build operator field
-> simulate coarse state trajectory
-> seal prediction hash
-> open holdout validation
```

## Trajectory Output Schema

Trajectory rows include `timepoint`, `residue_exposure`, `segment_compaction`, `contact_probability`, `operator_activation`, `frustration`, `state_basin_occupancy`, `interface_readiness`, `disorder_order_balance`, and `proteostasis_routing`.

## Failure Modes

Coordinate leakage, internal-runtime leakage, holdout-before-seal leakage, wrong grammar overpromotion, generic annotation-only promotion, one generic output for all regimes, and perturbation direction failure.
