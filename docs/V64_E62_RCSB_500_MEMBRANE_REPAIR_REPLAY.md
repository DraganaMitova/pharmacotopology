# V64 E62 RCSB 500 Membrane Repair Replay

V64 reruns the exact `V63_RCSB_500_DISCOVERY_BATCH` target list with the E62
engine line.

This is a paired repair replay:

- `batch_mode = paired_repair_replay`
- `engine_version_used = E62`
- `baseline_engine_version = E61`
- `same_targets_as_v63 = true`
- target count: `500`
- target order and 30% sequence-cluster representatives are inherited from V63
- coordinates, native contacts, ligand geometry, structure models, validation
  labels, and runtime artifacts are blocked before sealing
- engine not modified inside V64

The runner is:

```bash
python3 scripts/run_v64_e62_rcsb_500_membrane_repair_replay_v0.py
```

Core outputs:

- `data/protein_esperanto_engine/V64/v64_e62_rcsb_500_membrane_repair_certificate.json`
- `data/protein_esperanto_engine/V64/v64_e61_vs_e62_comparison.json`
- `data/protein_esperanto_engine/V64/v64_e62_rcsb_500_membrane_repair_scoring_report.json`
- `data/protein_esperanto_engine/V64/v64_e62_rcsb_500_membrane_repair_failure_report.json`
- `first_contact_clean_pharmacotopology_layer_run/V64_E62_RCSB_500_MEMBRANE_REPAIR_REPLAY/v64_e62_rcsb_500_membrane_repair_replay_certificate.json`
- `first_contact_clean_pharmacotopology_layer_run/V64_E62_RCSB_500_MEMBRANE_REPAIR_REPLAY/V64_E62_RCSB_500_MEMBRANE_REPAIR_REPLAY_REPORT.md`

V63/E61 baseline:

- supported: `238`
- failed accepted: `262`
- membrane_misread: `216`
- abstain: `0`
- accepted accuracy: `0.476`

V64/E62 replay:

- supported: `328`
- failed accepted: `172`
- membrane_misread: `123`
- abstain: `0`
- accepted accuracy: `0.656`
- raw accuracy: `0.656`
- coverage: `1.0`
- controls: `20/20`

Net movement:

- supported: `+90`
- failed accepted: `-90`
- membrane_misread: `-93`
- accepted accuracy: `+0.18000000000000005`
- coverage: `0.0`

Repair map:

- stable supported: `235`
- repaired: `93`
- persistent failure: `169`
- new failure: `3`

The V64 result is directional but incomplete. E62 repaired part of the V63
membrane wall without a large soluble false-membrane explosion, but abstention
did not return and `123` membrane failures remained. That forces the dedicated
V65 membrane topology panel before any larger V66 expansion.
