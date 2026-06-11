# V66 E63 Fast Membrane Repair Replay 200

V66 starts the adaptive 200-target learning loop. It does not rerun the whole
V63/V64 500-target set.

The 200-target repair set is:

- `70` targets from the V65 membrane topology panel
- `80` selected V64/V63 membrane_misread failures
- `3` V64 E62 oligomer_state_misread regressions
- `47` sentinels that E61/E62 had already handled correctly

The runner is:

```bash
python3 scripts/run_v66_e63_fast_membrane_repair_replay_200_v0.py
```

Core outputs:

- `data/protein_esperanto_engine/V66/v66_e63_fast_membrane_repair_certificate.json`
- `data/protein_esperanto_engine/V66/v66_e62_vs_e63_comparison.json`
- `data/protein_esperanto_engine/V66/v66_e63_fast_membrane_repair_scoring_report.json`
- `data/protein_esperanto_engine/V66/v66_e63_fast_membrane_repair_failure_report.json`
- `first_contact_clean_pharmacotopology_layer_run/V66_E63_FAST_MEMBRANE_REPAIR_REPLAY_200/v66_e63_fast_membrane_repair_replay_200_certificate.json`
- `first_contact_clean_pharmacotopology_layer_run/V66_E63_FAST_MEMBRANE_REPAIR_REPLAY_200/V66_E63_FAST_MEMBRANE_REPAIR_REPLAY_200_REPORT.md`

Same-200 E62 baseline:

- supported: `96`
- failed accepted: `104`
- abstain: `7`
- accepted accuracy: `0.46113989637305697`
- raw accuracy: `0.48`
- coverage: `0.965`

V66/E63 replay:

- supported: `130`
- failed accepted: `70`
- abstain: `28`
- accepted accuracy: `0.5930232558139535`
- raw accuracy: `0.65`
- coverage: `0.86`
- controls: `15/15`

Repair movement:

- supported: `+34`
- failed accepted: `-34`
- abstain: `+21`
- accepted accuracy: `+0.13188335944089657`
- coverage: `-0.10499999999999998`

False membrane repair:

- false membrane calls: `24 -> 0`
- peripheral-as-transmembrane: `7 -> 0`
- true transmembrane misses under E63: `0`
- oligomer regressions remaining under E63: `0`
- sentinel regressions under E63: `0`

Remaining failure mode:

- `membrane_misread`: `70`

V66 passes the fast membrane repair replay. It also shows the next remaining
positive grammar gap: true membrane targets without enough topology evidence
still fail as membrane_misread. The next batch is the adaptive discovery shard
`V67_RCSB_NONREDUNDANT_200_DISCOVERY_E63`.
