# V67 E64 Mixed Fast Discovery 200

V67 is the first E64 mixed fast discovery shard. It keeps the loop at `200`
targets and combines fresh discovery with direct replay and regression guards.

The runner is:

```bash
python3 scripts/run_v67_e64_mixed_fast_discovery_200_v0.py
```

Target composition:

- `100` new RCSB nonredundant discovery targets
- `70` V66 failed-accepted targets replayed under E64
- `30` E64-compatible sentinels from V62/V64/V65/V66

Core outputs:

- `data/protein_esperanto_engine/E64/e64_v66_tm_complex_lineage_revision_certificate.json`
- `data/protein_esperanto_engine/V67/v67_mixed_fast_discovery_target_manifest.json`
- `data/protein_esperanto_engine/V67/v67_mixed_fast_discovery_scoring_report.json`
- `data/protein_esperanto_engine/V67/v67_mixed_fast_discovery_failure_report.json`
- `data/protein_esperanto_engine/V67/v67_old_v66_failure_repair_report.json`
- `data/protein_esperanto_engine/V67/v67_mixed_fast_discovery_certificate.json`
- `first_contact_clean_pharmacotopology_layer_run/V67_RCSB_NONREDUNDANT_200_DISCOVERY_E64/v67_e64_mixed_fast_discovery_200_certificate.json`
- `first_contact_clean_pharmacotopology_layer_run/V67_RCSB_NONREDUNDANT_200_DISCOVERY_E64/V67_E64_MIXED_FAST_DISCOVERY_200_REPORT.md`

V67/E64 reports:

- status: `V67_E64_MIXED_FAST_DISCOVERY_FAILURES_MINED_REVIEW_REQUIRED`
- targets: `200`
- supported: `93`
- failed accepted: `106`
- abstain: `1`
- accepted accuracy: `0.46733668341708545`
- raw accuracy: `0.465`
- coverage: `0.995`
- controls: `18/18`
- sentinel regressions: `0`

Four answers:

1. E64 repaired `1` of the `70` V66 failed-accepted targets; `69` remain.
2. The dominant failure mode is `assembly_required_core_vs_membrane_topology`
   with `48` cases.
3. E64 created `0` sentinel regressions in the mixed shard.
4. The engine is still `over_accepting_relative_to_abstention`.

Failure mode table:

| failure_mode | count |
| --- | ---: |
| `assembly_required_core_vs_membrane_topology` | `48` |
| `cofactor_locked_basin_vs_membrane_topology` | `33` |
| `short_interface_context_missing` | `12` |
| `disorder_or_low_complexity_context` | `12` |
| `assembly_required_core` | `1` |
| `over_abstain_membrane_multidomain_folding_proteostasis` | `1` |

The next missing Esperanto word is the assembly/topology separation: the engine
accepts assembly-required core grammar where reality shows membrane topology.
That points to `E65_ASSEMBLY_REQUIRED_FOLDING_GRAMMAR` with a focused
oligomer/assembly panel before any larger batch.
