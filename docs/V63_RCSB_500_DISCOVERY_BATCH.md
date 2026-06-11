# V63 RCSB 500 Discovery Batch

V63 expands the saturation campaign from the V61/V62 100-target repair lane to
500 automatically selected RCSB/PDB protein polymer entities.

This is a discovery batch:

- `batch_mode = discovery`
- `engine_version_used = E61`
- `target_selection_manual = false`
- RCSB experimental protein polymer entities only
- 30% sequence-identity cluster representatives
- `N = 500`
- length filter `40-800` residues
- coordinates, native contacts, AlphaFold-style models, post-seal validation
  annotations, and runtime artifacts blocked before sealing
- engine not modified inside V63

The runner is:

```bash
python3 scripts/run_v63_rcsb_500_discovery_batch_v0.py --refresh-intake
```

After the intake cache exists, deterministic regeneration uses:

```bash
python3 scripts/run_v63_rcsb_500_discovery_batch_v0.py
```

Core outputs:

- `data/protein_esperanto_engine/V63/v63_rcsb_500_certificate.json`
- `data/protein_esperanto_engine/V63/v63_rcsb_500_scoring_report.json`
- `data/protein_esperanto_engine/V63/v63_rcsb_500_failure_report.json`
- `first_contact_clean_pharmacotopology_layer_run/V63_RCSB_500_DISCOVERY_BATCH/v63_rcsb_500_discovery_batch_certificate.json`
- `first_contact_clean_pharmacotopology_layer_run/V63_RCSB_500_DISCOVERY_BATCH/V63_RCSB_500_DISCOVERY_BATCH_REPORT.md`

V63/E61 reports:

- accepted: `500`
- supported: `238`
- failed accepted: `262`
- abstain: `0`
- accepted accuracy: `0.476`
- raw accuracy: `0.476`
- coverage: `1.0`
- controls: `19/19`

Top mined failure modes:

- `membrane_misread`: `216`
- `wrong_regime`: `29`
- `disorder_misread`: `15`
- `oligomer_state_misread`: `2`

The dominant V63 finding is that E61 over-prioritized incidental
ligand/cofactor and assembly/copy context over membrane topology context in
broad RCSB space. That failure class is promoted to the E62 engine revision.
