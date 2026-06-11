# V62 E61 Repair and Saturation Rerun

V62 reruns the exact `V61_RCSB_NONREDUNDANT_100_BATCH` target list with the
E61 engine line.

This is a repair batch, not a fresh broad-discovery claim:

- `batch_mode = repair`
- `engine_lineage = E60 -> E61`
- `same_targets_as_v61 = true`
- `known_failure_repair_probe = true`
- `engine_modified_during_batch = false`
- README remains manual-owned and untouched

The V62 runner is:

```bash
python3 scripts/run_v62_e61_repair_and_saturation_rerun_v0.py
```

Core outputs:

- `data/protein_esperanto_engine/V62/v62_e61_repair_certificate.json`
- `data/protein_esperanto_engine/V62/v62_e60_vs_e61_comparison.json`
- `data/protein_esperanto_engine/V62/v62_e61_repair_scoring_report.json`
- `first_contact_clean_pharmacotopology_layer_run/V62_E61_REPAIR_AND_SATURATION_RERUN/v62_e61_repair_and_saturation_rerun_certificate.json`
- `first_contact_clean_pharmacotopology_layer_run/V62_E61_REPAIR_AND_SATURATION_RERUN/V62_E61_REPAIR_AND_SATURATION_RERUN_REPORT.md`

V62 uses the committed V61/E60 result as the baseline:

- accepted: `81`
- supported: `8`
- failed accepted: `73`
- accepted accuracy: `0.09876543209876543`

V62/E61 reports:

- accepted: `96`
- supported: `90`
- failed accepted: `6`
- abstain: `4`
- accepted accuracy: `0.9375`
- raw accuracy: `0.9`
- coverage: `0.96`

Net movement:

- supported: `+82`
- failed accepted: `-67`
- abstain: `-15`
- accepted accuracy: `+0.8387345679012346`
- coverage: `+0.1499999999999999`

Failure-mode movement:

- `cofactor_ligand_missing`: `43 -> 0`
- `oligomer_state_misread`: `13 -> 0`
- `membrane_misread`: `11 -> 0`
- `weak_sequence_signal`: `19 -> 4`
- `disorder_misread`: `4 -> 4`
- `wrong_regime`: `2 -> 2`

Remaining failures are preserved as the next failure-mining material. V62 does
not mutate the engine inside the batch, and it does not claim that E61 has
saturated broad protein space.
