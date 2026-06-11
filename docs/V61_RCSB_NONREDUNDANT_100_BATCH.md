# V61 RCSB Nonredundant 100 Batch

V61 starts the `V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN`.

The batch freezes the Protein Esperanto engine at `E60`, whose engine source
last commit is `d927781`, and runs an automatically selected RCSB/PDB panel:

- `target_selection_manual = false`
- RCSB experimental protein polymer entities only
- 30% sequence-identity cluster representatives
- `N = 100`
- length filter `40-800` residues
- coordinates, native contacts, AlphaFold-style models, post-seal validation
  annotations, and internal runtime artifacts blocked before sealing
- no engine modification inside V61
- README remains manual-owned and untouched

The V61 runner is:

```bash
python3 scripts/run_v61_rcsb_nonredundant_100_batch_v0.py
```

The refresh path rebuilds the cached intake from the RCSB Search API grouped
polymer-entity query and RCSB Data API metadata:

```bash
python3 scripts/run_v61_rcsb_nonredundant_100_batch_v0.py --refresh-intake
```

Core outputs:

- `data/protein_esperanto_engine/V61/v61_rcsb_nonredundant_100_certificate.json`
- `first_contact_clean_pharmacotopology_layer_run/V61_RCSB_NONREDUNDANT_100_BATCH/v61_rcsb_nonredundant_100_batch_certificate.json`
- `data/protein_esperanto_engine/V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN/ledgers/protein_universe_ledger_v0.json`
- `data/protein_esperanto_engine/V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN/ledgers/engine_version_ledger_v0.json`
- `data/protein_esperanto_engine/V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN/ledgers/failure_grammar_ledger_v0.json`
- `data/protein_esperanto_engine/V61_TO_V70_PROTEIN_ESPERANTO_SATURATION_CAMPAIGN/ledgers/claim_ledger_v0.json`

The success target is not forced raw perfection. V61 reports:

- `accepted_accuracy`
- `raw_accuracy`
- `coverage`
- `controls_passed`
- `failed_accepted_count`
- `failure_modes`
- `missing_esperanto_candidates`

Any failed accepted prediction becomes a failure-grammar row for a later engine
revision. V61 itself does not repair those failures.
