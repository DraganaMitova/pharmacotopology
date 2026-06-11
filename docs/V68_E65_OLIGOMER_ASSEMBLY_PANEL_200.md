# V68 E65 Oligomer Assembly Panel 200

V68 tests `E65_ASSEMBLY_REQUIRED_FOLDING_GRAMMAR`, the missing Esperanto word
exposed by V67:

```text
assembly_required_core_vs_topology_provider
```

The runner is:

```bash
python3 scripts/run_v68_e65_oligomer_assembly_panel_200_v0.py
```

Target composition:

- `48` V67 assembly-required failure replays
- `40` biological oligomer / assembly-required positives
- `30` coiled-coil or helix-bundle assembly candidates
- `30` true transmembrane oligomer/channel/pore sentinels
- `25` soluble monomeric hydrophobic-core sentinels
- `15` cofactor/ligand soluble sentinels
- `12` generic-complex-but-not-assembly-required controls

Core outputs:

- `data/protein_esperanto_engine/E65/e65_assembly_required_folding_grammar_certificate.json`
- `data/protein_esperanto_engine/V68/v68_oligomer_assembly_panel_target_manifest.json`
- `data/protein_esperanto_engine/V68/v68_oligomer_assembly_panel_scoring_report.json`
- `data/protein_esperanto_engine/V68/v68_oligomer_assembly_panel_failure_report.json`
- `data/protein_esperanto_engine/V68/v68_oligomer_assembly_panel_certificate.json`
- `first_contact_clean_pharmacotopology_layer_run/V68_OLIGOMER_ASSEMBLY_PANEL_200/v68_e65_oligomer_assembly_panel_200_certificate.json`
- `first_contact_clean_pharmacotopology_layer_run/V68_OLIGOMER_ASSEMBLY_PANEL_200/V68_E65_OLIGOMER_ASSEMBLY_PANEL_200_REPORT.md`

V68/E65 reports:

- status: `V68_E65_OLIGOMER_ASSEMBLY_PANEL_PASSED_REVIEW_REQUIRED`
- targets: `200`
- supported: `200`
- failed accepted: `0`
- abstain: `12`
- accepted accuracy: `1.0`
- raw accuracy: `1.0`
- coverage: `0.94`
- controls: `17/17`
- sentinel regressions: `0`

Required V68 metrics:

- assembly_required_correct: `118`
- true_TM_preserved: `30`
- soluble_hydrophobic_not_called_assembly: `25`
- cofactor_not_called_assembly: `15`
- generic_complex_not_called_assembly: `12`
- coiled_coil_register_detected: `30`
- domain_swap_candidates_detected: `3`
- assembly_ambiguous_abstained: `0`

V67 dominant failure replay:

- replay count: `48`
- repaired by E65: `48`
- remaining: `0`

The generic-complex controls are the important calibration check: all `12`
abstained cleanly instead of becoming false assembly calls.

Next batch:

```text
V69_RCSB_NONREDUNDANT_200_DISCOVERY_E65
```

V68 is a specialist grammar panel. It does not make a broad solved-folding
claim; the next move is a fresh 200-target RCSB discovery shard on E65.
