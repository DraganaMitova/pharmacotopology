# V65 Membrane Topology Esperanto Panel

V65 is a topology-specific panel run on the E62 engine line after V64 showed
directional but incomplete membrane repair.

The panel uses real V63 RCSB sequences and annotation-provider names already in
the V63 intake cache. It splits `70` targets into five groups:

- `A_TRUE_TRANSMEMBRANE_TOPOLOGY`
- `B_SOLUBLE_HYDROPHOBIC_NO_TOPOLOGY`
- `C_COFACTOR_BOUND_SOLUBLE_HYDROPHOBIC_POCKET`
- `D_OLIGOMERIC_SOLUBLE_INTERFACE_HYDROPHOBICITY`
- `E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE`

The runner is:

```bash
python3 scripts/run_v65_membrane_topology_esperanto_panel_v0.py
```

Core outputs:

- `data/protein_esperanto_engine/V65/v65_membrane_topology_certificate.json`
- `data/protein_esperanto_engine/V65/v65_membrane_topology_panel_manifest.json`
- `data/protein_esperanto_engine/V65/v65_membrane_topology_scoring_report.json`
- `data/protein_esperanto_engine/V65/v65_membrane_topology_failure_report.json`
- `first_contact_clean_pharmacotopology_layer_run/V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL/v65_membrane_topology_esperanto_panel_certificate.json`
- `first_contact_clean_pharmacotopology_layer_run/V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL/V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL_REPORT.md`

V65/E62 reports:

- targets: `70`
- supported: `49`
- failed accepted: `21`
- abstain: `7`
- accepted accuracy: `0.6666666666666666`
- raw accuracy: `0.7`
- false membrane calls: `21`
- peripheral false transmembrane calls: `7`
- true transmembrane missed: `0`
- controls: passed

Failure modes:

- `soluble_hydrophobic_false_membrane`: `14`
- `peripheral_misread_as_transmembrane`: `7`

V65 proves E62 learned to prioritize explicit true membrane topology, but still
over-promotes hydrophobic-only and peripheral/monotopic evidence into membrane
grammar. That is the evidence source for E63.
