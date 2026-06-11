# V65 Membrane Topology Esperanto Panel

Status: `V65_MEMBRANE_TOPOLOGY_PANEL_FAILURES_REVIEW_REQUIRED`
Targets total: `70`
Supported: `49`
Failed accepted: `21`
Abstain: `7`
Accepted accuracy: `0.6666666666666666`
Raw accuracy: `0.7`
Controls: `13/13`
False membrane calls: `21`
Peripheral false transmembrane calls: `7`
True transmembrane missed: `0`

## Groups
- `A_TRUE_TRANSMEMBRANE_TOPOLOGY`: `14` targets, labels `{'supported': 14}`
- `B_SOLUBLE_HYDROPHOBIC_NO_TOPOLOGY`: `14` targets, labels `{'contradicted': 14}`
- `C_COFACTOR_BOUND_SOLUBLE_HYDROPHOBIC_POCKET`: `14` targets, labels `{'supported': 14}`
- `D_OLIGOMERIC_SOLUBLE_INTERFACE_HYDROPHOBICITY`: `14` targets, labels `{'supported': 14}`
- `E_PERIPHERAL_MEMBRANE_ASSOCIATED_NOT_TRANSMEMBRANE`: `14` targets, labels `{'contradicted': 7, 'supported': 7}`

## Failure Modes
- `soluble_hydrophobic_false_membrane`: `14`
- `peripheral_misread_as_transmembrane`: `7`

## Boundary
V65 is a membrane-topology grammar panel. It distinguishes topology evidence from hydrophobicity, cofactor pockets, oligomeric interfaces, and peripheral membrane association; it does not make a broad folding claim.
