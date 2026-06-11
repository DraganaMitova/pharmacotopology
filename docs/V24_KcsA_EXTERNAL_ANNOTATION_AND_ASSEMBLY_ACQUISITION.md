# V24 KcsA External Annotation and Assembly Acquisition

V24 is a zero-MD acquisition sprint for KcsA. It locks external membrane/pore/assembly-context source references and preserves V19 pressure evidence while keeping all high-risk claims disabled.

## Locked behavior

- KcsA is treated as a `membrane_pore_oligomer_object`.
- Pore/filter and potassium-ion diagnostic shell evidence are preserved.
- Transmembrane helix scaffold context is preserved.
- Chain/interface context is preserved.
- Soluble compact-core misclassification remains forbidden.
- Tetramer and whole-channel fold claims remain forbidden.
- No membrane MD is executed.
- External couplings are scanned but not synthesized and not required for V24.

## Claim boundary

`positive_pressure_evidence_found=true` is allowed only when evidence is contextual and guarded. `positive_folding_evidence_found` remains false and `claim_allowed` remains false.
