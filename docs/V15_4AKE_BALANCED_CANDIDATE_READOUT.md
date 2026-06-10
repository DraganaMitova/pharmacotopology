# V15 4AKE Balanced Candidate Readout

This layer is postprocess-only. It reads an existing 4AKE OpenMM run and evaluates the effective balanced candidate set under the V15 dynamic role grammar.

It does **not**:

- rerun MD,
- use native precision for selection,
- synthesize evidence from PDB-only or visual-only artifacts,
- use a fixed residue-count separation cutoff.

It asks whether 4AKE has DCA-lane-backed balanced candidates that persist in the trajectory tail across replicas under a dynamic frequency sweep. If a clean frequency band is found, the result can be bridged into the V15 unified grammar panel. Claims remain disabled.
