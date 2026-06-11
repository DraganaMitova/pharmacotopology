# V52 Coarse Operator Folding Simulator MVP

V52 is the first running engine. It simulates residue and segment beads, not atoms.

## Scope

The MVP builds a sequence field, gates evidence, selects a mechanism grammar, creates an operator field, runs a coarse trajectory, emits predicted interaction probabilities and state-basin occupancy, and seals the result before validation.

## Explicit Non-Goals

No full atomistic MD, no AlphaFold-style coordinate prediction, no PDB/AF coordinate input before seal, and no claim that `folding_problem_solved` is true.

## Acceptance Controls

- Random sequence does not create target-specific mechanism.
- Shuffled sequence weakens operator coherence.
- Wrong evidence does not validate the target.
- Forced wrong grammar fails or abstains.
- Coordinate leakage blocks prediction.
- Wild type and mutant differ in the expected direction.

## Output Artifacts

For each target, V52 emits input evidence manifest, selected grammar, operator field, sequence-field map, trajectory summary, contact/interaction probability map, state-basin occupancy, perturbation table, falsifiers, and sealed packet hash.
