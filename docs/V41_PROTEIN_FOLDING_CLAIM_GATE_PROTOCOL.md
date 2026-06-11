# V41 Protein Folding Claim Gate Protocol

## Purpose

V41 is a claim boundary audit. It decides the maximum honest scientific claim allowed by the locked V32-V40 evidence.

It is not another positive-result generator. It does not run MD, does not use coordinates, and does not mark protein folding solved.

## Claim Levels

- `C0_NO_CLAIM`: No scientific claim beyond code execution.
- `C1_CLEAN_EVIDENCE_PIPELINE`: The system cleanly imports external evidence and blocks leakage.
- `C2_MECHANISM_GRAMMAR_CLASSIFICATION`: The system classifies folding-problem grammar from non-coordinate evidence.
- `C3_FALSIFIABLE_OPERATOR_MECHANISM`: The system generates falsifiable, perturbation-sensitive operator-level mechanism predictions supported by independent non-coordinate holdouts.
- `C4_PROSPECTIVE_MECHANISM_GENERALIZATION`: The system generalizes prospectively to unseen targets under masking, decoy pressure, source separation, and independent holdout validation.
- `C5_PROTEIN_FOLDING_SOLVED`: The system can predict de novo 3D folds / functional conformational ensembles across broad protein classes with independent blind validation comparable to structure-prediction benchmarks.

## Current Boundary

The current gate is expected to allow at most `C3_FALSIFIABLE_OPERATOR_MECHANISM`.

`C5_PROTEIN_FOLDING_SOLVED` is blocked.

## Required Inputs

V41 audits certificates/reports from:

- V32
- V33
- V34
- V35
- V36
- V37
- V38
- V39
- V40

Runtime reports may be used only as audit evidence about pipeline behavior. They are not biological evidence.

## Required Outputs

V41 writes:

- `data/claim_gate/V41/claim_ladder.json`
- `data/claim_gate/V41/evidence_to_claim_matrix.csv`
- `data/claim_gate/V41/allowed_claim_text.md`
- `data/claim_gate/V41/forbidden_claim_text.md`
- `data/claim_gate/V41/scientist_attack_surface.md`
- `first_contact_clean_pharmacotopology_layer_run/V41_PROTEIN_FOLDING_CLAIM_GATE/v41_protein_folding_claim_gate_certificate.json`
- `first_contact_clean_pharmacotopology_layer_run/V41_PROTEIN_FOLDING_CLAIM_GATE/V41_PROTEIN_FOLDING_CLAIM_GATE_REPORT.md`

## Required Controls

V41 blocks or fails if:

1. Any earlier certificate has `claim_allowed=true`.
2. `folding_problem_solved=true` appears anywhere.
3. Coordinate-derived source count for non-coordinate claims is greater than zero.
4. Internal runtime source count for biological evidence is greater than zero.
5. MD was used as claim evidence.
6. V35 abstain is treated as positive evidence.
7. V38 answer key was used for assignment.
8. C5 is allowed without de novo structure/ensemble validation.
9. Allowed claim text omits `not solved protein folding`.
10. Forbidden claims omit `we solved protein folding`.
11. Scientist attack surface is empty.
12. Max claim is below C2 despite V36-V40 passing.

## Allowed Status

The intended safe pass status is:

`V41_MECHANISM_CLAIM_ALLOWED_C5_BLOCKED`

That status means the bounded mechanism claim is allowed, while solved-protein-folding and C5 wording remain blocked.
