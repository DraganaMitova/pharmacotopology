# V58 Real Sequence Time-Blind Folding Replication Gate

V58 freezes the Protein Esperanto engine and exposes it to automatically selected real protein sequences from a public recent-release source.

## Source Of Truth

The first V58 implementation uses RCSB PDB recent-release protein polymer entities as a CAMEO-style seal-then-open benchmark. The raw intake is cached in `data/protein_esperanto_engine/V58/intake/` so tests are deterministic after the network fetch.

This is not a live CASP/CAMEO submission. It is a local replication gate that imitates the discipline:

```text
select real sequences automatically
-> hide coordinate/native-contact answers from prediction
-> seal predictions
-> open post-seal validation labels
-> report failures
```

## Pre-Seal Inputs

Allowed:

- amino-acid sequence
- length
- organism/taxonomy
- basic sequence-derived marks
- pure non-coordinate metadata

Blocked:

- PDB/mmCIF coordinates
- native contacts
- AlphaFold/ESMFold/RoseTTAFold models
- structure-derived domains
- PDB-derived interface contacts
- post-seal validation papers or reports
- internal runtime artifacts as biological evidence

## Runs

V58 emits both:

- `sequence_only`: raw amino-acid sequence only
- `sequence_plus_annotation`: sequence plus pure metadata and sequence-derived marks

## Proof Levels

- Level 1: regime selection
- Level 2: region localization proxy
- Level 3: topology or observable behavior
- Level 4: process replication, explicitly not claimed until V59

## Required Certificate Fields

The certificate must include:

- `folding_problem_solved = false`
- `atomistic_md_executed = false`
- `engine_modified_after_target_selection = false`
- `coordinate_truth_used_before_seal = false`
- `alphafold_used_before_seal = false`
- `target_selection_manual = false`
- `failure_cases_reported = true`

## Claim Boundary

Allowed if passed:

```text
A frozen Protein Esperanto engine generalized to automatically selected real protein sequences and predicted folding regime, operator trajectory, and post-seal structural/experimental observables under leakage-controlled validation.
```

Forbidden:

```text
Universal protein folding is solved.
Coordinates were predicted de novo.
Atomistic folding was solved.
AlphaFold was used before sealing.
Process replication was proven from PDB structures alone.
External review is unnecessary.
```
