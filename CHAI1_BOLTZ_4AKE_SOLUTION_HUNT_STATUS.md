# 4AKE solution hunt status: Chai-1/Boltz ensemble path

## Honest status

The latest attempt did **not** produce a solved 4AKE structure inside this sandbox.

What was completed:

- Added a bounded Chai-1/Boltz-style conformational-ensemble adapter.
- Added a 4AKE probe script that can consume already-generated Chai-1/Boltz PDB/mmCIF samples.
- Added a hard-timeout local predictor-command path.
- Added provenance rejection for AlphaFold/ColabFold-like inputs by default.
- Added no-MSA command rejection unless explicitly overridden.
- Validated that the adapter selects a conformational sample without using 4AKE native truth.
- Validated that the full test suite does not hang.

What was **not** completed:

- No actual Chai-1/Boltz neural folding run completed inside the sandbox.
- `nvidia-smi` is unavailable, so this environment does not expose a CUDA GPU.
- `chai_lab` and `boltz` are not installed in the runtime.
- ESMFold import exists through `esm`, but loading `esm.pretrained.esmfold_v1()` fails because `openfold` is missing.
- Therefore no real non-AlphaFold PDB/mmCIF was generated here.

## Why this path was added

4AKE is a conformational/dynamics-heavy adenylate kinase case. A single static contact map is likely the wrong target when the predictor source can generate multiple conformational samples. Chai-1 is useful for this because its documented command-line inference generates multiple samples by default and uses embeddings without MSAs/templates unless MSA/template flags are supplied. Boltz also has an explicit single-sequence mode (`msa: empty`), although its documentation notes this reduces accuracy.

## New files

```text
src/pharmacotopology/folding_chai1_conformational_ensemble.py
scripts/run_4ake_chai1_conformational_ensemble_probe_v0.py
tests/test_chai1_conformational_ensemble_v0.py
CHAI1_BOLTZ_4AKE_SOLUTION_HUNT_STATUS.md
```

## Main command when real Chai-1/Boltz output exists

```bash
PYTHONPATH=src python3 scripts/run_4ake_chai1_conformational_ensemble_probe_v0.py \
  --source-accession 4AKE:A \
  --predicted-structure-dir /path/to/chai_or_boltz_4ake_outputs \
  --predicted-source-id chai1_single_sequence_conformational_ensemble \
  --predicted-chain A \
  --out-dir first_contact_clean_pharmacotopology_layer_run/4ake_chai1_conformational_ensemble_v0
```

## Bounded command mode

```bash
PYTHONPATH=src python3 scripts/run_4ake_chai1_conformational_ensemble_probe_v0.py \
  --source-accession 4AKE:A \
  --prediction-command "chai-lab fold {fasta} {out_dir}" \
  --predicted-source-id chai1_single_sequence_conformational_ensemble \
  --timeout-seconds 600 \
  --out-dir first_contact_clean_pharmacotopology_layer_run/4ake_chai1_conformational_ensemble_v0
```

The FASTA is written only inside a temporary directory and deleted when the command exits or times out.

## Safety / no-hang behavior

```text
GIF generation: not used
AlphaFold-like input: rejected by default
MSA-like command: rejected by default
Native 4AKE truth before selection: false
Coordinate truth before selection: false
Raw sequence persisted: false
Prediction subprocess timeout: configurable
Full suite default: bounded and no visual tests unless opt-in
```

## Current local validation

```text
9 passed, 2 skipped in 15.73s
```

## Current no-input Chai/Boltz probe result

```text
parsed_conformer_count = 0
selected_contact_count = 0
precision = 0.0
recall = 0.0
long_range_recall = 0.0
folding_problem_solved_after_native_audit = false
claim_rejection_reason = missing_or_insufficient_chai1_boltz_conformational_ensemble_output
```

## Bottom line

The project is now ready to test the most plausible next solution class: **MSA-free neural conformational ensemble prediction**.

But this sandbox did not generate a real neural PDB/mmCIF, so 4AKE is still not solved here.
