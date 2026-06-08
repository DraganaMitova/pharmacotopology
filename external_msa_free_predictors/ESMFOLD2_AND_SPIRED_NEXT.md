# 4AKE strongest single-sequence predictors to try next

This package now targets more than ESMFold v1.

## Current strongest no-MSA axis

1. **ESMFold2 / ESMC** — newly released Biohub model. This is the highest-aim route because it updates the old ESMFold v1 idea with ESMC representations and a diffusion-based structure predictor. Use this if you have Biohub Platform access, Modal, or a machine that can download the open weights.
2. **ESMFold v1 API / local** — easiest if the public ESM Atlas endpoint works.
3. **OmegaFold** — primary-sequence predictor, no MSA.
4. **SPIRED** — single-sequence predictor with strong speed/accuracy tradeoff.

## Expected truth

If any of these produces a PDB, run:

```bash
PYTHONPATH=src python3 scripts/run_msa_free_single_sequence_structure_probe_v0.py \
  --source-accession 4AKE:A \
  --predicted-pdb /path/to/predicted.pdb \
  --predicted-source-id esmfold2_single_sequence_4ake \
  --predicted-pdb-chain A \
  --include-sequence-physical-priors \
  --min-votes 1
```

The claim threshold remains:

```text
precision >= 0.70 and recall >= 0.70
```

No AlphaFold-like file is accepted unless explicitly overridden for debugging.
