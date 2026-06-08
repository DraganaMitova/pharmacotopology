# external_msa_free_predictors

This directory is the real 4AKE rescue harness.

## One command

```bash
./run_all_4ake_msa_free_tryhard.sh /path/to/project
```

It tries:

1. ESMFold/ESM Atlas API
2. Local ESMFold v1 if installed
3. OmegaFold if installed or cloneable
4. Any externally dropped PDB files under `tryhard_runs/`

Then it runs the project's no-AlphaFold probe on every usable PDB.

## Direct ESMFold API only

```bash
./esmfold_api.sh
```

or

```bash
./esmfold_api.sh 4ake.fasta 4ake_esmfold_api.pdb
```

## Drop-in mode

If you generate a PDB from ESMFold2, SPIRED, OmegaFold, Biohub Platform, Colab, or another single-sequence predictor:

```bash
mkdir -p tryhard_runs/manual
cp /path/to/your_4ake_prediction.pdb tryhard_runs/manual/
./run_all_4ake_msa_free_tryhard.sh /path/to/project
```

## Claim threshold

```text
precision >= 0.70 and recall >= 0.70
```
