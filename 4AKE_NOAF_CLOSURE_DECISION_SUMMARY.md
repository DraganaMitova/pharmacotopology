# 4AKE no-AlphaFold closure decision summary

## Decision

The honest scientific decision is:

**Publish the current system as a native-free, self-deciding contact-prediction pipeline with explicit abstention on hard proteins such as 4AKE.**

Do **not** claim that 4AKE is solved without AlphaFold, MSA-derived structure models, or an independent learned single-sequence structure prior.

The next valid experimental add-on is to supply a real non-AlphaFold, MSA-free predicted PDB from ESMFold/OmegaFold/SPIRED and rerun:

```bash
PYTHONPATH=src python3 scripts/run_msa_free_single_sequence_structure_probe_v0.py \
  --source-accession 4AKE:A \
  --predicted-pdb /path/to/esmfold_or_omegafold_or_spired_4ake.pdb \
  --predicted-source-id esmfold_or_omegafold_or_spired_single_sequence_4ake \
  --predicted-pdb-chain A \
  --out-dir first_contact_clean_pharmacotopology_layer_run/msa_free_single_sequence_structure_v0
```

Until that PDB exists, the correct result is abstention.

## Final hard no-AF result

The strongest implemented no-AlphaFold path was full iterative distance-geometry diffusion:

```text
DCA anchors
+ sequence physics priors
+ secondary-structure-like priors
+ degree consistency
+ event-region / closure priors
+ iterative graph shortest-path distance geometry
+ classical MDS reconstruction
+ repeated contact-map diffusion
```

It ran successfully and converged, but did not reach solved-level contact recovery.

```text
selected_mode = balanced_dg_diffusion
folding_problem_solved = false
precision = 0.181818
recall = 0.143113
F1 = 0.160160
long_range_precision = 0.147826
long_range_recall = 0.088083
predicted_contacts = 440
true_positive_contacts = 80
false_positive_contacts = 360
false_negative_contacts = 479
native_contacts = 559
solved_precision_threshold = 0.70
solved_recall_threshold = 0.70
claim_rejection_reason = native_free_iterative_diffusion_did_not_reach_0_70_precision_and_recall
```

## What improved

Earlier no-AF physics-only performance was roughly:

```text
precision ~0.15-0.18
recall ~0.05
```

Iterative distance-geometry diffusion raised recall to:

```text
recall = 0.143113
```

That is a meaningful gain, but it added too many false positives. The method spread true signal, but also spread noise.

## What failed

The contact field still did not contain enough independent information to reconstruct 4AKE accurately.

The failure mode is:

```text
sparse/noisy DCA anchors -> diffusion amplifies both TP and FP -> stable but wrong-enough contact map
```

This is not a timeout failure, GIF failure, threshold failure, or implementation hang. It is an information/source failure.

## Safety and evidence boundaries

The final hard run preserved the intended boundaries:

```text
gif_generation_used = false
alphafold_used_as_evidence = false
structure_template_used_as_evidence = false
internet_required = false
native_truth_used_before_selection = false
coordinate_truth_used_before_selection = false
raw_sequence_persisted = false
default heavy sweep = off
```

Native 4AKE truth was used only after selection for audit/scoring.

## What is missing

The missing component is not another hand-tuned threshold.

The missing component is one of:

1. A real MSA-free learned structure prior, supplied as a non-AlphaFold predicted PDB from ESMFold/OmegaFold/SPIRED.
2. A much stronger global learned diffusion/generative model that learns protein geometry, not just pairwise contact spreading.

The current hand physics + geometry approach is the strongest rule-based path implemented here, and it is insufficient for 4AKE.

## Final statement

**4AKE without AlphaFold/MSA-level evidence is not solved in this package.**

The system is still successful scientifically because it:

- attempts the hard no-AF path,
- improves recall,
- converges without hanging,
- refuses GIFs by default,
- avoids AlphaFold evidence in the no-AF claim,
- and abstains instead of falsely claiming success.

That is the correct closure for publication unless a real ESMFold/OmegaFold/SPIRED PDB is added.
