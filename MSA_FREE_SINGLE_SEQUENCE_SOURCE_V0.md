# MSA-free single-sequence structure source v0

This pass adds the missing source boundary for 4AKE/no-AlphaFold work:

- The existing DCA/self-consistency path can stabilize weak contact evidence, but it cannot create an independent structure signal when the MSA/contact signal is weak.
- Hand-coded sequence chemistry/secondary-structure/degree priors are useful as non-claim votes, but they are not an independent structure source.
- The new missing slot is a bounded MSA-free single-sequence structure source: ESMFold/OmegaFold/SPIRED/HelixFold-Single-style output, or any custom local predictor that writes a PDB.

The adapter is deliberately safe by default:

```bash
PYTHONPATH=src python3 scripts/run_msa_free_single_sequence_structure_probe_v0.py \
  --source-accession 4AKE:A \
  --out-dir first_contact_clean_pharmacotopology_layer_run/msa_free_single_sequence_structure_v0
```

With no predictor output it abstains:

```text
script_safety_rejection = missing_single_sequence_predictor_output
benchmark_claim_allowed_by_ensemble = false
folding_problem_solved = false
```

To use a precomputed non-AlphaFold single-sequence model:

```bash
PYTHONPATH=src python3 scripts/run_msa_free_single_sequence_structure_probe_v0.py \
  --source-accession 4AKE:A \
  --predicted-pdb /path/to/omegafold_or_esmfold_4ake.pdb \
  --predicted-source-id omegafold_single_sequence_4ake \
  --predicted-pdb-chain A \
  --out-dir first_contact_clean_pharmacotopology_layer_run/msa_free_single_sequence_structure_v0
```

To run a local predictor command, pass it explicitly; it is never automatic:

```bash
PYTHONPATH=src python3 scripts/run_msa_free_single_sequence_structure_probe_v0.py \
  --source-accession 4AKE:A \
  --prediction-command "your_predictor --fasta {fasta} --out {output_pdb}" \
  --predicted-source-id custom_single_sequence_structure_model \
  --timeout-seconds 180 \
  --out-dir first_contact_clean_pharmacotopology_layer_run/msa_free_single_sequence_structure_v0
```

Important boundary rules:

- AlphaFold-like source IDs are rejected by default in this MSA-free probe.
- Predictor execution is opt-in only.
- Every predictor subprocess has a timeout.
- Query FASTA is temporary by default and is not persisted.
- GIF generation remains off by default; this path does not render GIFs.
- Native coordinates are only attached after selection for audit metrics.
