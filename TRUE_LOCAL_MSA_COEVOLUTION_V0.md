# TRUE_LOCAL_MSA_COEVOLUTION_V0

## Purpose

This batch closes the previous half-step. The earlier `local_coevolution_expansion_v0`
experiment used a proxy because the uploaded archive did not contain a raw MSA.
This version implements the real local co-evolution experiment and refuses to
silently fall back to the proxy.

## Raw MSA status in this ZIP

`raw_msa_presence_scan.json` found no raw MSA files in the uploaded repository:

```json
{
  "candidate_count": 0,
  "likely_raw_msa_count": 0,
  "raw_msa_for_4ake_found": false
}
```

Therefore the true 4AKE local-MSA experiment is implemented, tested with a
synthetic aligned MSA fixture, and packaged as a runnable pipeline, but the real
4AKE result still requires a raw MSA file.

## What was added

- `src/pharmacotopology/folding_true_local_msa_coevolution.py`
- `scripts/run_true_local_msa_coevolution_v0.py`
- `scripts/build_4ake_msa_with_jackhmmer_v0.sh`
- `scripts/check_raw_msa_presence_v0.py`
- `tests/test_true_local_msa_coevolution_v0.py`
- `first_contact_clean_pharmacotopology_layer_run/true_local_msa_coevolution_v0/ACTION_REQUIRED_RAW_MSA.md`

## Scientific boundary

The implemented algorithm is the user-requested one:

1. Take top safe DCA anchors.
2. For each anchor `(i, j)`, keep only MSA sequences where residues `i` and `j`
   match the 4AKE query residues.
3. Inside the anchor-conditioned sub-MSA, compute normalized mutual information
   for local window pairs `(i + di, j + dj)`.
4. Select high-scoring long-range local pairs with a degree cap.
5. Attach native contacts only after the prediction is frozen for audit.

No native contacts, coordinates, ESMFold, AlphaFold, or structure template are
used before selection.

## How to run locally

First build an MSA using a local sequence database:

```bash
PYTHONPATH=src bash scripts/build_4ake_msa_with_jackhmmer_v0.sh \
  --database /path/to/uniprot_or_uniref.fasta \
  --out-dir external_msa/4ake_jackhmmer \
  --cpu 4
```

Then run the true local co-evolution benchmark:

```bash
PYTHONPATH=src python3 scripts/run_true_local_msa_coevolution_v0.py \
  --source-accession 4AKE:A \
  --msa external_msa/4ake_jackhmmer/4ake_jackhmmer.afa \
  --top-anchor-count 50 \
  --window 5 \
  --threshold 0.5 \
  --min-filtered-sequences 16 \
  --out-dir first_contact_clean_pharmacotopology_layer_run/true_local_msa_coevolution_v0
```

If only Stockholm is available, pass the `.sto` file directly:

```bash
PYTHONPATH=src python3 scripts/run_true_local_msa_coevolution_v0.py \
  --source-accession 4AKE:A \
  --msa external_msa/4ake_jackhmmer/4ake_jackhmmer.sto
```

## Validation performed here

Focused only; no full suite:

```text
PYTHONPATH=src python3 -m pytest -q tests/test_true_local_msa_coevolution_v0.py
3 passed in 1.60s
```

Missing-real-MSA runner path:

```json
{
  "status": "action_required_raw_msa_missing",
  "proxy_fallback_used": false,
  "raw_msa_available_for_true_local_mi": false
}
```

Compileall:

```text
OK
```

Final CPU cleanup:

```json
{"matched_count": 0}
```

## Decision

This is not a solved result yet, because the real 4AKE raw MSA was not present
in the archive and the sandbox cannot fetch/build one from a sequence database.
But the missing piece is now executable, not conceptual. The next local run with
`4ake_jackhmmer.afa` will produce the actual precision/recall/long-range audit.
