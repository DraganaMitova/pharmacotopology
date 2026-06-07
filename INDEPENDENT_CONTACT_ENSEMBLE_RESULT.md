# Independent Contact Ensemble Probe V0

## Decision

This batch implements the next honest branch after Phase-3:

- not another frontier threshold
- not a default selector replacement
- an independent contact evidence ingestion layer
- an ensemble verifier that refuses to make a benchmark claim unless a non-native independent structure/contact source is present

## Why this exists

Phase-3 proved that 4AKE is not blocked by candidate generation ceiling:

- broad candidate pool ceiling: 1.000
- aggressive frontier can reopen the map
- but collapse precision dies when only DCA + sequence closure are available

Therefore, the next valid branch is independent evidence. The system now has a native-free interface for AlphaFold/RoseTTAFold/contact-map style sources.

## Added

- `src/pharmacotopology/folding_independent_contact_evidence.py`
- `scripts/run_independent_contact_ensemble_probe_v0.py`
- `tests/test_independent_contact_ensemble_v0.py`

The probe combines evidence families:

1. `candidate_region` — current sequence-closure/candidate frontier support
2. `external_coupling` — exact DCA/coupling pair support
3. `independent_structure` — AlphaFold/RoseTTAFold/contact-map/PDB-like predicted structure contacts

Default selection requires:

- at least 2 source-family votes
- candidate-region support
- independent-structure support

That means DCA + sequence closure alone no longer pretends to rescue 4AKE.

## Real independent source missing in this sandbox

The current archive does not include an AlphaFold/RoseTTAFold predicted structure or independent contact JSON for 4AKE. The probe therefore produces two honest artifacts:

### No independent source

Command:

```bash
PYTHONPATH=src python3 scripts/run_independent_contact_ensemble_probe_v0.py \
  --report-output first_contact_clean_pharmacotopology_layer_run/independent_contact_ensemble_4ake_no_independent_source_v0.json
```

Result:

```text
benchmark_claim_allowed = false
claim_rejection_reason = missing_independent_structure_source
final_pair_count = 0
long_range_recall = 0.0
```

This is intentional. No independent source means no rescue claim.

### Native-coordinate positive control

Command:

```bash
PYTHONPATH=src python3 scripts/run_independent_contact_ensemble_probe_v0.py \
  --native-coordinate-positive-control \
  --report-output first_contact_clean_pharmacotopology_layer_run/independent_contact_ensemble_4ake_positive_control_v0.json
```

Result:

```text
benchmark_claim_allowed = false
claim_rejection_reason = independent_source_is_native_coordinate_leakage_positive_control
final_pair_count = 259
final_long_range_pair_count = 193
long_range_precision = 1.0
long_range_recall = 1.0
```

This proves the new ensemble gate can filter the broad 4AKE frontier if a strong independent source exists, but it does **not** count as a scientific win because it uses benchmark coordinates.

## How to run the real next test

Download or generate a real predicted structure/contact source for the same 4AKE sequence, for example an AlphaFold DB / RoseTTAFold / ColabFold PDB model. Then run:

```bash
PYTHONPATH=src python3 scripts/run_independent_contact_ensemble_probe_v0.py \
  --predicted-pdb data/independent_contact_sources/AF-P69441-F1-model_v4.pdb \
  --predicted-source-id alphafold_db_AF-P69441-F1 \
  --report-output first_contact_clean_pharmacotopology_layer_run/independent_contact_ensemble_4ake_alphafold_v0.json \
  --decisions-output first_contact_clean_pharmacotopology_layer_run/independent_contact_ensemble_4ake_alphafold_decisions_v0.csv \
  --selected-output first_contact_clean_pharmacotopology_layer_run/independent_contact_ensemble_4ake_alphafold_selected_contacts_v0.csv \
  --evidence-output first_contact_clean_pharmacotopology_layer_run/independent_contact_ensemble_4ake_alphafold_evidence_v0.json
```

Benchmark claim is allowed only if:

```text
benchmark_claim_allowed = true
coordinate_truth_used_before_selection = false
native_truth_used_before_selection = false
raw_sequence_exposed = false
```

## Current conclusion

The code path for independent evidence is now present and guarded. With no independent source it refuses to claim anything. With a native-coordinate positive control it shows the expected rescue but marks it as leakage.

The next real answer depends on feeding the probe an actual independent AlphaFold/RoseTTAFold/contact-map file.
