# 4AKE ESMFold Try-Hard Audit Summary

Date: 2026-06-08

## Verdict

4AKE was solved by the MSA-free learned structure path, specifically by the direct ESMFold single-sequence PDB.

This does **not** mean the handcrafted physics / DCA / iterative DG path alone solved 4AKE. That path had previously failed. The solved mode is:

```text
folding_solution_mode = direct_msa_free_single_sequence_structure
```

## Input / source

The Mac-safe try-hard run produced a real PDB from the ESMFold API:

```text
status = success
http_status = 200
pdb_exists = true
pdb_atom_lines = 1655
sequence_length = 214
timeout_seconds = 600.0
```

## Direct structure audit

The direct ESMFold contact map passed the solved thresholds:

```text
native_contact_precision = 0.880068
native_contact_recall = 0.932021
contact_map_f1 = 0.9053
predicted_contact_count = 592
true_positive_contacts = 521
false_positive_contacts = 71
false_negative_contacts = 38
long_range_contact_precision = 0.76652
long_range_contact_recall = 0.901554
solved_precision_threshold = 0.70
solved_recall_threshold = 0.70
claim_allowed = true
```

## Ensemble / collapse audit

The contact ensemble improved precision and long-range recall, but did not fully solve all-contact recall:

```text
ensemble_contact_collapse_solved = false
contact_precision = 0.797945
contact_recall = 0.416816
final_pair_count = 292
true_positive_contacts = 233
native_contact_count = 559
long_range_precision = 0.769912
long_range_recall = 0.901554
```

## Safety

```text
alphafold_used_by_this_script = false
msa_used_by_this_script = false
native_truth_used_before_selection = false
coordinate_truth_used_before_selection = false
native_truth_attached_after_selection_for_evaluation = true
script_safety_rejection = none
raw_sequence_exposed_in_persisted_artifacts = false
query_fasta_persisted = false
GIF generation = not used
```

## Exact interpretation

Yes: 4AKE was solved without AlphaFold and without MSA, **if the accepted goal is MSA-free sequence-to-structure prediction audited against native contacts**.

No: the handcrafted pharmacotopology contact-collapse engine alone did not solve all-contact recall; it still needs the learned ESMFold structure as the independent global geometry prior.

The missing ingredient was therefore not another threshold. It was a learned single-sequence global geometry prior.

## What is still needed

1. Preserve this result as the first real solved 4AKE no-AlphaFold/no-MSA run.
2. Run Chai-1/Boltz later as independent replication, not as a prerequisite for this verdict.
3. Update the paper/README wording to separate:
   - `direct_msa_free_structure_solved = true`
   - `contact_ensemble_collapse_solved = false`
   - `overall_msa_free_folding_solved = true`
