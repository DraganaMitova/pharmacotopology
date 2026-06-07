# Event-region contact collapse diagnostics V0

This patch separates two different claims that were previously being mixed:

1. **frontier/event-region recovery**: selected 8x8 segment regions contain native contacts;
2. **contact-map recovery**: selected residue-pair contacts are correct.

For 1CLL, `7/7` frontier events is not a perfect contact map. The uncollapsed frontier contains 432 residue-pair candidates, 21 true positives and 411 false positives. That is precision `0.048611`, recall `0.064220`, and `perfect_contact_map = false`.

## New collapse layer

Added `event_region_contact_collapse_v0` in:

- `src/pharmacotopology/folding_event_region_contact_collapse.py`
- `scripts/run_contact_collapse_diagnostics_v0.py`
- `tests/test_event_region_contact_collapse_v0.py`

The layer takes each selected event region and ranks residue pairs using only pre-evaluation signals:

- external coupling density / local DCA support;
- sequence-law support;
- ridge coherence;
- region-boundary coherence;
- residue degree cap.

Native/contact truth is attached only after collapse for evaluation.

## Generated artifacts

- `first_contact_clean_pharmacotopology_layer_run/contact_collapse_diagnostics_v0.json`
- `first_contact_clean_pharmacotopology_layer_run/1cll_contact_collapse_v0.json`
- `first_contact_clean_pharmacotopology_layer_run/contact_collapse_pairs_v0.csv`
- `first_contact_clean_pharmacotopology_layer_run/contact_collapse_event_rows_v0.csv`

## 1CLL result

Uncollapsed frontier:

- selected frontier events: `7`
- predicted region pairs: `432`
- true positives: `21`
- precision: `0.048611`
- recall over all native contacts: `0.064220`
- long-range recall: `0.361111`

Precision collapse (`frontier_precision`):

- collapsed pairs: `21`
- true positives: `9`
- precision: `0.428571`
- precision improvement: `8.816338x`
- total native recall: `0.027523`
- long-range recall: `0.250000`
- status: better precision, but not solved recall.

Recall collapse (`frontier_recall`):

- collapsed pairs: `111`
- true positives: `18`
- precision: `0.162162`
- frontier-native retention: `0.857143`
- frontier-long-native retention: `1.000000`
- long-range recall: `0.361111`, matching the uncollapsed frontier long-range recall.
- status: keeps almost all recoverable frontier signal while still removing 321/432 region pairs.

This means collapse is now explicit and measurable. The system no longer treats all 64 pairs inside an 8x8 event region as contacts.

## 4AKE / 1MBN note

The collapse diagnostic also exposes that one strategy does not fit every protein:

- `1MBN:A`, ridge-coupling collapse: precision `0.375000`, but low retention.
- `4AKE:A`, ridge-coupling collapse: precision `0.750000`, but low retention.
- `1TIM:A` remains a guard case: aggressive recall collapse does not help and should not be interpreted as solved.

So the next real front is not another global threshold. It is a row/event adaptive collapse strategy that decides when to use sparse-boundary collapse, ridge-coupling collapse, recall mode, or abstain.

## Honest conclusion

This patch does **not** claim the folding problem is solved.

It cracks the previous ambiguity: `7/7 frontier events` and `perfect contact map` are now separated in code, metrics, reports and tests. For 1CLL, the precision collapse can raise contact precision from `0.048611` to `0.428571`, while the recall collapse can retain `0.857143` of native contacts available inside the frontier. The remaining unsolved problem is the adaptive choice of collapse mode without using native truth.
