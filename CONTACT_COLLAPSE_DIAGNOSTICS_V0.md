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
- `tests/test_contact_collapse_pipeline_integration_v0.py`

The layer takes each selected event region and ranks residue pairs using only pre-evaluation signals:

- external coupling density / local DCA support;
- sequence-law support;
- ridge coherence;
- region-boundary coherence;
- residue degree cap.

Native/contact truth is attached only after collapse for evaluation. The production-facing strategy is `frontier_balanced`: it uses the internal collapse score and a fixed `6 pairs/event` budget. The six-pair budget came from the non-oracle fixed-budget probe on the locked 1CLL frontier, not from selecting individual true contacts.

## Pipeline integration

Collapse is no longer only a diagnostic script. It is wired into the main selector/trace-loop surfaces:

- `src/pharmacotopology/folding_coupling_nucleus_selector.py`
- `src/pharmacotopology/folding_external_coupling_trace_loop.py`

The integrated pipeline now exposes contact-collapse summaries and writes collapse artifacts next to the existing selector/frontier files.

New selector-side artifacts:

- `first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_contact_collapse_rows.csv`
- `first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_collapsed_contacts.csv`
- `first_contact_clean_pharmacotopology_layer_run/coupling_nucleus_selector_contact_collapse_events.csv`

New external-frontier artifacts:

- `first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_contact_collapse_rows.csv`
- `first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_collapsed_contacts.csv`
- `first_contact_clean_pharmacotopology_layer_run/external_coupling_trace_loop_contact_collapse_events.csv`

## Generated diagnostic artifacts

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
- long-range precision: `0.600000`
- long-range recall: `0.250000`
- long-range F1: `0.352941`
- status: strong precision, but too much recall loss.

Recall collapse (`frontier_recall`):

- collapsed pairs: `111`
- true positives: `18`
- precision: `0.162162`
- frontier-native retention: `0.857143`
- frontier-long-native retention: `1.000000`
- long-range precision: `0.200000`
- long-range recall: `0.361111`, matching the uncollapsed frontier long-range recall.
- long-range F1: `0.257426`
- status: keeps almost all recoverable frontier signal while still removing 321/432 region pairs.

Balanced collapse (`frontier_balanced`, fixed 6 pairs/event):

- collapsed pairs: `41`
- true positives: `13`
- false positives: `28`
- precision: `0.317073`
- precision improvement: `6.523229x`
- long-range precision: `0.458333`
- long-range recall: `0.305556`
- long-range F1: `0.366667`
- status: first non-oracle tradeoff where contact precision and long-range recall both clear the `0.30` line on 1CLL.

This means collapse is now explicit and measurable. The system no longer treats all 64 pairs inside an 8x8 event region as contacts.

## 4AKE / 1MBN note

The collapse diagnostic also exposes that one strategy does not fit every protein:

- `1MBN:A`, ridge-coupling collapse: precision `0.375000`, but low retention.
- `4AKE:A`, ridge-coupling collapse: precision `0.750000`, but low retention.
- `1TIM:A` remains a guard case: aggressive recall collapse does not help and should not be interpreted as solved.

So the next real front is not another global region threshold. It is row/event adaptive collapse policy: when to use balanced frontier collapse, sparse precision collapse, ridge-coupling collapse, recall mode, or abstain.

## Honest conclusion

This patch does **not** claim the folding problem is solved.

It cracks the previous ambiguity: `7/7 frontier events` and `perfect contact map` are now separated in code, metrics, reports and tests. For 1CLL, the balanced collapse reduces the candidate space from `432` to `41`, raises contact precision from `0.048611` to `0.317073`, and keeps long-range recall at `0.305556`. The remaining unsolved problem is generalizing the collapse policy across proteins without using native truth.
