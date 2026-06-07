# SELF_CONSISTENT_FOLDING_V0

This batch adds a native-free self-consistent contact-loop probe.

It does **not** treat self-generated contacts as external independent evidence. The old external-source contract stays closed: AlphaFold/PDB/RoseTTAFold-like external structures are still the only way to permit an external benchmark claim in DONE mode.

## What the new layer does

For a locked target, the runner:

1. Builds the candidate contact space from the existing coupling trace-loop event frontier.
2. Seeds the loop only with external evolutionary coupling contacts that land inside that candidate space.
3. Iteratively re-scores candidate contacts using rank-normalized internal signals:
   - direct coupling support,
   - prior loop presence,
   - endpoint graph anchoring,
   - shared-neighbor closure,
   - event score rank,
   - event recurrence rank.
4. Cuts the selected contact set by a row-local largest-gap boundary, with the seed evidence acting as the identity envelope.
5. Builds matched negative controls by scrambling the seed/coupling contacts while preserving sequence-separation buckets and confidence ordering.
6. Lets the run decide by rank/gap against those controls.

There is no static confidence threshold such as `0.7`.

## 4AKE result without AlphaFold/PDB

Command:

```bash
python3 scripts/run_self_consistent_contact_loop_v0.py \
  --source-accession 4AKE:A \
  --event-source aggressive_relative_frontier
```

Observed in this batch:

- candidate events: `1078`
- candidate pairs: `19824`
- external coupling seed pairs inside frontier: `197`
- final self-consistent pairs: `197`
- internal self-consistency status: `self_consistency_survived_matched_negative_controls`
- external independent claim allowed: `false`
- global folding claim allowed: `false`
- folding problem solved: `false`

Native audit after selection, not used before selection:

- final contact precision: `0.147208`
- final contact recall: `0.051878`
- final long-range recall: `0.093264`
- precision delta vs seed: `-0.030457`
- long-range recall delta vs seed: `-0.005182`

Interpretation: the loop can stabilize an internal self-source, but this self-source does **not** improve 4AKE folding/contact accuracy without external structure evidence. So this is a real harness, not a fake “solved” patch.

## Verification

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest -q
# 2 passed
```

The environment used for this batch had external pytest plugins that hung on plain `python3 -m pytest -q`; disabling plugin autoload verified the project tests themselves.
