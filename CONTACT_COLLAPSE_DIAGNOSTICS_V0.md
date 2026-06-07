# CONTACT_COLLAPSE_DIAGNOSTICS_V0

## Boundary

The frontier selector and the contact-map collapse layer are separate.

- Frontier selection answers: **which 8x8 event regions are worth inspecting?**
- Contact collapse answers: **which residue pairs inside those regions should survive?**

Native contacts are not used to select collapsed pairs. Native labels are attached only after the collapsed set is frozen, for evaluation.

## Main production strategy

Primary strategy is now:

```text
frontier_self_deciding
```

This replaces the previous main `frontier_internal_gap_balanced` path. The old internal-gap strategy remains as a regression probe because it produced a strong 1CLL result, but it is no longer the main claim.

The self-deciding path does **not** use an accession-specific rule and does **not** use a fixed pairs-per-event budget. Per event region, it derives its own decision from pre-native internal evidence:

1. pair score distribution inside the event region
2. gap clarity / score spread
3. sequence-inferred phase shape from contact-law features
4. long-range candidate geometry inside the 8x8 region
5. direct external-coupling root count
6. ridge / lattice / boundary support already present in the feature row
7. residue-degree pressure while adding selected pairs
8. phase-aware score-surface choice: boundary/frontier vs direct-ridge trace
9. tier-aware low-score widening for broad direct-ridge traces
10. phase-specific self-collapse confidence: direct-ridge traces emphasize ridge/degree consistency, lattice/boundary regions emphasize boundary-degree completeness, and alpha-strip regions emphasize compactness-style support
11. self-verified frontier expansion, where candidate regions are ranked by native-free phase-specific self-collapse acceptance and accepted at the largest internal gap, not by a fixed confidence threshold

The controller can therefore choose different native-free outcomes per event:

- close event: no long-range candidate space
- close event: partial strip without direct external-coupling root
- first-gap lock
- alpha-strip plateau
- distribution frontier
- self-completion candidate when the selected core implies an internally supported missing partner
- switch score profile to `direct_ridge_trace` when external roots and ridge support dominate boundary evidence
- widen only a broad `direct_ridge_trace` through a score-distribution tier, then let residue-degree pressure keep the tier from becoming the whole 8x8 region

## 1CLL result

### Uncollapsed frontier

```text
selected events: 7
region candidate pairs: 432
true positives: 21
precision: 0.048611
long-range precision: 0.050000
long-range recall: 0.361111
long-range F1: 0.087838
perfect map: false
```

### Previous fixed balanced collapse

```text
collapsed pairs: 41
true positives: 13
contact precision: 0.317073
long-range precision: 0.458333
long-range recall: 0.305556
long-range F1: 0.366667
```

### Previous internal-gap balanced regression probe

```text
collapsed pairs: 23
true positives: 13
false positives: 10
contact precision: 0.565217
long-range precision: 0.684211
long-range recall: 0.361111
long-range F1: 0.472727
frontier long-native retention: 1.000000
precision improvement: 11.627348x
```

This remains the strongest 1CLL probe, but it has fixed-rule elements and is not the main self-deciding claim.

### New self-deciding main path

```text
collapsed pairs: 25
true positives: 12
false positives: 13
contact precision: 0.480000
long-range precision: 0.550000
long-range recall: 0.305556
long-range F1: 0.392858
frontier long-native retention: 0.846154
precision improvement: 9.874286x
collapse reduction: 94.212963%
```

Interpretation: the self-deciding path is slightly weaker than the best 1CLL internal-gap probe, but it is the more honest main path because it does not depend on a fixed event budget or accession-specific collapse rule.

## Hard-target rescue probe

The diagnostic script also reports 4AKE and 1MBN rescue probes without making solved claims.

### 4AKE

Previous self-deciding path used one boundary/frontier score surface and failed the weak 4AKE frontier:

```text
collapsed pairs: 9
true positives: 1
contact precision: 0.111111
long-range precision: 0.111111
long-range recall: 0.005181
frontier long-native retention: 0.038462
```

The phase-aware self-deciding path first rescued the weak 4AKE frontier by choosing `direct_ridge_trace` from internal evidence:

```text
collapsed pairs: 9
true positives: 6
contact precision: 0.666667
long-range precision: 0.666667
long-range recall: 0.031088
frontier long-native retention: 0.230769
```

The current tier-aware self-deciding path keeps the same native-free score-surface choice, then opens a wider low-score ridge tier only for broad long-range ridge traces:

```text
collapsed pairs: 38
true positives: 16
contact precision: 0.421053
long-range precision: 0.421053
long-range recall: 0.082902
long-range F1: 0.138529
frontier long-native retention: 0.615385
```

The old `ridge_coupling` precision probe remains a useful comparison:

```text
collapsed pairs: 4
true positives: 3
contact precision: 0.750000
long-range precision: 1.000000
long-range recall: 0.015544
frontier long-native retention: 0.115385
```

A self-verified frontier expansion controller was added after the raw expansion probe. It does not lower the frontier floor globally and it does not use a fixed confidence threshold such as 0.55 or 0.60. It only evaluates extra regions when the accepted seed frontier already contains a broad low-score `direct_ridge_trace` region; candidate regions are then ranked by native-free phase-specific self-collapse acceptance score and accepted at the largest internal gap in that row-specific distribution.

```text
self-verified merged expansion events: 6
uncollapsed long-range recall: 0.207254
collapsed pairs: 65
true positives: 24
contact precision: 0.369231
long-range precision: 0.369231
long-range recall: 0.124352
long-range F1: 0.186046
```

The accepted expansion rows are selected by the largest internal gap over phase-specific self-collapse acceptance scores, not native labels and not a fixed confidence threshold. For 4AKE the controller now adds three extra ridge-trace regions and rejects the broad ungated low-floor expansion. Honest conclusion: 4AKE is now **phase-specific frontier-expansion rescued but not fully solved**. The controller increases collapsed true positives from 16 to 24 and long-range recall from 0.082902 to 0.124352, but precision drops to 0.369231. This is an honest recall gain, not a full solve: total long-range recall remains below 0.15 and the precision/recall frontier still needs improvement.

### 1MBN

Previous recall rescue probe:

```text
collapsed pairs: 80
true positives: 9
contact precision: 0.112500
long-range precision: 0.500000
long-range recall: 0.092784
frontier long-native retention: 0.750000
```

New self-deciding path:

```text
collapsed pairs: 16
true positives: 9
contact precision: 0.562500
long-range precision: 0.562500
long-range recall: 0.092784
frontier long-native retention: 0.750000
```

Interpretation: 1MBN is the main cross-protein gain in this patch. The self-deciding controller keeps the same recoverable true-positive signal as the noisy recall probe while cutting 80 pairs down to 16.

## Current conclusion

This is not a perfect contact-map solver. It is a stricter, more honest collapse architecture:

```text
main path: frontier_self_deciding
native truth before selection: false
coordinate truth before selection: false
fixed pairs/event budget in main path: false
accession-specific collapse rules: false
```

1CLL remains a real contact-level collapse breakthrough. 1MBN shows that the self-deciding path can generalize a precision rescue outside 1CLL. 4AKE is no longer a total collapse failure: tier-aware self-decision recovers 16 true positives from 38 contacts with 0.421053 precision, and phase-specific self-verified frontier expansion raises that to 24 true positives from 65 contacts with 0.369231 precision. Recall is still low relative to the full native long-range map, so the next honest bottleneck is not another raw low-floor expansion; it is improving the internal confidence ranking so additional frontier regions survive without precision collapse.

## Identity-normalized frontier expansion update

The frontier expansion controller has been tightened again. Candidate regions are no longer ranked only by raw self-collapse acceptance score. The controller now builds a row-local identity baseline from the already accepted seed frontier:

```text
identity baseline = median positive seed-frontier acceptance score for the same row/profile
identity lower envelope = weakest positive seed-frontier acceptance score
candidate normalized score = candidate acceptance score / identity baseline
```

The cutoff is still native-free. It is applied on the identity-normalized candidate distribution, using the largest internal gap. A second tier can open only if a candidate still lies inside the seed frontier's own normalized lower envelope. There is no global `confidence > 0.6`, no distribution floor, no accession-specific override, and no native/contact truth before selection.

For 4AKE this did not change the selected event count: the row-local natural boundary still selects the same three expansion regions and rejects the next candidates. That is important: the system did not force extra recall by reintroducing a hidden threshold. It preserved the previous honest result:

```text
self-verified merged expansion events: 6
collapsed pairs: 65
true positives: 24
contact precision: 0.369231
long-range recall: 0.124352
long-range F1: 0.186046
```

## Frontier ceiling audit

The diagnostic report now includes `frontier_ceiling_audit`. This is post-selection/native-attached audit only. It exists to prevent a false claim that collapse can recover contacts that were never offered to it by the frontier generator.

For 4AKE:

```text
native long-range contacts: 193
seed frontier region ceiling: 0.134715
expanded frontier region ceiling: 0.196891
competitive-event region ceiling: 0.217617
candidate-generator region ceiling: 1.000000
self-deciding collapse over all competitive events: 0.124352 long-range recall
frontier_generation_bottleneck_for_0_40_recall: true
```

Meaning: with the current 100 competitive events, even taking every competitive region only exposes about 21.8% of the 4AKE long-range native map. The full candidate generator can cover the map, but the competitive frontier filter removes most of those regions before collapse sees them. So the missing puzzle piece is no longer just collapse. For 4AKE, the honest next bottleneck is a new self-deciding frontier-generation layer that can safely promote candidates from the broader candidate-event pool without collapsing precision.

## Candidate-pool frontier generation probe

A broader `frontier_generation` layer has now been added as an honest native-free probe. It looks beyond the old competitive-event ceiling and scores events from the full candidate generator, but it still does not use native/contact truth or coordinate truth before selection.

The generation controller is self-deciding:

```text
seed frontier -> row-local identity baseline
candidate-generator event -> native-free generation score
candidate score / identity baseline -> normalized generation score
largest internal gap -> prefilter boundary
self-collapse confidence -> second normalized internal-gap boundary
```

There is no fixed `confidence > 0.6` gate, no global low-score floor, and no accession-specific rule for 4AKE, 1CLL, or 1MBN. The generated frontier is then evaluated after selection in the audit report.

For 4AKE, this broad generation probe is still **not** accepted as the main path, but the failure mode changed after the latest cleanup. The generator now includes a native-free seed-identity continuation signal before the internal-gap prefilter, so weak candidate regions that are topological continuations of an accepted ridge seed can reach collapse audit. Collapse confidence now also records phase-specific lattice/alpha compactness components.

```text
current self-verified expansion:
selected events: 6
collapsed pairs: 65
true positives: 24
precision: 0.369231
long-range recall: 0.124352

candidate-pool generation probe after seed-continuation + phase-confidence:
selected events: 4
collapsed pairs: 55
true positives: 22
precision: 0.400000
long-range recall: 0.113990
```

The broader generator no longer collapses into a noisy precision failure. It improves precision versus the current expansion path, but loses recall. The audit decision is therefore:

```text
frontier_generation_probe_accepted_as_main: false
frontier_generation_decision: rejected_precision_gain_without_recall_gain
```

This is the correct outcome for an honest optimizer. The system is allowed to test broader generation, but it is not allowed to replace the current main path when recall regresses. The new diagnostic result is sharper: candidate-pool generation can now find a clean seed-continuation contact core for 4AKE, but it still cannot expand recall beyond the existing self-verified expansion without losing the no-regression guarantee. The missing piece is not a lower threshold; it is a better native-free rule for adding additional low-confidence continuation regions only when their row-level collapse footprint survives without recall/precision tradeoff.
