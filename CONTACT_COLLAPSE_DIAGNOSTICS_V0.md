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

The controller can therefore choose different native-free outcomes per event:

- close event: no long-range candidate space
- close event: partial strip without direct external-coupling root
- first-gap lock
- alpha-strip plateau
- distribution frontier
- self-completion candidate when the selected core implies an internally supported missing partner
- switch score profile to `direct_ridge_trace` when external roots and ridge support dominate boundary evidence

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

The new phase-aware self-deciding path chooses `direct_ridge_trace` from internal evidence on the weak-ridge regions:

```text
collapsed pairs: 9
true positives: 6
contact precision: 0.666667
long-range precision: 0.666667
long-range recall: 0.031088
frontier long-native retention: 0.230769
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

Native-free frontier expansion was added as a probe, not as a solved main path:

```text
merged expansion events: 17
uncollapsed long-range recall: 0.160622
collapsed pairs: 49
true positives: 8
contact precision: 0.163265
long-range precision: 0.205128
long-range recall: 0.041451
long-range F1: 0.068966
```

Honest conclusion: 4AKE is **partially cracked at precision level** by phase-aware self-deciding collapse, but it is not fully solved. Expansion exposes slightly more recall, yet the system correctly keeps it as a probe because precision collapses when the frontier is widened.

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

1CLL remains a real contact-level collapse breakthrough. 1MBN shows that the self-deciding path can generalize a precision rescue outside 1CLL. 4AKE is no longer a total collapse failure: phase-aware self-decision recovers 6 true positives from 9 contacts with 0.666667 precision, but recall is still low. The next honest bottleneck is a frontier-expansion controller that widens candidate space only when its own post-collapse confidence does not collapse.
