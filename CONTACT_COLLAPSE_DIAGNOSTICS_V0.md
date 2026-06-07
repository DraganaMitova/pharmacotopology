# CONTACT_COLLAPSE_DIAGNOSTICS_V0

## Boundary

The frontier selector and the contact-map collapse layer are now separate.

- Frontier selection answers: **which 8x8 event regions are worth inspecting?**
- Contact collapse answers: **which residue pairs inside those regions should survive?**

Native contacts are not used to select collapsed pairs. Native labels are attached only after the collapsed set is frozen, for evaluation.

## Main production strategy

Primary strategy is now:

```text
frontier_internal_gap_balanced
```

It replaces the earlier fixed six-pair balanced budget in the main pipeline surface.

The strategy uses only internal pre-native signals:

1. ranked pair score inside each frontier event region
2. first score-gap lock for decisive one-pair events
3. short/sparse region closing to avoid reopening noisy local regions
4. five-pair internal-gap core when no decisive single-pair gap exists
5. one support-completion rectangle corner when the selected core implies a missing internally supported pair
6. one edge-rescue pair when a first-gap lock would otherwise discard a coherent boundary partner

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
precision improvement: 6.52x
```

### New internal-gap balanced collapse

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
collapse reduction: 94.6759%
```

Interpretation: the new strategy keeps the same long-range recall as the uncollapsed 1CLL frontier, but cuts the candidate map from 432 region pairs to 23 contact pairs.

## Hard-target rescue probe

The diagnostic script also reports 4AKE and 1MBN rescue probes without making them solved claims.

### 4AKE

Best current rescue probe is `ridge_coupling`:

```text
collapsed pairs: 4
true positives: 3
contact precision: 0.75
long-range precision: 1.0
long-range recall: 0.015544
```

This proves there is a precise contact signal inside the 4AKE frontier, but recall remains extremely low.

### 1MBN

Best current recall rescue probe is `frontier_recall`:

```text
collapsed pairs: 80
true positives: 9
contact precision: 0.112500
long-range precision: 0.500000
long-range recall: 0.092784
frontier long-native retention: 0.750000
```

This proves 1MBN has recoverable long-range frontier signal, but the collapse is still too noisy.

## Current conclusion

This is not a perfect contact map. It is a real contact-level collapse breakthrough for 1CLL:

```text
432 candidates -> 23 collapsed pairs
21 frontier true contacts -> 13 retained contacts
long-range recall preserved at 0.361111
contact precision raised from 0.048611 to 0.565217
```

The next unsolved frontier is not region detection for 1CLL. It is cross-protein generalization of the collapse rules, especially for 4AKE and 1MBN.
