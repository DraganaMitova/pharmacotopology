# LOCAL_COEVOLUTION_EXPANSION_V0

This batch tests the local co-evolution hypothesis around sparse DCA anchors:

```text
safe DCA anchors
→ +/- window around each anchor
→ local co-evolution support
→ selected contacts
→ native audit only after selection
```

Raw MSA files are not bundled in this ZIP. Therefore the implemented default is an honest proxy channel, not a new MI calculation:

```text
raw_msa_available_for_true_local_mi = false
local_mi_channel_is_proxy_not_new_msa_calculation = true
```

The proxy uses only safe external coupling rows and sequence-only chemistry/secondary-structure terms:

```text
nearby external APC/raw support
anchor-conditioned local support
local patch support
chemistry/secondary-structure support
```

No native contacts, native coordinates, ESMFold, AlphaFold, or learned geometry prior are used before selection.

## Focused 4AKE result

```text
precision = 0.296875
recall = 0.033989
long_range_recall = 0.098446
F1 = 0.060995
accepted_local_pair_count = 14
safe_anchor_count = 50
folding_problem_solved = false
```

## 8-row benchmark result

```text
mean_precision = 0.167227
mean_recall = 0.031057
mean_long_range_recall = 0.138154
mean_F1 = 0.051568
mean_accepted_local_pair_count = 19.375
folding_problem_solved = false
```

Interpretation: local co-evolution proxy around DCA anchors improves precision compared with loose dense-geodesic expansion, but it does not create enough true long-range coverage. Because raw MSA is absent, this is not a definitive rejection of true local-MSA MI; it is a rejection of the available safe proxy channel.
