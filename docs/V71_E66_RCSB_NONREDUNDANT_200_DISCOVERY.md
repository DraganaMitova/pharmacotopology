# V71 E66 RCSB Nonredundant 200 Discovery

V71 was a fresh discovery blade under E66. It used 200 fresh 30% representative RCSB protein entities not used by V61 through V70.

Shard composition:

```text
V71A broad RCSB nonredundant: 50
V71B disorder / low-complexity / flexible-region enriched: 50
V71C beta-barrel / beta-propeller / repeat / solenoid enriched: 50
V71D coiled-coil / helix-bundle / multidomain enriched: 50
```

Result:

```text
status: V71_E66_RCSB_NONREDUNDANT_200_DISCOVERY_FAILURES_MINED_REVIEW_REQUIRED
targets_total: 200
accepted_count: 181
accepted_supported: 111
clean_abstain_supported: 19
failed_accepted: 70
accepted_accuracy: 0.6132596685082873
coverage: 0.905
controls_passed: true
sentinel_regressions: 0
```

Failure taxonomy:

```text
disorder_misread: 31
closed_beta_topology: 30
membrane_topology_missed_or_misread: 5
other: 3
multidomain_allostery: 1
```

Dominant missing word:

```text
disorder_misread
```

The autopsy was specific: E66 accepted many disorder/low-complexity/flexible-region cases as `oligomerization_controlled_folding` when generic copy or assembly metadata was present. That became E67:

```text
E67_DISORDER_BOUNDARY_AND_FOLD_UPON_BINDING_GRAMMAR
```

The near-tie signal remains:

```text
closed_beta_topology: 30
```

That is not erased. It is carried into V72 as tracking-only and becomes the next repair target after E67 passes.
