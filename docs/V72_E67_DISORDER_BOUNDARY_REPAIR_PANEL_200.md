# V72 E67 Disorder Boundary Repair Panel 200

V72 was the repair panel for V71's dominant failure:

```text
disorder_misread: 31
```

E67 added:

```text
disorder_boundary_and_fold_upon_binding
```

New words:

```text
IDR_boundary
structured_domain_plus_IDR_tail
fold_upon_binding_region
phase_prone_low_complexity
flexible_loop_not_disorder
disorder_with_local_motif
```

Panel composition:

```text
V71 disorder failure replay: 31
disorder-boundary positive expansion: 60
generic oligomer controls: 25
true TM / disorder conflict sentinels: 20
assembly-required / disorder priority sentinels: 14
metal-ligand / disorder conflict sentinels: 20
V71 closed-beta tracking replay: 30
```

Result:

```text
status: V72_E67_DISORDER_BOUNDARY_REPAIR_PANEL_PASSED_REVIEW_REQUIRED
targets_total: 200
accepted_count: 200
accepted_supported: 170
failed_accepted: 30
targeted_failed_accepted: 0
accepted_accuracy: 0.85
coverage: 1.0
controls_passed: true
sentinel_regressions: 0
```

Repair result:

```text
V71 disorder failures repaired: 31 / 31
disorder expansion supported: 60 / 60
generic oligomer controls preserved: 25 / 25
true TM preserved: 20 / 20
assembly-required preserved: 14 / 14
metal-ligand preserved: 20 / 20
```

The 30 failed accepted rows are intentionally tracking-only:

```text
closed_beta_tracking_remaining: 30
```

They are not E67 disorder repair failures. They define the next repair target:

```text
V73_BETA_CLOSURE_TOPOLOGY_REPAIR_PANEL_200
```
