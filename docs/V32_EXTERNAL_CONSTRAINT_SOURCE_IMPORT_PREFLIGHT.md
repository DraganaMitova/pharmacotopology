# V32 External Constraint Source Import Preflight

V31 clean-abstained because XCL1 and KcsA candidate files in the selected scope were generated internal reports, not real external constraints. V32 is therefore **not** a constraint-backed readout. It is an import/acquisition gate for real external constraint files.

## Boundary

V32 does not run MD, does not tune thresholds, does not use native metrics for selection, and does not allow folding claims.

Internal runtime artifacts remain audit-only:

```text
first_contact_clean_pharmacotopology_layer_run/* -> generated_internal_report -> audit only
```

External annotations remain context-only, not constraints:

```text
annotation-only -> role context only
```

Only real imported files under `data/external_constraints/<target>/...`, with explicit source/provenance rows in `data/external_constraints/v32_external_constraint_source_import_manifest.json`, can support a future constraint-backed readout.

## Target requirements

### XCL1_lymphotactin

Requires both state buckets before a state-switch constraint-backed readout:

- state_A / chemokine-like / monomer constraint
- state_B / dimer / beta-sandwich-like constraint
- no mixed-state pooling

### KcsA

Requires both pore/coupling and assembly/interface context before a KcsA constraint-backed readout:

- pore/filter/selectivity or external coupling/contact support
- assembly/interface/tetramer/chain-interface context constraint
- no whole-channel/tetramer fold claim from annotation alone

## Outcomes

If no real external constraints are imported:

```text
V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED
```

If real external constraints satisfy a selected target requirement:

```text
V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED
```

If internal generated reports are supplied as external evidence:

```text
V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_BLOCKED
```
