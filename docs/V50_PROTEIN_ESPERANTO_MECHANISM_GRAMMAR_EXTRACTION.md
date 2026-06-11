# V50 Protein Esperanto Mechanism Grammar Extraction

Goal: turn V44-V49 from separate mechanism packets into one formal grammar.

## Universal Marks

`charge`, `hydrophobicity`, `aromatic_density`, `low_complexity`, `disorder_tendency`, `secondary_structure_tendency`, `motif_presence`, `proline_glycine_cysteine_effects`, `membrane_tendency`, `domain_boundary_tendency`, `interface_motif_tendency`, and `evolutionary_conservation_if_allowed`.

## Universal Operators

`closure_operator`, `repulsion_operator`, `frustration_operator`, `disorder_operator`, `phase_operator`, `interface_operator`, `membrane_pressure_operator`, `dual_basin_switch_operator`, `proteostasis_operator`, and `host_hijack_operator`.

## Mechanism Classes

`globular_closure`, `intrinsic_disorder_phase_separation`, `membrane_multidomain_folding_proteostasis`, `metamorphic_fold_switching`, `short_region_host_interface_hijacking`, `fold_upon_binding_disorder`, `cofactor_ligand_assisted_stabilization`, and `oligomerization_controlled_folding`.

## Grammar Form

Every mechanism must be expressible as:

```text
MARK + PRESSURE + OPERATOR -> STATE CHANGE -> TESTABLE EFFECT
```

V50 passes only if each hard class from V44-V48 maps into that form without special pleading.

## Evidence Boundary

Allowed before seal: pure non-coordinate evidence and explicitly tagged spatial-proxy prediction inputs.

Forbidden before seal: coordinate-derived evidence, internal runtime artifacts as biological evidence, and holdouts opened before the prediction hash exists.

## Simulation Readiness

The frozen grammar is ready for V51/V52 only when it defines transition rules, perturbation prediction rules, falsification rules, and null controls.
