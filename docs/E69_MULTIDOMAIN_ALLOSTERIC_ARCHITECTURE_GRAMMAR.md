# E69 Multidomain Allosteric Architecture Grammar

E69 adds the `multidomain_allosteric_architecture` mechanism after V74 exposed `multidomain_allostery` as the dominant missing Esperanto word.

E69 also repairs a token-boundary bug in the generic-complex guard: `complex` now has to appear as a standalone word, so explicit `low complexity` disorder/phase evidence is no longer flattened into a generic-complex abstention.

V75 proved the E69 repair under the self-decision cortex:

```text
failed_accepted = 0
accepted_accuracy = 1.0
accepted_supported = 187 / 187
V74 multidomain/allostery failures repaired = 33 / 33
real physical calibration inputs loaded = 8 locked RCSB coordinate rows
physical basis claim allowed = false
```

The same self-decision judge cleanly abstained on disulfide, coiled-coil, and repeat/solenoid candidate words instead of accepting wrong predictions.

V75 also keeps one multidomain positive as a conservative abstain because the endogenous operator-basis probe marked it assignment-sensitive. That is treated as honest non-acceptance, not a failed accepted prediction.

## New Words

| word | role |
| --- | --- |
| `multidomain_allostery` | overall coupled-domain basin |
| `domain_boundary` | boundary separating modular domains |
| `hinge_region` | flexible coupling axis |
| `interdomain_lock` | interface state that locks domain geometry |
| `allosteric_basin_shift` | basin shift caused by domain coupling |
| `domain_reorientation` | orientation change between domains |
| `modular_architecture` | coupled-domain packing state |
| `domain_swapping` | swapped-domain subtype |

## Priority

E69 is below the already repaired priority classes:

| protected class | reason |
| --- | --- |
| `membrane_multidomain_folding_proteostasis` | true TM topology must remain primary |
| `assembly_required_folding` | partner-completed core stays distinct from multidomain allostery |
| `metal_cluster_and_ligand_locked_basin` | metal/ligand basin locking stays primary |
| `disorder_boundary_and_fold_upon_binding` | IDR boundary stays primary |
| `beta_closure_topology` | closed-beta topology stays primary |
| `intrinsic_disorder_phase_separation` | explicit low-complexity/phase evidence must not be blocked by the generic-complex guard |

## Next

Run `V76_SECRETORY_DISULFIDE_REPAIR_PANEL_200` after adding `E70_SECRETORY_DISULFIDE_REDOX_TOPOLOGY_GRAMMAR`, because V75 mined `disulfide_secretory_redox_context` as the loudest remaining missing word.
