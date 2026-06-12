# E69 Multidomain Allosteric Architecture Grammar

E69 adds the `multidomain_allosteric_architecture` mechanism after V74 exposed `multidomain_allostery` as the dominant missing Esperanto word.

E69 also repairs a token-boundary bug in the generic-complex guard: `complex` now has to appear as a standalone word, so explicit `low complexity` disorder/phase evidence is no longer flattened into a generic-complex abstention.

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

Run `V75_MULTIDOMAIN_ALLOSTERY_REPAIR_PANEL_200` to test whether E69 closes the V74 multidomain/allostery failures without absorbing the remaining coiled/repeat and disulfide/secretory signals.
