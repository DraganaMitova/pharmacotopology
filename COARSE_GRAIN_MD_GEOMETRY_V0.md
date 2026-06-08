# COARSE_GRAIN_MD_GEOMETRY_V0

This batch tests the structure-first pivot:

```text
sequence / safe DCA restraints -> global C-alpha MD-style geometry -> extracted contacts -> native audit
```

It is intentionally not:

```text
contacts -> factor graph -> structure
```

## Boundary

- No AlphaFold/ESMFold.
- No template structure.
- No trained geometry prior.
- No coordinate/native truth before selection.
- Native contacts are attached only after the final geometry emits a contact map.
- No OpenMM/GROMACS dependency in this bounded implementation.
- No GIF generation.
- No predictor subprocesses.
- No full suite requirement.

## What this is

A dependency-free coarse-grain MD-style falsifier. It moves the entire C-alpha trace globally using:

- chain springs,
- lightweight secondary-structure springs,
- DCA distance restraints when available,
- excluded volume,
- weak compaction,
- damped velocity updates and a deterministic temperature schedule.

Then it extracts contacts from the final geometry with an 8 Å cutoff and audits against native contacts.

## Claim gate

The method may claim a solve only if every row clears:

- precision >= 0.70,
- recall >= 0.70,
- long-range recall >= 0.70,
- matched-control F1 and long-range margins >= 0.15,
- no coordinate/native/template/learned-prior taint.

Pure sequence mode is required for a universal physical law claim.

## Bounded run result

Command:

```bash
PYTHONPATH=src timeout 180s python3 scripts/run_coarse_grain_md_geometry_v0.py --mode all --md-steps 72 --restarts 3 --max-restraints 96
```

Summary:

| mode | precision | recall | long-range recall | F1 | solved |
|---|---:|---:|---:|---:|---|
| pure sequence global MD geometry | 0.372736 | 0.380634 | 0.000000 | 0.372930 | false |
| external DCA restrained global MD geometry | 0.108041 | 0.441618 | 0.238876 | 0.172420 | false |
| external DCA multistart global MD geometry | 0.108425 | 0.444780 | 0.235371 | 0.173142 | false |

4AKE row under multistart DCA MD geometry:

| protein | extracted contacts | extracted LR contacts | precision | recall | long-range recall | F1 | solved |
|---|---:|---:|---:|---:|---:|---:|---|
| 4AKE:A | 1738 | 768 | 0.139240 | 0.432916 | 0.134021 | 0.210722 | false |

## Interpretation

The structure-first global relaxation changed the failure mode:

- pure sequence global geometry gives decent local/core contact F1, but zero long-range topology;
- DCA restraints inject long-range coverage, but the extracted contact map becomes too noisy and precision collapses;
- multistart does not rescue the precision/long-range tradeoff.

So the missing object is not merely "move everything globally." The missing object is a global geometry field that is simultaneously:

```text
globally relaxed + long-range aware + precision preserving
```

This bounded dependency-free MD-style layer did not replace the learned ESMFold geometry prior.
