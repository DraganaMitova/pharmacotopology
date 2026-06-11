# V56 Operator To Custom Force Bridge

V56 defines the physical-simulation bridge without making it the first proof.

## Bridge Rule

```text
Protein Esperanto operator
-> custom bias / restraint / field term
-> OpenMM comparison against unbiased baseline
```

## Example Mapping

- `closure_operator` -> weak attractive bias between predicted regions.
- `disorder_operator` -> penalty against over-collapse.
- `interface_operator` -> partner-proximity or surface-readiness bias.
- `dual_basin_switch_operator` -> two-state basin potential.
- `membrane_pressure_operator` -> environment-dependent burial/exposure term.
- `host_hijack_operator` -> short motif host-interface bias.

## Guard

Custom forces must not encode holdout coordinates or force the answer. The coarse sealed operator trajectory remains the first proof layer.
