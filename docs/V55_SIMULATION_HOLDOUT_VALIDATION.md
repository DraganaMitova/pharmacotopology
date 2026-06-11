# V55 Simulation Holdout Validation

V55 validates sealed operator trajectories after prediction hashes exist.

## Allowed Post-Seal Evidence

Crosslinks, FRET, NMR chemical shifts, SAXS, EPR/DEER, DCA couplings, mutagenesis, biochemical rescue data, and coordinate-derived sources only after sealing.

## Validation Question

```text
Did the simulated operator trajectory predict the same regions, states, and effects that independent evidence later supports?
```

The validation question is not whether a perfect RMSD structure was produced.

## Required Guard

If a holdout is opened before the sealed hash exists, the run blocks. Failed predictions remain failed; they are not repaired after holdout opening.
