# V16 Transfer Data Preflight

V16 transfer testing starts with data/material preflight, not MD.

This layer checks whether the locked V16 pressure targets have enough clean public target/context material to proceed to a zero-MD role-transfer readout.

It does not:

- run molecular dynamics
- tune thresholds
- change the V15 grammar
- use native metrics for selection
- claim solved folding

## Locked targets

- `p53_TAD_MDM2` — disorder / partner-induced binding pressure
- `KcsA` — membrane / pore / oligomer pressure
- `XCL1_lymphotactin` — metamorphic / fold-switching pressure

## Expected statuses

If target material has not been downloaded yet, the expected status is:

```text
V16_TRANSFER_DATA_PREFLIGHT_BLOCKED_MISSING_REQUIRED_MATERIAL
```

After public RCSB structures are downloaded and provenanced, the expected status is:

```text
V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT
```

This still does not mean protein folding is solved. It only means the next role-readout layer can be run without mixing missing-data failures with grammar failures.
