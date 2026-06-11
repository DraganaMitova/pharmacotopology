# V16 Zero-MD Role Transfer Readout

This readout tests whether the locked V15/V16 role grammar can assign pressure-class roles from clean public target/context material before any MD.

It is not a new tuning panel and not a proof of solved folding.

## Locked policies

- `claim_allowed = false`
- `new_md_executed = false`
- `target_specific_threshold_tuning_allowed = false`
- `fixed_residue_cutoff_used = false`
- `native_metrics_used_for_selection = false`

## Pressure targets

- `p53_TAD_MDM2`: disorder/partner-induced binding context.
- `KcsA`: membrane/pore/complex context; biological tetramer and membrane annotation remain later requirements.
- `XCL1_lymphotactin`: two-state metamorphic/fold-switch context.

## Valid success

Success is not “folded.” Success is conservative role assignment or clean abstain without forbidden misclassification.
