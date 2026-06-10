# V13a 1UBQ Adaptive Chemical Policy V0

This layer is a claim-safe postprocess readout for completed V13a 1UBQ trajectories.
It exists because the strict/provenance repair showed that 1UBQ no longer fails due to missing target data, but the legacy global chemical hard gate can suppress a repeatable balanced-core signal.

## Constraints

- No synthetic target fallback.
- No partial target fallback.
- No native fallback.
- No protein-name-specific threshold.
- No native precision used for threshold selection.
- No selector role ontology rewrite.
- `claim_allowed=false` until an official purpose-aware runtime is reviewed.

## Policy

For `single_domain_compact` targets, chemical evidence is recorded and reported, but the adaptive soft guard does not allow a single global chemical cutoff to kill a DCA-supported, topology-compatible, replica-supported balanced core when the other guards are clean.

The readout still requires:

- nonzero balanced core,
- DCA support/enrichment,
- zero noise outside effective anchor set,
- no long-range evidence pollution,
- complete classification coverage,
- enough replica support.

## Outputs

Main command:

```bash
bash "$REPO_ROOT/scripts/run_v13a_1ubq_adaptive_chemical_policy_readout.sh"
python3 "$REPO_ROOT/scripts/print_v13a_adaptive_chemical_policy_summary.py"
```

Main output directory:

```text
first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_ADAPTIVE_CHEMICAL_POLICY_READOUT/
```

Main certificate:

```text
v13a_1ubq_purpose_gate_readout_certificate.json
```

The certificate includes:

- selected adaptive soft frequency band,
- hard chemical threshold sweep,
- selected hard chemical admissibility band if any,
- chemical score by selected pair,
- corrected mean frequency by selected pair, guaranteed to be a contact-presence frequency in [0, 1].

## Interpretation

If the adaptive soft guard finds a core while the hard chemical sweep only passes at very low or zero threshold, the correct status is:

```text
chemical_gate_policy_mismatch_with_positive_adaptive_soft_signal
claim_allowed=false
```

This is not a solved biological transfer claim. It is a clean localization of a policy mismatch.
