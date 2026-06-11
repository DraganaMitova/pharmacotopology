# V19 KcsA Membrane/Pore Evidence Readout

V19 is a zero-MD pressure-evidence readout for KcsA. It checks whether the locked grammar can read KcsA as a membrane/pore object without treating it as a soluble compact core.

It does not run membrane MD, does not tune target-specific thresholds, does not use native metrics for selection, and does not claim the whole channel fold is solved.

Expected positive status:

```text
V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED
```

Required behavior:

```text
membrane_pore_role_evidence_found = true/false
soluble_core_misclassification_avoided = true
positive_folding_evidence_found = false
claim_allowed = false
new_md_executed = false
fixed_residue_cutoff_used = false
native_metrics_used_for_selection = false
```

The pore/filter, helix, ion, and chain-interface probes are report-only context/evidence readouts. They are not single fixed selection thresholds.
