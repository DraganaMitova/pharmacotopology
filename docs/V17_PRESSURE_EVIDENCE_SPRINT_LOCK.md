# V17 Pressure Evidence Sprint Lock

V17 is a one-shot sprint stage. It combines:

1. evidence manifest
2. evidence file preflight
3. zero-MD evidence-readiness readout
4. one selected V18 target decision

It does **not** run MD, tune thresholds, use native metrics for selection, use a fixed residue cutoff, or upgrade V16 role classification into folding evidence.

Expected decision when the V16 pressure structures are present:

```text
selected_V18_target = p53_TAD_MDM2
selected_V18_test = p53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST
```

p53/MDM2 is selected first because it is the fastest clean evidence test for the disorder/partner-induced role class: isolated autonomous fold claims remain forbidden, while partner-bound interface/helix evidence can be read out without heavy MD.
