# V37 Mechanism Question Probes Protocol

V37 tests whether V36 non-coordinate evidence dossiers can assign the correct
folding-problem mechanism grammar without using target names as the only signal.

The expected classes are:

- KcsA: `membrane_pore_filter_oligomeric_ion_selectivity`
- XCL1 / lymphotactin: `metamorphic_two_state_fold_switch`
- Alpha-synuclein / SNCA: `intrinsic_disorder_contextual_ensemble`

V37 is not de novo folding prediction, structure prediction, MD, or a
folding-solved claim. It maps the class of mechanism problem and the claims that
must remain forbidden.

## Controls

The probe must pass baseline assignment, target-name masking, swapped-dossier
detection, partial evidence controls, forced wrong-grammar controls, generic
annotation controls, coordinate leakage controls, internal runtime leakage
controls, and placeholder source controls.

## Allowed Inputs

Production V37 reads only:

- `data/external_evidence_dossiers/KcsA/evidence_dossier.json`
- `data/external_evidence_dossiers/XCL1_lymphotactin/evidence_dossier.json`
- `data/external_evidence_dossiers/alpha_synuclein_SNCA/evidence_dossier.json`

It does not read coordinate contacts, V33/V34 KcsA coordinate CSVs, native
structure metrics, predicted coordinate files, or runtime reports as evidence.
