# V36 Real Evidence Dossier Protocol

V36 moves from gate-only abstention into small, auditable scientific dossiers
for three hard protein classes:

- KcsA: membrane potassium-channel pore/filter/interface grammar.
- XCL1 / lymphotactin: metamorphic two-state state-switch grammar.
- Alpha-synuclein / SNCA: intrinsically disordered ensemble grammar.

The dossiers are evidence maps, not folding predictions. They may use external
sequence/function/disorder/family/state annotations and literature-derived state
labels. They must not use PDB coordinate-derived contact tables, RCSB coordinate
contacts, AlphaFold/ESMFold/RoseTTAFold predicted coordinates, native coordinate
metrics before selection, generated runtime reports, or MD outputs as evidence.

## Required Buckets

KcsA requires:

- `sequence_or_family_identity`
- `ion_selectivity_context`
- `membrane_topology_context`
- `filter_or_signature_context`

XCL1 / lymphotactin requires:

- `state_A_function_context`
- `state_B_function_context`
- `metamorphic_two_state_context`
- `no_mixed_state_pooling_rule`

Alpha-synuclein / SNCA requires:

- `intrinsic_disorder_context`
- `ensemble_context`
- `disorder_to_order_context`
- `no_single_native_fold_rule`

## Interpretation

V36 readiness means each target has enough real external non-coordinate dossier
evidence to define the type of folding problem it represents. It does not mean
that a structure was predicted, MD was run, native metrics were used, or protein
folding was solved.
