# Evidence Taxonomy

V49 freezes the evidence boundary for reviewer-readable source separation. Every new source row must say what kind of evidence it is, what role it plays, and whether it can be used before prediction sealing.

## Categories

### `pure_non_coordinate`

Sequence and annotation evidence that does not directly encode distances, contacts, coordinates, or structural answer keys.

Examples: sequence, UniProt features, Pfam/InterPro, DisProt, disorder annotations, PTM annotations, function annotations, motif labels, composition, charge, hydrophobicity, low-complexity, and secondary-structure propensity.

Allowed roles: `prediction_input`, `holdout_validation`, `baseline_only`, or `blocked`, depending on target protocol.

### `spatial_proxy_non_coordinate`

Evidence that is not coordinate data but can encode spatial information.

Examples: DCA/evolutionary couplings, FRET, crosslinking, NMR chemical shifts, SAXS, EPR/DEER, distance restraints, and contact-like biochemical constraints.

Rule: `spatial_proxy_non_coordinate` sources must be explicitly tagged. They must not be represented as `pure_non_coordinate` evidence. Each source must declare whether it is used as `prediction_input`, `holdout_validation`, `baseline_only`, or `blocked`.

### `coordinate_derived`

Evidence derived from explicit structures, coordinate models, native contacts, or structural alignments.

Examples: PDB/mmCIF coordinates, contact maps from structures, AlphaFold/ESMFold/RoseTTAFold coordinates, native distance matrices, structural alignments, and PDB-derived interface contacts.

Rule: `coordinate_derived` sources are blocked before prediction sealing. They are allowed only after sealing as scoring holdouts when explicitly marked `holdout_validation`.

### `internal_runtime`

Project-generated runtime artifacts.

Examples: `first_contact_clean_pharmacotopology_layer_run` reports, certificates, generated summaries, previous prediction outputs, and validation summaries from this repository.

Rule: `internal_runtime` artifacts are allowed only as audit metadata. They are never biological prediction evidence and never biological claim evidence.

## Required Source Role Fields

All new source manifest rows must include these fields:

```text
source_class
source_role
spatial_proxy
coordinate_derived
internal_runtime
allowed_for_prediction
allowed_for_holdout
allowed_for_claim
leakage_risk
rationale
```

`source_role` must be one of:

```text
prediction_input
holdout_validation
baseline_only
blocked
audit_metadata
```

## Example Row

```json
{
  "source_id": "EXAMPLE_UNIPROT_SEQUENCE",
  "source_class": "pure_non_coordinate",
  "source_role": "prediction_input",
  "spatial_proxy": false,
  "coordinate_derived": false,
  "internal_runtime": false,
  "allowed_for_prediction": true,
  "allowed_for_holdout": false,
  "allowed_for_claim": true,
  "leakage_risk": "low",
  "rationale": "Sequence and feature annotations do not expose coordinates, native contacts, or post-seal holdouts."
}
```

## Boundary Checks

- A source with `spatial_proxy: true` must have `source_class: spatial_proxy_non_coordinate`.
- A source with `coordinate_derived: true` must not be allowed for prediction before sealing.
- A source with `internal_runtime: true` must not be allowed for biological prediction or biological claim support.
- A source used as `holdout_validation` must not be opened before the sealed prediction hash exists.
