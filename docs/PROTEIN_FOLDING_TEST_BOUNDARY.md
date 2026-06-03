# Protein Folding Test Boundary

This branch adds a small protein-folding hypothesis lab to the project. The
goal is not to solve folding. The goal is to create a safe benchmark shell
where topology signatures can be compared against folding-relevant reference
summaries.

The test claim is deliberately narrow:

```text
topology signatures extracted from protein sequences and structures
can be compared against folding-relevant benchmark summaries
under explicit uncertainty and failure labels
```

## Modeled Chain

```text
protein sequence
-> predicted topology signature
-> reference topology signature
-> contact-map proxy similarity
-> fold-class match
-> uncertainty radius
-> evidence readiness
```

The signature dimensions are:

```text
sequence_complexity
secondary_structure_balance
contact_map_closure
hydrophobic_core_closure
loop_disorder_pressure
domain_boundary_stability
long_range_contact_order
conformational_flexibility
knot_or_entanglement_signature
uncertainty_radius
```

## Current Benchmark Status

The default benchmark rows are placeholders. They are useful for exercising the
code path, output schema, tests, and safety boundary, but they are not external
validation.

Real calibration requires reference rows derived from external structure
resources such as contact maps, fold-class annotations, structure ensembles, or
curated folding benchmark sets.

The repository now includes a first locked 10-row external benchmark slice:

```text
data/folding_benchmarks_real_10.locked.json
```

It contains real sequences and source accessions from RCSB PDB/CATH-style
labels plus DisProt disorder references:

```text
2 alpha-rich rows
2 beta-rich rows
2 alpha-beta mixed rows
2 multidomain / boundary-sensitive rows
2 disordered / flexible rows
```

Each row includes:

```text
real sequence
source_database
source_accession
reference_structure_source
reference_label_source
is_external_reference = true
curation_notes
```

The predictor receives the sequence only. Labels, accessions, and reference
topology signatures remain benchmark-side scoring metadata.

## External Benchmark Adapter

The benchmark runner can load externally derived rows:

```bash
python3 scripts/run_folding_topology_benchmark.py --benchmark-file data/folding_benchmarks_real.json --require-external
```

The adapter validates that:

```text
the file contains reference rows
the sequence uses supported amino-acid symbols
the reference source looks external, such as pdb:, afdb:, casp:, cath:, scop:, disprot:, or external:
placeholder/example/template source labels are rejected when --require-external is used
all topology signature dimensions are present
all topology signature values are inside 0..1
no drug-design or clinical-use boundary is opened
```

The repo includes:

```text
data/folding_benchmarks_real.example.json
data/folding_benchmarks_real_10.locked.json
data/folding_benchmarks_real_500.locked.json
```

The example file is a schema template, not evidence. Copy it to
`data/folding_benchmarks_real.json` and replace the row with externally derived
structure summaries before using `--require-external`.

The 10-file is a small locked benchmark slice, not a general validation set. On
the current model it reports:

```text
external rows = 10
fold-class matches = 4
fold-class mismatches = 6
accuracy = 0.4
folding_problem_solved = false
```

The 500-file is a target shell, not a completed benchmark. It records the
intended proof ladder and current lock blockers until real rows exist.

```bash
python3 scripts/build_real_folding_benchmark_500.py \
  --size 500 \
  --output data/folding_benchmarks_real_500.locked.json \
  --lock
```

The intended stratification is:

```text
100 mostly-alpha / compact domains
100 mostly-beta / long-range-contact domains
100 alpha-beta mixed domains
100 multidomain / boundary-sensitive proteins
100 disordered or flexible proteins
```

Run the real-10 slice with:

```bash
python3 scripts/run_folding_topology_benchmark.py \
  --benchmark-file data/folding_benchmarks_real_10.locked.json \
  --require-external \
  --report-output first_contact_clean_pharmacotopology_layer_run/real_folding_10_report.json \
  --csv-output first_contact_clean_pharmacotopology_layer_run/real_folding_10_rows.csv
```

And visualize it with:

```bash
python3 scripts/render_folding_benchmark_dashboard.py \
  --report first_contact_clean_pharmacotopology_layer_run/real_folding_10_report.json \
  --csv first_contact_clean_pharmacotopology_layer_run/real_folding_10_rows.csv \
  --output first_contact_clean_pharmacotopology_layer_run/real_folding_10_dashboard.html
```

When the rows are real and locked, the benchmark can be run with:

```bash
python3 scripts/run_folding_topology_benchmark.py \
  --benchmark-file data/folding_benchmarks_real_500.locked.json \
  --require-external \
  --report-output first_contact_clean_pharmacotopology_layer_run/real_folding_500_report.json \
  --csv-output first_contact_clean_pharmacotopology_layer_run/real_folding_500_rows.csv
```

And visualized with:

```bash
python3 scripts/render_folding_benchmark_dashboard.py \
  --report first_contact_clean_pharmacotopology_layer_run/real_folding_500_report.json \
  --csv first_contact_clean_pharmacotopology_layer_run/real_folding_500_rows.csv \
  --output first_contact_clean_pharmacotopology_layer_run/real_folding_500_dashboard.html
```

## Output Fields

Each benchmark comparison writes:

```text
protein_id
sequence_length
source_database
source_accession
reference_structure_source
reference_label_source
predicted_topology_signature
reference_topology_signature
predicted_internal_fold_class
contact_map_similarity
fold_class_match
uncertainty_radius
evidence_readiness
failure_reason
is_external_reference
reference_topology_signature_kind
curation_notes
```

The default `failure_reason` is expected to say that an external structure
benchmark is not attached. That is a useful failure: it prevents placeholder
numbers from being mistaken for evidence.

External benchmark runs also write sidecar files beside the report:

```text
*_certificate.json
*_failures.csv
*_confusion_matrix.csv
```

When a benchmark file is loaded, the JSON report also includes:

```text
reference_dataset_validation
```

This records how many rows were loaded, how many had external reference-source
labels, whether external rows were required, and any violations or warnings.

The dashboard makes the following proof surfaces visible:

```text
confusion matrix
similarity distribution
per-class accuracy
worst mismatch table
radar overlay
locked benchmark certificate
dataset hash
commit hash
lock blockers
```

## Boundary

This module does not:

```text
solve protein folding
predict atomic structure
generate protein sequences
design molecules
recommend compounds
make clinical claims
```

It is a benchmark adapter and hypothesis-review shell. Its best use is to make
the next evidence step obvious.
