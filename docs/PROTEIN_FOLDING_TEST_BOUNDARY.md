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

That result is frozen as:

```text
real_external_label_benchmark_v0
```

It is real in the sense that the rows contain real sequences and external
accessions. It is still label/prototype-based in the topology-reference channel:

```text
sequence -> predicted topology
external label -> locked broad-class prototype
compare
```

It is not yet a structure-derived folding test.

## Structure-Derived Topology Benchmark

The next benchmark channel is:

```text
structure_derived_topology_benchmark
```

It keeps three evidence channels separate:

```text
A. sequence-derived prediction
B. structure-derived topology evidence
C. external human/database label
```

The report therefore shows:

```text
prediction_vs_structure_score
prediction_vs_label_score
structure_vs_label_agreement
```

For folded PDB rows, the structure evidence is extracted from local PDB
coordinate files using a C-alpha contact graph plus HELIX/SHEET/B-factor
signals. For DisProt rows, the evidence channel is explicitly
`disorder_reference`, not coordinates.

The compact checked-in evidence file is:

```text
data/folding_benchmarks_real_10_structure_evidence.json
```

It records:

```text
8 coordinate_contact_graph rows
2 disorder_reference rows
folding_problem_solved = false
```

To rebuild that evidence file from local PDB files:

```bash
python3 scripts/extract_structure_topology_signatures.py \
  --benchmark-file data/folding_benchmarks_real_10.locked.json \
  --pdb-dir /path/to/pdb_files \
  --output data/folding_benchmarks_real_10_structure_evidence.json
```

To run the structure benchmark:

```bash
python3 scripts/run_structure_folding_topology_benchmark.py
```

The current checked-in structure report says:

```text
prediction_vs_structure_accuracy = 0.3
prediction_vs_label_accuracy = 0.4
structure_vs_label_agreement_rate = 0.9
sequence_order_sensitivity_score = 0.0
composition_only_warning = true
revision_required = true
claim_allowed = false
folding_problem_solved = false
```

That is a useful contradiction/revision signal, not a model success claim.

## Sequence-Order Controls

The structure benchmark also runs internal perturbation controls:

```text
real sequence
composition-preserving shuffled sequence
reversed sequence
local-window shuffled sequence
hydrophobic-cluster destroyed sequence
charge-pattern scrambled sequence
```

Generated control sequences are not written to the output files. The output
records only control type, whether composition was preserved, predicted class
change, and signature delta.

The most important metric is:

```text
sequence_order_sensitivity_score
```

The structure benchmark value is `0.0`, meaning the original recipe did not yet
respond to sequence order changes. That is the point of the test.

## Order-Aware Recipe Layer

The next layer is:

```text
order_aware_folding_topology_benchmark
```

It still receives sequence only. It does not receive CATH labels, PDB classes,
DisProt labels, structure-derived signatures, or reference fold classes during
prediction.

It extracts:

```text
hydrophobic cluster topology
charge pattern topology
proline/glycine breaker distribution
cysteine spacing / bridge potential
windowed local structure pressure
segment boundary contrast
long-range closure potential
predicted contact-prior graph
```

Run it with:

```bash
python3 scripts/run_order_aware_folding_topology_benchmark.py
```

It writes:

```text
real_folding_10_order_aware_report.json
real_folding_10_order_aware_rows.csv
real_folding_10_contact_prior.csv
real_folding_10_control_separation.csv
real_folding_10_order_aware_dashboard.html
real_folding_10_motif_alignment_report.json
real_folding_10_motif_alignment_rows.csv
real_folding_10_failure_diagnosis.csv
real_folding_10_evidence_conflicts.csv
real_folding_10_motif_alignment_dashboard.html
```

The current checked-in order-aware report says:

```text
sequence_order_sensitivity_score = 0.278079
real_vs_shuffled_separation_mean = 0.278079
contact_prior_signal_seen = true
recipe_order_blind = false
composition_only_warning = false
prediction_vs_structure_accuracy = 0.2
prediction_vs_label_accuracy = 0.1
revision_required = true
claim_allowed = false
folding_problem_solved = false
```

This means the recipe now reacts to sequence order and contact-prior topology,
but it still has not earned a folding claim. Scaling beyond the 10-row slice
should wait until the order-aware signal remains stable and the
prediction-vs-structure behavior improves without threshold tuning.

## Motif-to-Structure Alignment Layer

The next diagnostic layer is:

```text
motif_to_structure_alignment_benchmark
```

It still receives sequence only during prediction. It does not receive CATH
labels, DisProt labels, PDB classes, structure-derived signatures, or reference
fold classes until after the motif evidence vector has already been produced.

It makes the topology evidence vector the first output:

```text
alpha periodicity
beta alternation
compact core
disorder run
domain boundary
long-range closure
breaker / turn
charge frustration
```

Only after that evidence vector exists does the layer map it to a provisional
broad class, compare it against structure/label channels, and emit per-protein
failure diagnosis. The practical output is not higher accuracy yet; it is a
clearer explanation of why the current recipe fails.

Run it with:

```bash
python3 scripts/run_motif_alignment_benchmark.py
```

It writes:

```text
real_folding_10_motif_alignment_report.json
real_folding_10_motif_alignment_rows.csv
real_folding_10_failure_diagnosis.csv
real_folding_10_evidence_conflicts.csv
real_folding_10_motif_alignment_dashboard.html
```

The current checked-in motif alignment report says:

```text
prediction_vs_structure_accuracy = 0.1
prediction_vs_label_accuracy = 0.1
sequence_order_sensitivity_score = 0.278079
real_vs_shuffled_separation_mean = 0.278079
contact_prior_signal_seen = true
motif_signal_seen = true
evidence_conflict_mean = 0.895054
uncertainty_gating_used = true
forced_prediction_count = 2
abstained_prediction_count = 8
high_confidence_wrong_count = 0
revision_required = true
claim_allowed = false
folding_problem_solved = false
```

This means the layer has become more cautious rather than more accurate. That
is intentional for this revision. A high-conflict motif vector now usually
abstains, and high-confidence wrong predictions are gated to zero on this
locked 10-row slice.

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

Structure-derived benchmark runs write:

```text
real_folding_10_structure_report.json
real_folding_10_structure_rows.csv
real_folding_10_structure_dashboard.html
real_folding_10_order_controls.csv
real_folding_10_falsification_report.json
real_folding_10_order_aware_report.json
real_folding_10_order_aware_rows.csv
real_folding_10_contact_prior.csv
real_folding_10_control_separation.csv
real_folding_10_order_aware_dashboard.html
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
