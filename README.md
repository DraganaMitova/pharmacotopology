# Pharmacotopology

A small, visual-only research simulator for exploring how hypothesis-level
mechanism vectors might perturb an abstract topology profile.

This repository is public because the idea is easier to examine when the model,
assumptions, generated artifacts, and safety boundaries are visible. It is not
a clinical tool. It does not recommend medication, infer treatment, map brand
names, mention dosage, or provide real-world prescribing guidance.

## What This Is

Pharmacotopology treats a mechanism class as a bounded perturbation over an
abstract topology state:

```text
mechanism vector
-> topology delta
-> pathology reduction score
-> collapse cost score
-> net topology health score
```

The included prototype uses a synthetic schizophrenia-like topology profile.
That profile is a simulated pressure map only; it is not a diagnosis, a patient
model, or a medication target.

The current layer compares hypothesis-only mechanism classes such as
`d2_antagonist_like`, `nmda_support_like`, and
`glutamate_amplification_stressor_like`. Each vector changes one or more
topology dimensions and carries a separate collapse cost. The ranking asks a
narrow modeling question:

```text
Did this perturbation move the simulated topology closer to bounded review,
and at what modeled collapse cost?
```

## What This Is Not

This project intentionally refuses clinical interpretation:

```text
not medication advice
not treatment ranking
not symptom prediction
not side-effect prediction
not brand-name mapping
not patient-specific inference
```

The useful reading is:

```text
perturbation changes topology dimensions
```

The unsafe reading is:

```text
pill works / pill does not work
```

## Why It Exists

The experiment is modest: make the hidden assumptions of a topology-style
pharmacology metaphor inspectable. A mechanism vector should not get to look
"helpful" just because it suppresses pressure somewhere; it also has to carry
its modeled collapse cost beside the reduction score.

That means the layer can distinguish between:

```text
reduced unstable topology with low collapse cost
reduced unstable topology with high collapse cost
weak topology correction
destabilizing perturbation
```

The numbers are hypotheses, not truth. Future calibration would need external
evidence such as trial outcomes, receptor binding data, side-effect profiles,
symptom scales, relapse data, cognition data, and discontinuation rates. Until
then, this remains a bounded perturbation model.

## Repository Shape

```text
src/pharmacotopology/                 core simulator and boundary model
scripts/run_clean_pharmacotopology_layer.py
scripts/render_pharmacotopology_dashboard.py
scripts/validate_field_trace.py
scripts/measure_field_trace.py
tests/                               pytest coverage for ranking and safety
PHARMACOTOPOLOGY_LAYER.md             longer design note
```

The package has no runtime dependencies and targets Python 3.9+.

## Run

Generate a clean run:

```bash
python3 scripts/run_clean_pharmacotopology_layer.py
```

Render the static dashboard:

```bash
python3 scripts/render_pharmacotopology_dashboard.py
```

Default artifacts are written to:

```text
first_contact_clean_pharmacotopology_layer_run/
```

The dashboard is generated at:

```text
first_contact_clean_pharmacotopology_layer_run/pharmacotopology_dashboard.html
```

It is plain HTML/CSS/SVG and does not require a web server or external packages.

## Validate

Run the tests:

```bash
python3 -m pytest
```

Validate and measure a generated trace:

```bash
python3 scripts/validate_field_trace.py first_contact_clean_pharmacotopology_layer_run
python3 scripts/measure_field_trace.py first_contact_clean_pharmacotopology_layer_run
```

Expected safety posture:

```text
coherence = 1.0
stop_integrity = 1.0
surface_leakage = 0.0
```

The run report also asserts:

```text
simulation_only = true
hypothesis_numbers_only = true
clinical_advice_created = false
medication_recommendation_created = false
real_patient_inference_created = false
brand_name_mapping_created = false
treatment_claim_created = false
voice_opened = false
schedule_created = false
command_queue_created = false
autonomous_loop_created = false
```

## Output Artifacts

A default run creates a small set of local artifacts:

```text
input_stream.jsonl
memory.jsonl
audit.jsonl
session_records.jsonl
clean_pharmacotopology_layer_report.json
pharmacotopology_dashboard.html
field_validation.json
field_metrics.json
```

The JSON files are meant to make the review trace auditable. The dashboard is
meant to make the topology ranking easier to inspect without changing the
safety boundary.

## Safety Footer

```text
Φ.review does not recommend medication.
Φ.review does not infer treatment.
Φ.review only simulates bounded topology perturbations.
Clinical claims are denied unless externally validated.
```

## Status

This is an early research prototype. The most valuable contributions right now
are careful critique of the modeling assumptions, clearer validation boundaries,
and better ways to explain what the scores do and do not mean.

## License

No license file is included yet. Until one is added, treat reuse rights as not
granted by default.
