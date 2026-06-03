# Pharmacotopology

Standalone KNOT `Φ.review` pharmacotopology perturbation simulator.

This project visualizes hypothesis-only topology perturbations. It does not
recommend medication, infer clinical treatment, mention dosage, or provide
real-world prescribing guidance.

## Run

```bash
python3 scripts/run_clean_pharmacotopology_layer.py
python3 scripts/render_pharmacotopology_dashboard.py
```

Default run output:

```text
first_contact_clean_pharmacotopology_layer_run/
```

Generated local dashboard:

```text
first_contact_clean_pharmacotopology_layer_run/pharmacotopology_dashboard.html
```

The dashboard is static HTML/CSS/SVG and requires no web server or external
packages.

## Validate

```bash
python3 -m pytest
python3 scripts/validate_field_trace.py first_contact_clean_pharmacotopology_layer_run
python3 scripts/measure_field_trace.py first_contact_clean_pharmacotopology_layer_run
```

Expected metric posture:

```text
coherence = 1.0
stop_integrity = 1.0
surface_leakage = 0.0
```

## Doctrine

Visualize this as:

```text
perturbation changes topology dimensions
```

Do not visualize it as:

```text
pill works / pill does not work
```

Safety footer:

```text
Φ.review does not recommend medication.
Φ.review does not infer treatment.
Φ.review only simulates bounded topology perturbations.
Clinical claims are denied unless externally validated.
```
