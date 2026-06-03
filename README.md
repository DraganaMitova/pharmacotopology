# Pharmacotopology

A bounded hypothesis workbench and visual-only research simulator for exploring
how hypothesis-level mechanism vectors might perturb an abstract topology
profile.

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

The practical use is not clinical decision support. It is structured hypothesis
work: make assumptions explicit, compare them under the same scoring rules, add
evidence weights when evidence exists, and keep uncertainty visible until a
claim can be externally calibrated.

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

The repository now exposes that gap directly through calibration readiness
metadata: evidence weight, uncertainty radius, score intervals, blockers, and
next calibration steps. That makes the project useful as research
infrastructure without pretending to be useful as medicine.

## Repository Shape

```text
src/pharmacotopology/                 core simulator and boundary model
scripts/run_clean_pharmacotopology_layer.py
scripts/render_pharmacotopology_dashboard.py
scripts/render_profile_comparison_dashboard.py
scripts/export_pharmacotopology_csv.py
scripts/run_sensitivity_analysis.py
scripts/explore_sensitivity.py
scripts/validate_field_trace.py
scripts/measure_field_trace.py
tests/                               pytest coverage for ranking and safety
docs/ADDING_PROFILES.md              small guide for adding synthetic profiles
PHARMACOTOPOLOGY_LAYER.md             longer design note
LICENSE                              MIT license
```

The package has no runtime dependencies and targets Python 3.9+.

## Run

Generate a clean run:

```bash
python3 scripts/run_clean_pharmacotopology_layer.py
```

Run against another synthetic topology profile:

```bash
python3 scripts/run_clean_pharmacotopology_layer.py --profile anxiety_like
```

Render the static dashboard:

```bash
python3 scripts/render_pharmacotopology_dashboard.py
```

Render the multi-profile dashboard:

```bash
python3 scripts/render_profile_comparison_dashboard.py
```

Default artifacts are written to:

```text
first_contact_clean_pharmacotopology_layer_run/
```

The dashboard is generated at:

```text
first_contact_clean_pharmacotopology_layer_run/pharmacotopology_dashboard.html
```

The multi-profile dashboard is generated at:

```text
first_contact_clean_pharmacotopology_layer_run/multi_profile_dashboard.html
```

It is plain HTML/CSS/SVG and does not require a web server or external packages.

## Profiles

The built-in source profiles are synthetic pressure maps:

```text
schizophrenia_like
depression_like
mania_like
anxiety_like
mixed_state_like
```

They are for hypothesis comparison only. They are not diagnoses, patient
models, or medication targets.

To add another bounded synthetic profile, see [docs/ADDING_PROFILES.md](docs/ADDING_PROFILES.md).

## Export CSV

Export machine-readable rankings and per-dimension deltas:

```bash
python3 scripts/export_pharmacotopology_csv.py
```

Default CSV outputs:

```text
first_contact_clean_pharmacotopology_layer_run/pharmacotopology_rankings.csv
first_contact_clean_pharmacotopology_layer_run/pharmacotopology_deltas.csv
```

The ranking CSV includes score intervals, evidence weight, uncertainty radius,
and evidence readiness labels. The delta CSV includes every topology dimension
plus collapse cost.

## Calibration Readiness Table

The dashboard includes a per-vector calibration table with:

```text
mechanism_id
evidence_stage
evidence_weight
uncertainty_radius
evidence_readiness_label
primary evidence sources
evidence refs
calibration blockers
interval kind
pathology_reduction interval
collapse_cost interval
net_topology_health interval
```

The current default vectors remain `uncalibrated_hypothesis`, which is the
honest starting point. The table is meant to show exactly where evidence would
attach later.

The dashboard also includes a synthetic profile comparison table that runs the
same mechanism set against every built-in profile and reports top mechanism,
fit label, net interval, and destabilizing count. It is a ranking robustness
view, not a clinical comparison.

## Multi-Profile Dashboard

Render a dedicated cross-profile dashboard:

```bash
python3 scripts/render_profile_comparison_dashboard.py
```

It includes:

```text
baseline topology maps for every built-in profile
combined ranking table with rank / net score per profile
mechanism/profile net-score heatmap
consistently positive vs highly profile-specific summary
```

## Sensitivity Analysis

Run a local pressure sweep to see whether rankings are stable when source
profile pressures move up or down:

```bash
python3 scripts/run_sensitivity_analysis.py
```

Sweep a smaller set of dimensions:

```bash
python3 scripts/run_sensitivity_analysis.py --profile anxiety_like --dimensions threat_propagation,sleep_instability
```

Default sensitivity outputs:

```text
first_contact_clean_pharmacotopology_layer_run/sensitivity_analysis_report.json
first_contact_clean_pharmacotopology_layer_run/sensitivity_rankings.csv
```

This is a robustness review, not a clinical interpretation.

## Sensitivity Explorer

Explore a mechanism assumption with deterministic Monte Carlo perturbations:

```bash
python3 scripts/explore_sensitivity.py --profile schizophrenia_like --mechanism nmda_support_like --vary collapse_cost 0.05:0.25
```

You can also vary a delta dimension:

```bash
python3 scripts/explore_sensitivity.py --mechanism nmda_support_like --vary delta:cognitive_fragmentation -0.30:-0.05
```

Default explorer outputs:

```text
first_contact_clean_pharmacotopology_layer_run/sensitivity_explorer_report.json
first_contact_clean_pharmacotopology_layer_run/sensitivity_explorer_samples.csv
```

The explorer reports a net-score distribution for the selected mechanism and a
robustness score for every mechanism under sampled perturbations.

## Interactive Wrapper

A small Streamlit or Gradio wrapper would be a natural later layer for local
exploration: profile selection, pressure sliders, ranking preview, CSV export,
and the same safety warnings. It is intentionally not included in the core
package yet, so the current project remains dependency-free and easy to audit.

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
calibration_status = uncalibrated_hypothesis_workbench
practical_use = bounded_hypothesis_comparison_and_falsification
clinical_use_allowed = false
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
calibration_readiness_report.json
pharmacotopology_rankings.csv
pharmacotopology_deltas.csv
sensitivity_analysis_report.json
sensitivity_rankings.csv
pharmacotopology_dashboard.html
multi_profile_dashboard.html
sensitivity_explorer_report.json
sensitivity_explorer_samples.csv
field_validation.json
field_metrics.json
```

The JSON files are meant to make the review trace auditable. The dashboard is
meant to make the topology ranking, evidence readiness, and uncertainty
intervals easier to inspect without changing the safety boundary.

## Safety Footer

```text
Φ.review does not recommend medication.
Φ.review does not infer treatment.
Φ.review only simulates bounded topology perturbations.
Clinical claims are denied unless externally validated.
```

## Status

This is an early research prototype with a practical non-clinical target:
bounded hypothesis comparison and falsification. The most valuable
contributions right now are careful critique of the modeling assumptions,
evidence-source mapping, calibration datasets, clearer validation boundaries,
and better ways to explain what the scores do and do not mean.

## License

MIT. See [LICENSE](LICENSE).
