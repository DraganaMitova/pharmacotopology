# KNOT Pharmacotopology Layer

The pharmacotopology layer is a simulated research layer and bounded
hypothesis workbench.

It treats a protein-centered mechanism as a topology perturbation vector:

```text
mechanism vector
-> topology delta
-> pathology reduction score
-> collapse cost score
-> net topology health score
```

The explicit reframe is:

```text
abstract_compound_class
-> protein_family
-> protein_state_shift
-> pathway/network perturbation
-> topology_delta
-> collapse_cost
-> evidence_readiness
```

It does not map brand-name medications.

It does not make treatment recommendations.

It does not infer anything about a real patient.

It does not propose molecules, optimize compounds, or make drug-design claims.

Its practical use is assumption inspection: compare mechanism hypotheses under
the same topology rules, keep collapse cost visible, expose uncertainty, and
record what evidence would be needed before any external claim could be made.

## First Prototype

The first admitted profile is a schizophrenia-like topology profile, represented
only as simulated pressure values:

```text
salience amplification
recurrence overbinding
symbolic closure pressure
threat propagation
falsification weakness
boundary instability
agency confusion
sensory intrusion
cognitive fragmentation
negative shutdown
sleep instability
```

Additional synthetic source profiles are available for bounded comparison:

```text
depression-like topology profile
mania-like topology profile
anxiety-like topology profile
mixed-state-like topology profile
```

All profiles are synthetic pressure maps. They are not diagnoses, patient
models, or medication targets.

Lower pressure means closer to bounded pattern review.

The first mechanism vectors are hypothesis-only receptor/mechanism classes:

```text
D2 antagonist-like
D2 partial agonist-like
5-HT2A antagonist-like
NMDA support-like
GABA stabilizing-like
muscarinic modulation-like
adrenergic dampening-like
histamine sedation load-like
anticholinergic burden-like
glutamate amplification stressor-like
```

Each vector changes topology dimensions and carries a separate collapse cost.
The pathology reduction score is computed over the dimensions affected by that
mechanism vector; the full resulting topology state is still written for review.
Each result also carries an evidence readiness label, evidence weight,
uncertainty radius, and score interval.
The rendered dashboard includes a per-vector calibration readiness table so the
protein mechanism labels, evidence posture, primary source slots, blockers, and
model uncertainty intervals are inspectable beside the ranking.
The dashboard also includes a synthetic profile comparison view that applies
the same mechanism vectors to every built-in profile.

The layer therefore can distinguish:

```text
reduced unstable topology with low collapse cost
reduced unstable topology with high collapse cost
weak topology correction
destabilizing perturbation
```

## Run

```bash
python3 scripts/run_clean_pharmacotopology_layer.py
python3 scripts/run_clean_pharmacotopology_layer.py --profile anxiety_like
python3 scripts/render_profile_comparison_dashboard.py
python3 scripts/export_pharmacotopology_csv.py
python3 scripts/run_sensitivity_analysis.py --profile anxiety_like --dimensions threat_propagation,sleep_instability
python3 scripts/explore_sensitivity.py --profile schizophrenia_like --mechanism nmda_support_like --vary collapse_cost 0.05:0.25
python3 scripts/validate_field_trace.py first_contact_clean_pharmacotopology_layer_run
python3 scripts/measure_field_trace.py first_contact_clean_pharmacotopology_layer_run
```

The report is written to:

```text
first_contact_clean_pharmacotopology_layer_run/clean_pharmacotopology_layer_report.json
```

Calibration readiness is written to:

```text
first_contact_clean_pharmacotopology_layer_run/calibration_readiness_report.json
```

CSV exports are written to:

```text
first_contact_clean_pharmacotopology_layer_run/pharmacotopology_rankings.csv
first_contact_clean_pharmacotopology_layer_run/pharmacotopology_deltas.csv
```

Sensitivity analysis is written to:

```text
first_contact_clean_pharmacotopology_layer_run/sensitivity_analysis_report.json
first_contact_clean_pharmacotopology_layer_run/sensitivity_rankings.csv
first_contact_clean_pharmacotopology_layer_run/sensitivity_explorer_report.json
first_contact_clean_pharmacotopology_layer_run/sensitivity_explorer_samples.csv
```

The multi-profile dashboard is written to:

```text
first_contact_clean_pharmacotopology_layer_run/multi_profile_dashboard.html
```

Expected safety posture:

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
stop_integrity = 1.0
surface_leakage = 0.0
```

## Calibration Boundary

The initial numbers are hypotheses, not truth.

Later calibration may compare the simulated deltas against:

```text
clinical trial outcomes
side-effect profiles
receptor binding data
patient symptom scales
relapse rates
negative symptom data
cognition data
real-world discontinuation rates
```

Until that calibration exists, the layer remains a bounded perturbation model,
not clinical guidance.
