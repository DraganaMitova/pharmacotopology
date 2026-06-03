# KNOT Pharmacotopology Layer

The pharmacotopology layer is a simulated research layer.

It treats a medication mechanism as a topology perturbation vector:

```text
mechanism vector
-> topology delta
-> pathology reduction score
-> collapse cost score
-> net topology health score
```

It does not map brand-name medications.

It does not make treatment recommendations.

It does not infer anything about a real patient.

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
python3 scripts/validate_field_trace.py first_contact_clean_pharmacotopology_layer_run
python3 scripts/measure_field_trace.py first_contact_clean_pharmacotopology_layer_run
```

The report is written to:

```text
first_contact_clean_pharmacotopology_layer_run/clean_pharmacotopology_layer_report.json
```

Expected safety posture:

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
