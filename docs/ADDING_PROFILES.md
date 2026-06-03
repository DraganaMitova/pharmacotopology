# Adding a Synthetic Topology Profile

Profiles in this repository are synthetic pressure maps for bounded hypothesis
comparison. They are not diagnoses, patient models, medication targets, or
clinical cohorts.

## Checklist

1. Add a `TopologyProfile` constant in `src/pharmacotopology/layer.py`.
2. Use a `profile_id` ending in `_topology_profile`.
3. Include every `PATHOLOGY_DIMENSIONS` key exactly once.
4. Keep each pressure value inside the simulator bounds, currently `0.0` to
   `1.25`.
5. Include this safety phrase in the description:

```text
not a diagnosis, patient model, or medication target
```

6. Register the profile in `DEFAULT_TOPOLOGY_PROFILES`.
7. Export the constant from `src/pharmacotopology/__init__.py` if it should be
   public.
8. Add or update a test that confirms the profile is available and bounded.

## Minimal Example

```python
DEFAULT_EXAMPLE_LIKE_PROFILE = TopologyProfile(
    profile_id="example_like_topology_profile",
    description=(
        "Synthetic example pressure profile only; not a diagnosis, patient "
        "model, or medication target."
    ),
    dimensions={
        "salience_amplification": 0.50,
        "recurrence_overbinding": 0.50,
        "symbolic_closure_pressure": 0.50,
        "threat_propagation": 0.50,
        "falsification_weakness": 0.50,
        "boundary_instability": 0.50,
        "agency_confusion": 0.50,
        "sensory_intrusion": 0.50,
        "cognitive_fragmentation": 0.50,
        "negative_shutdown": 0.50,
        "sleep_instability": 0.50,
    },
)
```

## Verify

```bash
python3 -m pytest
python3 scripts/run_clean_pharmacotopology_layer.py --profile example_like
python3 scripts/render_pharmacotopology_dashboard.py
python3 scripts/validate_field_trace.py first_contact_clean_pharmacotopology_layer_run
python3 scripts/measure_field_trace.py first_contact_clean_pharmacotopology_layer_run
```

The expected boundary stays unchanged:

```text
clinical_use_allowed = false
clinical_advice_created = false
medication_recommendation_created = false
real_patient_inference_created = false
```
