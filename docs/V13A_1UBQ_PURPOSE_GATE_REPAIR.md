# V13a 1UBQ purpose-gate readout

This is a postprocess-only repair for the V13a 1UBQ transfer run. It does not rerun OpenMM and does not change the locked V12/V13 role ontology. It reads the completed V13a trajectories and asks whether a single-domain compact target can be evaluated through a balanced-core purpose policy instead of a 4AKE/V5-style strict/rescue preflight.

Guardrails:

- no hardcoded `1UBQ -> threshold X`
- no hardcoded `0.62`
- no native precision used for threshold choice
- no selector rule tuning
- no partial/synthetic/native target fallback
- claim_allowed stays false

Run after `V13a_1UBQ_REPAIR_FIXED` has real trajectories:

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/src:$REPO_ROOT/scripts"

python3 -m compileall -q src scripts
PYTHONPATH=src:scripts python3 -m pytest -q

bash "$REPO_ROOT/scripts/run_v13a_1ubq_purpose_gate_readout.sh"
python3 "$REPO_ROOT/scripts/print_v13a_purpose_gate_summary.py"
```

Output artifacts:

- `first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_PURPOSE_GATE_READOUT/v13a_1ubq_purpose_gate_readout_certificate.json`
- `first_contact_clean_pharmacotopology_layer_run/V13a_1UBQ_PURPOSE_GATE_READOUT/v13a_1ubq_purpose_gate_frequency_sweep.csv`

Interpretation:

- `purpose_readout_core_found`: postprocess found an internally supported balanced compact core; this is a readout, not a biological claim. Review before any official purpose-aware runtime rerun.
- `clean_abstain_no_stable_frequency_band`: trajectories/provenance were present, but no internal frequency band passed the purpose-fit checks.

The old V5 preflight may remain blocked by `effective_strict_count = 0`; for `single_domain_compact`, this script reports that strict/rescue absence is not automatically fatal if balanced evidence exists and noise/pollution checks are clean.
