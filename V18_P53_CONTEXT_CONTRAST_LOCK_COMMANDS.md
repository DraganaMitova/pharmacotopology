# V18 p53 Context Contrast Lock Commands

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"
```

## Tests

```bash
python3 -m compileall -q src scripts

python3 -m pytest -q \
  tests/test_v13a_purpose_gate_readout.py \
  tests/test_v13b_hierarchical_purpose_topology_readout.py \
  tests/test_v14_unified_protein_esperanto_grammar_readout.py \
  tests/test_v15_dynamic_separation_grammar_readout.py \
  tests/test_v15_4ake_dynamic_grammar_bridge.py \
  tests/test_v15_dynamic_role_grammar_panel_lock.py \
  tests/test_v15_4ake_balanced_candidate_readout.py \
  tests/test_v16_target_manifest_and_role_expectation_lock.py \
  tests/test_v16_transfer_data_preflight.py \
  tests/test_v16_zero_md_role_transfer_readout.py \
  tests/test_v16_pressure_evidence_gap_lock.py \
  tests/test_v17_pressure_evidence_sprint_lock.py \
  tests/test_v18_p53_tad_mdm2_partner_induced_evidence_test.py \
  tests/test_v18_p53_partner_induced_evidence_lock.py \
  tests/test_v18b_p53_isolated_tad_abstain_context.py
```

## Run V18 evidence test and lock

```bash
bash "$REPO_ROOT/scripts/download_v16_pressure_target_structures.sh"

bash "$REPO_ROOT/scripts/run_v17_pressure_evidence_sprint_lock.sh"
python3 "$REPO_ROOT/scripts/print_v17_pressure_evidence_sprint_lock_summary.py"

bash "$REPO_ROOT/scripts/run_v18_p53_tad_mdm2_partner_induced_evidence_test.sh"
python3 "$REPO_ROOT/scripts/print_v18_p53_tad_mdm2_partner_induced_evidence_summary.py"

bash "$REPO_ROOT/scripts/run_v18_p53_partner_induced_evidence_lock.sh"
python3 "$REPO_ROOT/scripts/print_v18_p53_partner_induced_evidence_lock_summary.py"

bash "$REPO_ROOT/scripts/run_v18b_p53_isolated_tad_abstain_context.sh"
python3 "$REPO_ROOT/scripts/print_v18b_p53_isolated_tad_abstain_context_summary.py"
```

Expected final status:

```text
V18_P53_PARTNER_INDUCED_EVIDENCE_LOCKED
V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT_PASSED_CLAIM_DISABLED
claim_allowed = false
positive_folding_evidence_found = false
new_md_executed = false
```
