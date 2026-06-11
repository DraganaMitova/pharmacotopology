# V19 KcsA Schema/Guard Repair Commands

This patch does not change the science claim. It fixes a schema adapter issue in the V19 KcsA readout.

The V16 zero-MD certificate stores target rows under `target_rows`, while V19 v26 looked for `target_results` only. Because of that, V19 observed KcsA pore/filter, membrane, helix, and chain/interface evidence, but did not carry the soluble-core leakage guard forward and clean-abstained.

The repair reads both schemas and preserves the no-claim policy.

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"
```

Run tests:

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
  tests/test_v18b_p53_isolated_tad_abstain_context.py \
  tests/test_v19_kcsa_membrane_pore_evidence_readout.py
```

Run V19 again:

```bash
bash "$REPO_ROOT/scripts/run_v19_kcsa_membrane_pore_evidence_readout.sh"
python3 "$REPO_ROOT/scripts/print_v19_kcsa_membrane_pore_evidence_summary.py"
```

Expected status after repair:

```text
V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED
positive_pressure_evidence_found = True
membrane_pore_role_evidence_found = True
positive_folding_evidence_found = False
claim_allowed = False
new_md_executed = False
membrane_md_executed = False
soluble_core_misclassification_avoided = True
```
