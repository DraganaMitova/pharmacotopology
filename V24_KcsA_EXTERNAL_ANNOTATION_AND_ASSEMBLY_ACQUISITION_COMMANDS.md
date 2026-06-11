# V24 KcsA External Annotation and Assembly Acquisition

This is a zero-MD acquisition/readout sprint for KcsA after V19. It locks external membrane/pore/assembly-context annotations, preserves the pore/filter and TM scaffold evidence, scans coupling/MSA availability, and decides the next panel.

It does **not** run membrane MD, does **not** claim a tetramer/whole-channel fold, does **not** use fixed residue cutoffs, and keeps `claim_allowed=false`.

## Run

```bash
export REPO_ROOT="/Users/draganamitova/My Projects/pharmacotopology"
cd "$REPO_ROOT" || exit 1
export PYTHONPATH="$REPO_ROOT/scripts:$REPO_ROOT/src"

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
  tests/test_v19_kcsa_membrane_pore_evidence_readout.py \
  tests/test_v19_kcsa_membrane_pore_evidence_lock.py \
  tests/test_v20_xcl1_state_specific_evidence_readout.py \
  tests/test_v21_pressure_evidence_panel_summary.py \
  tests/test_v22_external_evidence_and_annotation_acquisition_panel.py \
  tests/test_v23_p53_tad_mdm2_external_evidence_contrast_test.py \
  tests/test_v24_kcsa_external_annotation_and_assembly_acquisition.py

bash "$REPO_ROOT/scripts/run_v24_kcsa_external_annotation_and_assembly_acquisition.sh"
python3 "$REPO_ROOT/scripts/print_v24_kcsa_external_annotation_and_assembly_acquisition_summary.py"
```

Expected status:

```text
V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION_PASSED_CLAIM_DISABLED
```

## Requested artifacts

```text
v24_kcsa_external_annotation_manifest.json
v24_kcsa_assembly_preflight.json
v24_kcsa_pore_filter_annotation.json
v24_kcsa_membrane_topology_annotation.json
v24_kcsa_external_coupling_availability.json
v24_kcsa_next_decision.json
```
