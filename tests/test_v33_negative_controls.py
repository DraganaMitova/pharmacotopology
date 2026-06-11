from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v33_negative_controls_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v33_negative_controls", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_csv(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_stub_scripts(tmp_path: Path) -> None:
    scripts = tmp_path / "scripts"
    scripts.mkdir(parents=True, exist_ok=True)
    (scripts / "run_v33_constraint_backed_operator_readout_v0.py").write_text(
        "from pathlib import Path\n"
        "REPO_ROOT = Path('.')\n"
        "def build_v33(v32):\n"
        "    failed=[]\n"
        "    if v32.get('preflight_status') != 'V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED': failed.append('V32_status_not_ready_for_V33')\n"
        "    selected=v32.get('selected_V33_target')\n"
        "    if selected and selected != 'KcsA': failed.append('selected_target_not_supported_by_V33_v0')\n"
        "    rows=[]\n"
        "    for tr in v32.get('target_rows', []):\n"
        "        if tr.get('target') == 'KcsA': rows = tr.get('valid_real_external_constraint_rows', [])\n"
        "    types={r.get('evidence_type') for r in rows}\n"
        "    missing=[]\n"
        "    if 'pore_filter_coupling' not in types: missing.append('pore_filter_external_constraint_rows')\n"
        "    if 'assembly_interface_constraint' not in types: missing.append('assembly_interface_external_constraint_rows')\n"
        "    if failed: status='V33_BLOCKED_PRECONDITION_OR_PROVENANCE_FAILURE_CLAIM_DISABLED'\n"
        "    elif missing: status='V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED'\n"
        "    else: status='V33_CONSTRAINT_BACKED_OPERATOR_READOUT_PASSED_CLAIM_DISABLED'\n"
        "    return {'readout_status': status, 'operator_readout_found': status.endswith('PASSED_CLAIM_DISABLED'), 'missing_operator_evidence': missing, 'failed_checks': failed, 'claim_allowed': False, 'new_MD_allowed': False}\n",
        encoding="utf-8",
    )
    (scripts / "run_v32_external_constraint_source_import_preflight_v0.py").write_text(
        "from pathlib import Path\n"
        "REPO_ROOT = Path('.')\n"
        "def build_v32(v31, manifest):\n"
        "    rows=manifest.get('rows', [])\n"
        "    failed=[]\n"
        "    valid=[]\n"
        "    role=[]\n"
        "    for row in rows:\n"
        "        fp=row.get('file_path','')\n"
        "        et=row.get('evidence_type')\n"
        "        if fp.startswith('first_contact_clean_pharmacotopology_layer_run/'):\n"
        "            failed.append('internal_generated_artifact_supplied_as_external_source')\n"
        "        elif et in {'pore_filter_coupling','assembly_interface_constraint'}: valid.append(row)\n"
        "        elif et == 'annotation_context': role.append(row)\n"
        "    types={r.get('evidence_type') for r in valid}\n"
        "    ready=('pore_filter_coupling' in types and 'assembly_interface_constraint' in types and not failed)\n"
        "    if failed: status='V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_BLOCKED'\n"
        "    elif ready: status='V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED'\n"
        "    else: status='V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED'\n"
        "    return {'preflight_status': status, 'selected_V33_target': 'KcsA' if ready else None, 'provenance_clean': not failed, 'preflight_failed_checks': failed, 'claim_allowed': False, 'new_MD_allowed': False, 'target_rows': [{'target':'KcsA','valid_real_external_constraint_rows':valid,'role_context_rows':role}]}\n",
        encoding="utf-8",
    )


def _v31():
    return {
        "preflight_status": "V31_CLEAN_ABSTAIN_NO_REAL_EXTERNAL_CONSTRAINTS_FOR_SELECTED_TARGETS",
        "selected_V31_targets": ["XCL1_lymphotactin", "KcsA"],
        "claim_allowed": False,
        "new_md_executed": False,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "native_metrics_used_for_selection": False,
    }


def _v32(tmp_path: Path):
    pore = tmp_path / "data" / "external_constraints" / "KcsA" / "pore_filter" / "pore.csv"
    iface = tmp_path / "data" / "external_constraints" / "KcsA" / "assembly_interface" / "iface.csv"
    _write_csv(pore, "constraint_id,chain,residue_number,constraint_class\nP1,A,75,pore_filter_potassium_coordination_contact\n")
    _write_csv(iface, "constraint_id,chain_a,residue_i,chain_b,residue_j,constraint_class\nI1,A,10,B,10,assembly_interface_heavy_atom_contact\n")
    return {
        "preflight_status": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED",
        "selected_V33_target": "KcsA",
        "selected_V33_panel": "V33_CONSTRAINT_BACKED_OPERATOR_READOUT",
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "provenance_clean": True,
        "target_rows": [
            {"target": "KcsA", "valid_real_external_constraint_rows": [
                {"target": "KcsA", "file_path": "data/external_constraints/KcsA/pore_filter/pore.csv", "evidence_type": "pore_filter_coupling"},
                {"target": "KcsA", "file_path": "data/external_constraints/KcsA/assembly_interface/iface.csv", "evidence_type": "assembly_interface_constraint"},
            ]}
        ],
    }


def test_v33_negative_controls_pass_expected_blocks_and_abstains(tmp_path: Path) -> None:
    mod = _load_module()
    mod.REPO_ROOT = tmp_path
    mod.RUN_ROOT = tmp_path / "first_contact_clean_pharmacotopology_layer_run"
    _write_stub_scripts(tmp_path)
    v33 = {"readout_status": "V33_CONSTRAINT_BACKED_OPERATOR_READOUT_PASSED_CLAIM_DISABLED", "operator_readout_found": True}
    cert = mod.build_negative_controls(_v31(), _v32(tmp_path), v33)
    assert cert["negative_control_status"] == "V33_NEGATIVE_CONTROLS_PASSED_CLAIM_DISABLED"
    assert cert["passed_control_count"] == cert["control_count"]
    assert cert["claim_allowed"] is False
    assert cert["new_MD_allowed"] is False
    assert cert["positive_folding_evidence_found"] is False
    assert cert["folding_problem_solved"] is False
    by_id = {row["control_id"]: row for row in cert["controls"]}
    assert by_id["internal_runtime_source_poison_blocked_by_v32"]["passed"] is True
    assert by_id["annotation_only_context_does_not_select_v33"]["passed"] is True
    assert by_id["missing_pore_filter_bucket_abstains"]["passed"] is True
    assert by_id["missing_assembly_interface_bucket_abstains"]["passed"] is True


def test_v33_negative_control_writer_outputs_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    cert = {
        "negative_control_status": "V33_NEGATIVE_CONTROLS_PASSED_CLAIM_DISABLED",
        "control_count": 1,
        "passed_control_count": 1,
        "controls": [{"control_id": "dummy", "observed_status": "ok", "expected_status": "ok", "passed": True, "reason": "test"}],
        "claim_allowed": False,
        "new_MD_allowed": False,
        "folding_problem_solved": False,
        "next_action": "run_discriminative_or_randomized_constraint_controls_before_any_MD",
        "locked_interpretation": "test",
    }
    paths = mod.write_outputs(tmp_path / "out", cert)
    for path in paths.values():
        assert path.exists()
    written = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert written["next_decision"]["claim_allowed"] is False
