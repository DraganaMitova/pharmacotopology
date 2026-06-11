from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v33_constraint_backed_operator_readout_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v33_readout", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_csv(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _v32_ready(tmp_path: Path) -> dict:
    pore = tmp_path / "data" / "external_constraints" / "KcsA" / "pore_filter" / "kcsa_pore.csv"
    iface = tmp_path / "data" / "external_constraints" / "KcsA" / "assembly_interface" / "kcsa_iface.csv"
    _write_csv(
        pore,
        "constraint_id,chain,residue_number,residue_name,filter_motif,ion_name,min_distance_angstrom,constraint_class\n"
        "P1,A,75,THR,TVGYG,K,2.800,pore_filter_potassium_coordination_contact\n",
    )
    _write_csv(
        iface,
        "constraint_id,chain_a,residue_i,chain_b,residue_j,min_distance_angstrom,constraint_class\n"
        "I1,A,10,B,10,3.700,assembly_interface_heavy_atom_contact\n",
    )
    return {
        "preflight_status": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED",
        "selected_V33_target": "KcsA",
        "selected_V33_panel": "V33_CONSTRAINT_BACKED_OPERATOR_READOUT",
        "claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "provenance_clean": True,
        "target_rows": [
            {
                "target": "KcsA",
                "valid_real_external_constraint_rows": [
                    {
                        "target": "KcsA",
                        "file_path": "data/external_constraints/KcsA/pore_filter/kcsa_pore.csv",
                        "evidence_type": "pore_filter_coupling",
                        "state_or_context": "TVGYG_selectivity_filter_or_pore_filter_external_coordinate_contact",
                        "counts_as_real_external_constraint": True,
                        "valid_import_row": True,
                    },
                    {
                        "target": "KcsA",
                        "file_path": "data/external_constraints/KcsA/assembly_interface/kcsa_iface.csv",
                        "evidence_type": "assembly_interface_constraint",
                        "state_or_context": "tetramer_chain_interface_context_external_coordinate_contact",
                        "counts_as_real_external_constraint": True,
                        "valid_import_row": True,
                    },
                ],
            }
        ],
    }


def test_v33_passes_kcsa_constraint_backed_operator_readout(tmp_path: Path) -> None:
    mod = _load_module()
    mod.REPO_ROOT = tmp_path
    cert = mod.build_v33(_v32_ready(tmp_path))
    assert cert["readout_status"] == "V33_CONSTRAINT_BACKED_OPERATOR_READOUT_PASSED_CLAIM_DISABLED"
    assert cert["selected_V33_target"] == "KcsA"
    assert cert["operator_readout_found"] is True
    assert cert["positive_folding_evidence_found"] is False
    assert cert["folding_problem_solved"] is False
    assert cert["claim_allowed"] is False
    assert cert["new_MD_allowed"] is False
    readout = cert["constraint_backed_operator_readout"]
    assert readout["pore_filter_row_count"] == 1
    assert readout["assembly_interface_row_count"] == 1
    assert "pore_filter_operator" in readout["operator_buckets_assigned"]
    assert "assembly_interface_operator" in readout["operator_buckets_assigned"]


def test_v33_blocks_when_v32_not_ready(tmp_path: Path) -> None:
    mod = _load_module()
    mod.REPO_ROOT = tmp_path
    v32 = _v32_ready(tmp_path)
    v32["preflight_status"] = "V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED"
    v32["selected_V33_target"] = None
    cert = mod.build_v33(v32)
    assert cert["readout_status"] == "V33_BLOCKED_PRECONDITION_OR_PROVENANCE_FAILURE_CLAIM_DISABLED"
    assert "V32_status_not_ready_for_V33" in cert["failed_checks"]
    assert cert["claim_allowed"] is False


def test_v33_abstains_if_interface_bucket_missing(tmp_path: Path) -> None:
    mod = _load_module()
    mod.REPO_ROOT = tmp_path
    v32 = _v32_ready(tmp_path)
    v32["target_rows"][0]["valid_real_external_constraint_rows"] = [v32["target_rows"][0]["valid_real_external_constraint_rows"][0]]
    cert = mod.build_v33(v32)
    assert cert["readout_status"] == "V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED"
    assert "assembly_interface_external_constraint_rows" in cert["missing_operator_evidence"]
    assert cert["positive_folding_evidence_found"] is False


def test_v33_writer_outputs_certificate_report_and_decision(tmp_path: Path) -> None:
    mod = _load_module()
    mod.REPO_ROOT = tmp_path
    cert = mod.build_v33(_v32_ready(tmp_path))
    paths = mod.write_outputs(tmp_path / "out", cert)
    for key in ["certificate", "readout", "report", "decision"]:
        assert paths[key].exists(), key
    written = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert written["claim_allowed"] is False
    assert written["folding_problem_solved"] is False
    assert written["next_decision"]["next_action"] == "run_V33_negative_controls_before_any_claim_or_MD"
