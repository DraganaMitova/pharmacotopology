from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v34_kcsa_discriminative_content_controls_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v34_controls", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _pore_rows():
    rows = []
    for chain in ["A", "B", "C", "D"]:
        for res in [75, 76, 77, 78]:
            rows.append({
                "chain": chain,
                "residue_number": str(res),
                "filter_motif": "TVGYG",
                "ion_name": "K",
                "constraint_class": "pore_filter_potassium_coordination_contact",
            })
    return rows


def _iface_rows():
    pairs = [("A", "B"), ("A", "C"), ("A", "D"), ("B", "C"), ("B", "D"), ("C", "D")]
    return [{"chain_a": a, "chain_b": b, "constraint_class": "assembly_interface_heavy_atom_contact"} for a, b in pairs]


def test_validate_kcsa_content_passes_clean_signature() -> None:
    mod = _load_module()
    out = mod.validate_kcsa_content(_pore_rows(), _iface_rows())
    assert out["valid"] is True
    assert out["failures"] == []
    assert out["summary"]["interface_chain_pair_count"] == 6


def test_validate_kcsa_content_rejects_adversarial_damage() -> None:
    mod = _load_module()
    pore = _pore_rows()
    for row in pore:
        row["ion_name"] = "NA"
        row["filter_motif"] = "FAKE"
    iface = _iface_rows()
    for row in iface:
        row["chain_a"] = "A"
        row["chain_b"] = "A"
    out = mod.validate_kcsa_content(pore, iface)
    assert out["valid"] is False
    assert "missing_potassium_ion_identity" in out["failures"]
    assert "missing_TVGYG_motif_label" in out["failures"]
    assert "insufficient_distinct_chain_pairs" in out["failures"]


def test_build_v34_runs_controls_from_v33_certificate(tmp_path: Path) -> None:
    mod = _load_module()
    mod.REPO_ROOT = tmp_path
    mod.RUN_ROOT = tmp_path / "first_contact_clean_pharmacotopology_layer_run"
    pore_path = tmp_path / "data" / "external_constraints" / "KcsA" / "pore_filter" / "pore.csv"
    iface_path = tmp_path / "data" / "external_constraints" / "KcsA" / "assembly_interface" / "iface.csv"
    pore_path.parent.mkdir(parents=True, exist_ok=True)
    iface_path.parent.mkdir(parents=True, exist_ok=True)
    pore_path.write_text("chain,residue_number,filter_motif,ion_name,constraint_class\n" + "\n".join(
        f"{r['chain']},{r['residue_number']},{r['filter_motif']},{r['ion_name']},{r['constraint_class']}" for r in _pore_rows()
    ) + "\n", encoding="utf-8")
    iface_path.write_text("chain_a,chain_b,constraint_class\n" + "\n".join(
        f"{r['chain_a']},{r['chain_b']},{r['constraint_class']}" for r in _iface_rows()
    ) + "\n", encoding="utf-8")
    v33 = {
        "readout_status": "V33_CONSTRAINT_BACKED_OPERATOR_READOUT_PASSED_CLAIM_DISABLED",
        "selected_V33_target": "KcsA",
        "claim_allowed": False,
        "new_MD_allowed": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
        "selected_target_import_rows": [
            {"evidence_type": "pore_filter_coupling", "file_path": "data/external_constraints/KcsA/pore_filter/pore.csv"},
            {"evidence_type": "assembly_interface_constraint", "file_path": "data/external_constraints/KcsA/assembly_interface/iface.csv"},
        ],
    }
    cert = mod.build_v34(v33)
    assert cert["control_status"] == "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_PASSED_CLAIM_DISABLED"
    assert cert["passed_control_count"] == cert["control_count"]
    assert cert["claim_allowed"] is False
    assert cert["new_MD_allowed"] is False
    assert cert["folding_problem_solved"] is False


def test_v34_writer_outputs_artifacts(tmp_path: Path) -> None:
    mod = _load_module()
    cert = {
        "control_status": "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_PASSED_CLAIM_DISABLED",
        "control_count": 1,
        "passed_control_count": 1,
        "controls": [{"control_id": "dummy", "passed": True, "observed_valid": True, "expected_valid": True, "observed_failures": [], "expected_failure_any": [], "reason": "test"}],
        "baseline_content_signature": {},
        "claim_allowed": False,
        "new_MD_allowed": False,
        "folding_problem_solved": False,
        "next_action": "move_to_non_coordinate_external_evolutionary_or_blind_holdout_tests",
        "locked_interpretation": "test",
    }
    paths = mod.write_outputs(tmp_path / "out", cert)
    for path in paths.values():
        assert path.exists()
    written = json.loads(paths["certificate"].read_text(encoding="utf-8"))
    assert written["next_decision"]["claim_allowed"] is False
