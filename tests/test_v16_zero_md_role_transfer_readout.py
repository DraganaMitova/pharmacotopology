from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PREFLIGHT_SCRIPT = REPO_ROOT / "scripts" / "run_v16_transfer_data_preflight_v0.py"
READOUT_SCRIPT = REPO_ROOT / "scripts" / "run_v16_zero_md_role_transfer_readout_v0.py"
ROLE_MANIFEST = REPO_ROOT / "data" / "v16_locked_grammar_transfer_target_manifest.json"
MATERIAL_MANIFEST = REPO_ROOT / "data" / "v16_transfer_data_material_manifest.json"


def _v16_lock() -> dict:
    return {
        "kind": "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCK_v0",
        "lock_status": "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCKED",
        "source_v15_lock_status": "V15_THREE_PROTEIN_DYNAMIC_GRAMMAR_PANEL_LOCKED",
        "claim_allowed": False,
    }


def _write_mock_pdb(path: Path, chain_counts: dict[str, int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    atom_id = 1
    for chain, count in chain_counts.items():
        for resi in range(1, count + 1):
            lines.append(
                f"ATOM  {atom_id:5d}  CA  ALA {chain:1s}{resi:4d}    "
                f"{float(resi):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C\n"
            )
            atom_id += 1
    path.write_text("".join(lines) + "END\n", encoding="utf-8")


def _write_prov(path: Path, target_id: str, pdb_id: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "kind": "V16_PUBLIC_RCSB_STRUCTURE_PROVENANCE_v0",
                "target_id": target_id,
                "pdb_id": pdb_id,
                "source": "RCSB PDB public download",
                "usage_boundary": "target_context_only_not_native_selection_not_folding_solved_claim",
                "claim_allowed": False,
            }
        ),
        encoding="utf-8",
    )


def _ready_preflight_cert(tmp_path: Path) -> Path:
    material = json.loads(MATERIAL_MANIFEST.read_text())
    for target in material["targets"]:
        for mat in target["required_material"]:
            mat_path = tmp_path / mat["path"]
            prov_path = tmp_path / mat["provenance_path"]
            mat["path"] = str(mat_path)
            mat["provenance_path"] = str(prov_path)
            if target["target_id"] == "p53_TAD_MDM2":
                _write_mock_pdb(mat_path, {"A": 85, "B": 13})
            elif target["target_id"] == "KcsA":
                # V16b repair: 1K4C-style three-chain coordinate context is enough for
                # zero-MD role context readout, but not enough for tetramer/folding claims.
                _write_mock_pdb(mat_path, {"A": 219, "B": 212, "C": 103})
            else:
                _write_mock_pdb(mat_path, {"A": 80})
            _write_prov(prov_path, target["target_id"], mat["pdb_id"])
    material_path = tmp_path / "material.json"
    material_path.write_text(json.dumps(material), encoding="utf-8")
    lock_path = tmp_path / "v16_lock.json"
    lock_path.write_text(json.dumps(_v16_lock()), encoding="utf-8")
    out_dir = tmp_path / "preflight"
    subprocess.run(
        [
            sys.executable,
            str(PREFLIGHT_SCRIPT),
            "--role-manifest",
            str(ROLE_MANIFEST),
            "--material-manifest",
            str(material_path),
            "--v16-lock-cert",
            str(lock_path),
            "--out-dir",
            str(out_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    cert_path = out_dir / "v16_transfer_data_preflight_certificate.json"
    cert = json.loads(cert_path.read_text())
    assert cert["data_preflight_status"] == "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT"
    assert set(cert["ready_targets"]) == {"p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}
    return cert_path


def test_v16_kcsa_material_manifest_uses_three_chain_context_floor() -> None:
    manifest = json.loads(MATERIAL_MANIFEST.read_text())
    kcsa = next(t for t in manifest["targets"] if t["target_id"] == "KcsA")
    mat = kcsa["required_material"][0]
    assert mat["min_chains"] == 3
    assert "not a tetramer proof" in mat["material_sanity_note"]
    assert manifest["claim_allowed"] is False
    assert manifest["new_md_allowed"] is False
    assert manifest["target_specific_threshold_tuning_allowed"] is False


def test_v16_zero_md_role_transfer_readout_completed_with_ready_preflight(tmp_path: Path) -> None:
    preflight_cert = _ready_preflight_cert(tmp_path)
    out_dir = tmp_path / "readout"
    subprocess.run(
        [
            sys.executable,
            str(READOUT_SCRIPT),
            "--role-manifest",
            str(ROLE_MANIFEST),
            "--preflight-cert",
            str(preflight_cert),
            "--out-dir",
            str(out_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    cert = json.loads((out_dir / "v16_zero_md_role_transfer_readout_certificate.json").read_text())
    assert cert["role_transfer_status"] == "V16_ZERO_MD_ROLE_TRANSFER_READOUT_COMPLETED_CLAIM_DISABLED"
    assert set(cert["positive_role_targets"]) == {"p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["fixed_residue_cutoff_used"] is False
    assert cert["native_metrics_used_for_selection"] is False
    assert cert["forbidden_misclassification_violations"] == {}


def test_v16_zero_md_kcsa_does_not_claim_soluble_or_tetramer_fold(tmp_path: Path) -> None:
    preflight_cert = _ready_preflight_cert(tmp_path)
    out_dir = tmp_path / "readout"
    subprocess.run(
        [sys.executable, str(READOUT_SCRIPT), "--preflight-cert", str(preflight_cert), "--out-dir", str(out_dir)],
        check=True,
        capture_output=True,
        text=True,
    )
    cert = json.loads((out_dir / "v16_zero_md_role_transfer_readout_certificate.json").read_text())
    kcsa = next(row for row in cert["target_rows"] if row["target_id"] == "KcsA")
    assert kcsa["selected_core_or_clean_abstain"] == "membrane_pore_roles_detected_without_soluble_core_misclassification"
    assert "tetramer_interface_support_later_biological_assembly_required" in kcsa["monitor_only_roles"]
    assert "no_whole_fold_claim" in kcsa["limitations"]
    assert kcsa["forbidden_misclassification_violations"] == []
