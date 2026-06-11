from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "run_v16_transfer_data_preflight_v0.py"
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
            # Fixed-width-enough PDB ATOM line with CA atom and chain at column 22.
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


def _run(tmp_path: Path, material_manifest: Path) -> dict:
    v16_lock = tmp_path / "v16_lock.json"
    out_dir = tmp_path / "out"
    v16_lock.write_text(json.dumps(_v16_lock()), encoding="utf-8")
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--role-manifest",
            str(ROLE_MANIFEST),
            "--material-manifest",
            str(material_manifest),
            "--v16-lock-cert",
            str(v16_lock),
            "--out-dir",
            str(out_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads((out_dir / "v16_transfer_data_preflight_certificate.json").read_text())


def test_v16_material_manifest_is_preflight_only() -> None:
    manifest = json.loads(MATERIAL_MANIFEST.read_text())
    assert manifest["kind"] == "V16_TRANSFER_DATA_MATERIAL_MANIFEST_v0"
    assert manifest["claim_allowed"] is False
    assert manifest["data_preflight_only"] is True
    assert manifest["new_md_allowed"] is False
    assert manifest["grammar_changes_allowed"] is False
    assert manifest["target_specific_threshold_tuning_allowed"] is False
    assert manifest["native_metrics_not_used_for_selection"] is True
    assert manifest["fixed_residue_cutoff_used"] is False
    ids = {target["target_id"] for target in manifest["targets"]}
    assert ids == {"p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}


def test_v16_preflight_reports_missing_material_without_downloading(tmp_path: Path) -> None:
    # Use the repo manifest as-is; target files are intentionally not packaged.
    cert = _run(tmp_path, MATERIAL_MANIFEST)
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["download_attempted_by_this_script"] is False
    assert cert["data_preflight_status"] in {
        "V16_TRANSFER_DATA_PREFLIGHT_BLOCKED_MISSING_REQUIRED_MATERIAL",
        "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT",
    }
    assert cert["fixed_residue_cutoff_used"] is False


def test_v16_preflight_ready_with_mock_pdb_and_provenance(tmp_path: Path) -> None:
    material = json.loads(MATERIAL_MANIFEST.read_text())
    for target in material["targets"]:
        for mat in target["required_material"]:
            mat_path = tmp_path / mat["path"]
            prov_path = tmp_path / mat["provenance_path"]
            mat["path"] = str(mat_path)
            mat["provenance_path"] = str(prov_path)
            if target["target_id"] == "p53_TAD_MDM2":
                _write_mock_pdb(mat_path, {"A": 30, "B": 10})
            elif target["target_id"] == "KcsA":
                _write_mock_pdb(mat_path, {"A": 30, "B": 30, "C": 30, "D": 30})
            else:
                _write_mock_pdb(mat_path, {"A": 50})
            _write_prov(prov_path, target["target_id"], mat["pdb_id"])
    material_path = tmp_path / "material.json"
    material_path.write_text(json.dumps(material), encoding="utf-8")

    cert = _run(tmp_path, material_path)
    assert cert["data_preflight_status"] == "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT"
    assert set(cert["ready_targets"]) == {"p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}
    assert cert["blocked_targets"] == []
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False


def test_v16_preflight_blocks_if_v16_lock_missing(tmp_path: Path) -> None:
    v16_lock = tmp_path / "v16_lock_bad.json"
    out_dir = tmp_path / "out"
    v16_lock.write_text(json.dumps({"lock_status": "BROKEN", "claim_allowed": False}), encoding="utf-8")
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--role-manifest",
            str(ROLE_MANIFEST),
            "--material-manifest",
            str(MATERIAL_MANIFEST),
            "--v16-lock-cert",
            str(v16_lock),
            "--out-dir",
            str(out_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    cert = json.loads((out_dir / "v16_transfer_data_preflight_certificate.json").read_text())
    assert cert["data_preflight_status"] == "V16_TRANSFER_DATA_PREFLIGHT_BLOCKED_MANIFEST_OR_LOCK_FAILURE"
    assert "v16_manifest_lock_present" in cert["preflight_failed_checks"]
    assert cert["claim_allowed"] is False
