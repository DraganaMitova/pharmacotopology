from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PREFLIGHT_SCRIPT = REPO_ROOT / "scripts" / "run_v16_transfer_data_preflight_v0.py"
ZERO_MD_SCRIPT = REPO_ROOT / "scripts" / "run_v16_zero_md_role_transfer_readout_v0.py"
GAP_LOCK_SCRIPT = REPO_ROOT / "scripts" / "run_v16_pressure_evidence_gap_lock_v0.py"
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


def _zero_md_cert(tmp_path: Path) -> Path:
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
                _write_mock_pdb(mat_path, {"A": 219, "B": 212, "C": 103})
            else:
                _write_mock_pdb(mat_path, {"A": 120, "B": 120})
            _write_prov(prov_path, target["target_id"], mat["pdb_id"])
    material_path = tmp_path / "material.json"
    material_path.write_text(json.dumps(material), encoding="utf-8")
    lock_path = tmp_path / "v16_lock.json"
    lock_path.write_text(json.dumps(_v16_lock()), encoding="utf-8")
    preflight_dir = tmp_path / "preflight"
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
            str(preflight_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    zero_dir = tmp_path / "zero_md"
    subprocess.run(
        [
            sys.executable,
            str(ZERO_MD_SCRIPT),
            "--role-manifest",
            str(ROLE_MANIFEST),
            "--preflight-cert",
            str(preflight_dir / "v16_transfer_data_preflight_certificate.json"),
            "--out-dir",
            str(zero_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return zero_dir / "v16_zero_md_role_transfer_readout_certificate.json"


def _gap_lock(tmp_path: Path) -> dict:
    zero_cert = _zero_md_cert(tmp_path)
    out_dir = tmp_path / "gap"
    subprocess.run(
        [sys.executable, str(GAP_LOCK_SCRIPT), "--zero-md-cert", str(zero_cert), "--out-dir", str(out_dir)],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads((out_dir / "v16_pressure_evidence_gap_lock_certificate.json").read_text())


def test_v16_pressure_evidence_gap_lock_is_not_folding_evidence(tmp_path: Path) -> None:
    cert = _gap_lock(tmp_path)
    assert cert["gap_lock_status"] == "V16_PRESSURE_EVIDENCE_GAP_LOCKED"
    assert set(cert["role_classification_passed_targets"]) == {"p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}
    assert cert["positive_folding_evidence_targets"] == []
    assert cert["positive_contact_evidence_targets"] == []
    assert cert["evidence_claim_allowed_targets"] == []
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["fixed_residue_cutoff_used"] is False
    assert cert["native_metrics_used_for_selection"] is False
    assert cert["global_missing_layer"] == "evidence_layer_not_grammar_or_role_classification_layer"


def test_v16_pressure_evidence_gap_lock_target_specific_missing_layers(tmp_path: Path) -> None:
    cert = _gap_lock(tmp_path)
    rows = {row["target_id"]: row for row in cert["target_gap_rows"]}
    assert "isolated_TAD_disorder_or_clean_abstain_context" in rows["p53_TAD_MDM2"]["missing_for_evidence_test"]
    assert rows["p53_TAD_MDM2"]["next_evidence_test"] == "partner_induced_interface_or_helix_readout"
    assert "membrane_topology_annotation" in rows["KcsA"]["missing_for_evidence_test"]
    assert rows["KcsA"]["next_evidence_test"] == "membrane_pore_role_evidence_readout"
    assert "leakage_guard_preventing_mixed_state_fake_core" in rows["XCL1_lymphotactin"]["missing_for_evidence_test"]
    assert rows["XCL1_lymphotactin"]["next_evidence_test"] == "state_specific_role_separation_readout"


def test_v16_pressure_evidence_gap_lock_blocks_if_zero_md_claims_folding(tmp_path: Path) -> None:
    zero_cert = json.loads(_zero_md_cert(tmp_path).read_text())
    zero_cert["positive_folding_evidence_targets"] = ["KcsA"]
    bad_path = tmp_path / "bad_zero_md.json"
    bad_path.write_text(json.dumps(zero_cert), encoding="utf-8")
    out_dir = tmp_path / "gap_bad"
    subprocess.run(
        [sys.executable, str(GAP_LOCK_SCRIPT), "--zero-md-cert", str(bad_path), "--out-dir", str(out_dir)],
        check=True,
        capture_output=True,
        text=True,
    )
    cert = json.loads((out_dir / "v16_pressure_evidence_gap_lock_certificate.json").read_text())
    assert cert["gap_lock_status"] == "V16_PRESSURE_EVIDENCE_GAP_LOCK_BLOCKED"
    assert "positive_folding_evidence_targets_empty" in cert["gap_lock_failed_checks"]
    assert cert["claim_allowed"] is False


def test_v16_zero_md_certificate_uses_classification_naming(tmp_path: Path) -> None:
    zero = json.loads(_zero_md_cert(tmp_path).read_text())
    assert set(zero["role_classification_passed_targets"]) == {"p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}
    assert set(zero["pressure_role_transfer_passed_targets"]) == {"p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}
    assert zero["positive_folding_evidence_targets"] == []
    assert "positive_role_targets" not in zero
