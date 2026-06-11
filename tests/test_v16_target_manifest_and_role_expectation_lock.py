from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def _base_v15_lock() -> dict:
    return {
        "kind": "V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCKED_v0",
        "lock_status": "V15_THREE_PROTEIN_DYNAMIC_GRAMMAR_PANEL_LOCKED",
        "claim_allowed": False,
        "locked_claim": "unified_role_aware_evidence_grammar_three_protein_panel_locked_not_folding_solved",
    }


def test_v16_manifest_file_is_locked_transfer_manifest() -> None:
    manifest = json.loads((REPO_ROOT / "data" / "v16_locked_grammar_transfer_target_manifest.json").read_text())
    assert manifest["kind"] == "V16_LOCKED_GRAMMAR_TRANSFER_TARGET_MANIFEST_v0"
    assert manifest["claim_allowed"] is False
    policy = manifest["panel_policy"]
    assert policy["not_new_tuning_panel"] is True
    assert policy["not_proof_of_solved_folding"] is True
    assert policy["new_md_allowed_in_manifest_lock"] is False
    assert policy["data_download_allowed_in_manifest_lock"] is False
    assert policy["grammar_changes_allowed"] is False
    assert policy["target_specific_threshold_tuning_allowed"] is False
    assert policy["fixed_threshold_policy"] == "forbidden"
    assert policy["native_metrics_not_used_for_selection"] is True
    assert policy["fixed_residue_cutoff_used"] is False
    ids = {target["target_id"] for target in manifest["targets"]}
    assert ids == {"p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}


def test_v16_lock_script_locks_manifest_without_md_or_preflight(tmp_path: Path) -> None:
    manifest_path = REPO_ROOT / "data" / "v16_locked_grammar_transfer_target_manifest.json"
    v15_path = tmp_path / "v15_lock.json"
    out_dir = tmp_path / "out"
    v15_path.write_text(json.dumps(_base_v15_lock()))

    proc = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "run_v16_target_manifest_and_role_expectation_lock_v0.py"),
            "--manifest",
            str(manifest_path),
            "--v15-lock-cert",
            str(v15_path),
            "--out-dir",
            str(out_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCKED" in proc.stdout
    cert = json.loads((out_dir / "v16_target_manifest_and_role_expectation_lock_certificate.json").read_text())
    assert cert["lock_status"] == "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCKED"
    assert cert["claim_allowed"] is False
    assert cert["data_preflight_executed"] is False
    assert cert["new_md_executed"] is False
    assert cert["lock_failed_checks"] == []
    assert cert["fixed_threshold_policy"] == "forbidden"
    assert cert["native_metrics_not_used_for_selection"] is True


def test_v16_lock_blocks_if_v15_three_protein_lock_missing(tmp_path: Path) -> None:
    manifest_path = REPO_ROOT / "data" / "v16_locked_grammar_transfer_target_manifest.json"
    v15_path = tmp_path / "v15_partial.json"
    out_dir = tmp_path / "out_blocked"
    v15 = _base_v15_lock()
    v15["lock_status"] = "V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCKED"
    v15_path.write_text(json.dumps(v15))

    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "run_v16_target_manifest_and_role_expectation_lock_v0.py"),
            "--manifest",
            str(manifest_path),
            "--v15-lock-cert",
            str(v15_path),
            "--out-dir",
            str(out_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    cert = json.loads((out_dir / "v16_target_manifest_and_role_expectation_lock_certificate.json").read_text())
    assert cert["lock_status"] == "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCK_BLOCKED"
    assert "v15_three_protein_panel_locked" in cert["lock_failed_checks"]
    assert cert["claim_allowed"] is False


def test_v16_targets_have_abstain_and_forbidden_misclassification() -> None:
    manifest = json.loads((REPO_ROOT / "data" / "v16_locked_grammar_transfer_target_manifest.json").read_text())
    for target in manifest["targets"]:
        assert target["claim_allowed"] is False
        assert target.get("clean_abstain_allowed") is True or any(
            state.get("clean_abstain_allowed") is True for state in target.get("states", [])
        )
        assert target["allowed_evidence_roles"]
        assert target["forbidden_misclassification"]
        assert target["required_inputs"]
