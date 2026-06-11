from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_v15_lock_script_locks_claim_disabled_partial_panel(tmp_path: Path) -> None:
    cert = tmp_path / "v15.json"
    out = tmp_path / "out"
    cert.write_text(json.dumps({
        "global_status": "dynamic_separation_grammar_positive_on_1UBQ_1CLL_4AKE_missing_claim_disabled",
        "fixed_residue_cutoff_used": False,
        "claim_allowed": False,
        "positive_evidence_proteins": ["1UBQ", "1CLL"],
        "missing_artifacts": ["4AKE"],
        "coherence_checks": {
            "no_fixed_residue_cutoff_used_anywhere": True,
            "no_claim_allowed_anywhere": True,
            "positive_evidence_proteins": ["1UBQ", "1CLL"],
            "missing_artifacts": ["4AKE"],
        },
        "protein_rows": [],
    }))
    proc = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "run_v15_dynamic_role_grammar_panel_lock_v0.py"), "--v15-cert", str(cert), "--out-dir", str(out)],
        check=True,
        text=True,
        capture_output=True,
    )
    assert "V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCKED" in proc.stdout
    locked = json.loads((out / "v15_dynamic_role_grammar_panel_locked_certificate.json").read_text())
    assert locked["lock_status"] == "V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCKED"
    assert locked["claim_allowed"] is False
    assert locked["positive_evidence_proteins"] == ["1UBQ", "1CLL"]


def test_v15_lock_script_locks_three_protein_panel(tmp_path: Path) -> None:
    cert = tmp_path / "v15_three.json"
    out = tmp_path / "out_three"
    cert.write_text(json.dumps({
        "global_status": "dynamic_separation_grammar_coherent_across_three_object_types_claim_disabled",
        "fixed_residue_cutoff_used": False,
        "claim_allowed": False,
        "positive_evidence_proteins": ["4AKE", "1UBQ", "1CLL"],
        "missing_artifacts": [],
        "coherence_checks": {
            "no_fixed_residue_cutoff_used_anywhere": True,
            "no_claim_allowed_anywhere": True,
            "positive_evidence_proteins": ["4AKE", "1UBQ", "1CLL"],
            "missing_artifacts": [],
            "bridge_pending_artifacts": [],
        },
        "protein_rows": [],
    }))
    proc = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts" / "run_v15_dynamic_role_grammar_panel_lock_v0.py"), "--v15-cert", str(cert), "--out-dir", str(out)],
        check=True,
        text=True,
        capture_output=True,
    )
    assert "V15_THREE_PROTEIN_DYNAMIC_GRAMMAR_PANEL_LOCKED" in proc.stdout
    locked = json.loads((out / "v15_dynamic_role_grammar_panel_locked_certificate.json").read_text())
    assert locked["lock_status"] == "V15_THREE_PROTEIN_DYNAMIC_GRAMMAR_PANEL_LOCKED"
    assert locked["locked_claim"] == "unified_role_aware_evidence_grammar_three_protein_panel_locked_not_folding_solved"
    assert locked["claim_allowed"] is False
