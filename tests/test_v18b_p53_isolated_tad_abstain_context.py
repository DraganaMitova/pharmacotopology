from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v18b_p53_isolated_tad_abstain_context_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v18b_abstain", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _v18_cert() -> dict:
    return {
        "test_status": "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_PASSED_CLAIM_DISABLED",
        "positive_pressure_evidence_found": True,
        "partner_induced_role_evidence_found": True,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
    }


def _lock_cert() -> dict:
    return {
        "lock_status": "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCKED",
        "claim_allowed": False,
        "new_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
    }


def test_v18b_isolated_tad_abstain_preserves_bound_evidence() -> None:
    mod = _load_module()
    cert = mod.build_abstain_context(_v18_cert(), _lock_cert())
    assert cert["test_status"] == "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT_PASSED_CLAIM_DISABLED"
    assert cert["bound_complex_partner_induced_evidence_preserved"] is True
    assert cert["partner_induced_role_evidence_found"] is True
    assert cert["isolated_TAD_autonomous_core_selected"] is False
    assert cert["isolated_TAD_fold_claim_made"] is False
    assert cert["isolated_TAD_clean_abstain_valid"] is True
    assert cert["isolated_TAD_status"] == "clean_abstain_no_isolated_TAD_material_no_autonomous_core_evidence"
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["forbidden_misclassification_violations"] == []
    assert cert["failed_checks"] == []
    assert cert["contrast_status"] == "isolated_TAD_clean_abstain_bound_partner_induced_role_evidence_preserved"


def test_v18b_isolated_tad_abstain_writer(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_abstain_context(_v18_cert(), _lock_cert())
    out = tmp_path / "out"
    mod.write_outputs(out, cert)
    cpath = out / "v18b_p53_isolated_tad_abstain_context_certificate.json"
    rpath = out / "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT_REPORT.md"
    assert cpath.exists()
    assert rpath.exists()
    written = json.loads(cpath.read_text(encoding="utf-8"))
    assert written["test_status"].endswith("PASSED_CLAIM_DISABLED")
    assert written["positive_folding_evidence_found"] is False
