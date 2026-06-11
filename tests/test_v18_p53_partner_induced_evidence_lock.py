from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v18_p53_partner_induced_evidence_lock_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v18_lock", SCRIPT)
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
        "isolated_TAD_autonomous_fold_status": "clean_abstain_no_isolated_TAD_material_no_autonomous_fold_claim",
        "selected_core_or_clean_abstain": "partner_induced_interface_or_helix_signal_found",
        "available_evidence": ["partner_bound_complex_context", "interface_or_contact_evidence"],
        "missing_evidence": ["isolated_TAD_disorder_or_clean_abstain_context"],
        "forbidden_misclassification_violations": [],
        "interface_contact_readout": {
            "interface_contact_evidence_present": True,
            "selection_threshold_used": False,
            "chain_pair": "A-B",
            "min_ca_distance": 4.884,
            "multi_radius_contact_counts": {"6.0": 7, "8.0": 28, "10.0": 96},
            "contact_probe_policy": "multi_radius_report_only_no_single_fixed_selection_threshold",
        },
        "partner_bound_helix_proxy_readout": {"selection_threshold_used": False, "small_chain": "B"},
    }


def test_v18_partner_induced_evidence_lock_passes_without_claim_upgrade() -> None:
    mod = _load_module()
    cert = mod.build_lock(_v18_cert())
    assert cert["lock_status"] == "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCKED"
    assert cert["positive_pressure_evidence_found"] is True
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["fixed_residue_cutoff_used"] is False
    assert cert["native_metrics_used_for_selection"] is False
    assert cert["autonomous_TAD_fold_claim"] == "forbidden_not_made"
    assert cert["lock_failed_checks"] == []
    assert cert["next_step"] == "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT"


def test_v18_partner_induced_evidence_lock_writer(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_lock(_v18_cert())
    out = tmp_path / "out"
    mod.write_outputs(out, cert)
    cpath = out / "v18_p53_partner_induced_evidence_lock_certificate.json"
    rpath = out / "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCK_REPORT.md"
    assert cpath.exists()
    assert rpath.exists()
    written = json.loads(cpath.read_text(encoding="utf-8"))
    assert written["lock_status"] == "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCKED"
    assert written["positive_folding_evidence_found"] is False
