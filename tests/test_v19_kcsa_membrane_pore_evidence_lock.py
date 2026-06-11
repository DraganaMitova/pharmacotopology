from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "run_v19_kcsa_membrane_pore_evidence_lock_v0.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("v19_kcsa_lock", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _passing_v19() -> dict:
    return {
        "test_status": "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED",
        "positive_pressure_evidence_found": True,
        "membrane_pore_role_evidence_found": True,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "soluble_core_misclassification_avoided": True,
        "whole_channel_fold_claim_made": False,
        "tetramer_claim_made": False,
        "pore_filter_annotation_present": True,
        "role_buckets_assigned": ["pore_selectivity_filter_core", "transmembrane_helix_scaffold"],
        "available_evidence": ["leakage_guard_against_soluble_core_misread"],
        "missing_evidence": ["external_couplings_if_available"],
        "forbidden_misclassification_violations": [],
        "pore_filter_readout": {
            "sequence_motif_probe": {"selection_threshold_used": False},
            "potassium_ion_probe": {"selection_threshold_used": False},
        },
        "transmembrane_helix_readout": {"selection_threshold_used": False},
        "chain_interface_readout": {"selection_threshold_used": False},
    }


def test_v19_lock_accepts_passed_kcsa_pressure_evidence() -> None:
    mod = _load_module()
    cert = mod.build_lock(_passing_v19())
    assert cert["lock_status"] == "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCKED"
    assert cert["positive_pressure_evidence_found"] is True
    assert cert["membrane_pore_role_evidence_found"] is True
    assert cert["positive_folding_evidence_found"] is False
    assert cert["claim_allowed"] is False
    assert cert["new_md_executed"] is False
    assert cert["membrane_md_executed"] is False
    assert cert["lock_failed_checks"] == []


def test_v19_lock_blocks_if_folding_claim_appears() -> None:
    mod = _load_module()
    bad = _passing_v19()
    bad["positive_folding_evidence_found"] = True
    cert = mod.build_lock(bad)
    assert cert["lock_status"] == "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCK_BLOCKED"
    assert "positive_folding_evidence_forbidden" in cert["lock_failed_checks"]
    assert cert["positive_folding_evidence_found"] is False


def test_v19_lock_writer_outputs_certificate_and_report(tmp_path: Path) -> None:
    mod = _load_module()
    cert = mod.build_lock(_passing_v19())
    out = tmp_path / "out"
    mod.write_outputs(out, cert)
    cpath = out / "v19_kcsa_membrane_pore_evidence_lock_certificate.json"
    rpath = out / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCK_REPORT.md"
    assert cpath.exists()
    assert rpath.exists()
    written = json.loads(cpath.read_text(encoding="utf-8"))
    assert written["lock_status"] == "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCKED"
    assert written["claim_allowed"] is False
