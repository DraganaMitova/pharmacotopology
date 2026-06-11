#!/usr/bin/env python3
from __future__ import annotations

"""Lock the V19 KcsA membrane/pore pressure-evidence milestone.

Reporting/locking only. It consumes the V19 zero-MD KcsA readout and freezes the
interpretation: membrane/pore pressure evidence was found, soluble compact-core
misclassification was avoided, but no whole-channel fold/tetramer/membrane-MD
claim is allowed.
"""

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V19_CERT = RUN_ROOT / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT" / "v19_kcsa_membrane_pore_evidence_readout_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCK"


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def build_lock(v19: dict[str, Any]) -> dict[str, Any]:
    failed: list[str] = []
    if v19.get("test_status") != "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED":
        failed.append("v19_membrane_pore_evidence_readout_passed")
    if v19.get("positive_pressure_evidence_found") is not True:
        failed.append("positive_pressure_evidence_found")
    if v19.get("membrane_pore_role_evidence_found") is not True:
        failed.append("membrane_pore_role_evidence_found")
    if v19.get("positive_folding_evidence_found") is not False:
        failed.append("positive_folding_evidence_forbidden")
    if v19.get("claim_allowed") is not False:
        failed.append("claim_disabled")
    if v19.get("new_md_executed") is not False:
        failed.append("no_new_md")
    if v19.get("membrane_md_executed") is not False:
        failed.append("no_membrane_md")
    if v19.get("fixed_residue_cutoff_used") is not False:
        failed.append("fixed_residue_cutoff_retired")
    if v19.get("native_metrics_used_for_selection") is not False:
        failed.append("native_metric_selection_forbidden")
    if v19.get("soluble_core_misclassification_avoided") is not True:
        failed.append("soluble_core_misclassification_avoided")
    if v19.get("whole_channel_fold_claim_made") is not False:
        failed.append("whole_channel_fold_claim_forbidden")
    if v19.get("tetramer_claim_made") is not False:
        failed.append("tetramer_claim_forbidden")
    if v19.get("forbidden_misclassification_violations") not in ([], None):
        failed.append("no_forbidden_misclassification")

    pore = v19.get("pore_filter_readout") if isinstance(v19.get("pore_filter_readout"), dict) else {}
    motif = pore.get("sequence_motif_probe") if isinstance(pore.get("sequence_motif_probe"), dict) else {}
    ion = pore.get("potassium_ion_probe") if isinstance(pore.get("potassium_ion_probe"), dict) else {}
    helix = v19.get("transmembrane_helix_readout") if isinstance(v19.get("transmembrane_helix_readout"), dict) else {}
    iface = v19.get("chain_interface_readout") if isinstance(v19.get("chain_interface_readout"), dict) else {}
    if v19.get("pore_filter_annotation_present") is not True:
        failed.append("pore_filter_annotation_or_diagnostic_present")
    if motif.get("selection_threshold_used") is not False:
        failed.append("motif_probe_not_fixed_threshold")
    if ion.get("selection_threshold_used") is not False:
        failed.append("ion_probe_not_fixed_threshold")
    if helix.get("selection_threshold_used") is not False:
        failed.append("helix_probe_not_fixed_threshold")
    if iface.get("selection_threshold_used") is not False:
        failed.append("interface_probe_not_fixed_threshold")

    lock_status = (
        "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCKED"
        if not failed
        else "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCK_BLOCKED"
    )
    return {
        "kind": "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCK_v0",
        "run_mode": "lock_only_no_new_simulation_no_threshold_tuning",
        "lock_status": lock_status,
        "source_v19_status": v19.get("test_status"),
        "target_id": "KcsA",
        "role_class": "membrane_pore_oligomer_object",
        "locked_claim": "positive_membrane_pore_pressure_context_evidence_not_folding_claim",
        "positive_pressure_evidence_found": v19.get("positive_pressure_evidence_found") is True,
        "membrane_pore_role_evidence_found": v19.get("membrane_pore_role_evidence_found") is True,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "soluble_core_misclassification_avoided": v19.get("soluble_core_misclassification_avoided") is True,
        "whole_channel_fold_claim_made": False,
        "tetramer_claim_made": False,
        "role_buckets_assigned": v19.get("role_buckets_assigned", []),
        "available_evidence": v19.get("available_evidence", []),
        "missing_evidence": v19.get("missing_evidence", []),
        "pore_filter_readout": pore,
        "transmembrane_helix_readout": helix,
        "chain_interface_readout": iface,
        "forbidden_misclassification_violations": v19.get("forbidden_misclassification_violations", []),
        "lock_failed_checks": failed,
        "locked_interpretation": (
            "KcsA shows zero-MD membrane/pore pressure-context evidence: pore/filter diagnostic, potassium-ion shell, "
            "transmembrane helix scaffold, and chain/interface context are detected while soluble compact-core, "
            "whole-channel fold, tetramer, and membrane-MD claims remain forbidden."
        ),
        "next_step": "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT",
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V19 KcsA Membrane/Pore Evidence Lock",
        "",
        f"Lock status: `{cert.get('lock_status')}`",
        f"Source V19 status: `{cert.get('source_v19_status')}`",
        f"Positive pressure evidence found: `{cert.get('positive_pressure_evidence_found')}`",
        f"Membrane/pore role evidence found: `{cert.get('membrane_pore_role_evidence_found')}`",
        f"Positive folding evidence found: `{cert.get('positive_folding_evidence_found')}`",
        f"Soluble-core misclassification avoided: `{cert.get('soluble_core_misclassification_avoided')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        f"Membrane MD executed: `{cert.get('membrane_md_executed')}`",
        "",
        "## Interpretation",
        cert.get("locked_interpretation", ""),
        "",
        "## Evidence buckets",
        f"`{cert.get('role_buckets_assigned')}`",
        "",
        "## Failed checks",
        f"`{cert.get('lock_failed_checks')}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v19_kcsa_membrane_pore_evidence_lock_certificate.json",
        "report": out_dir / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCK_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v19-cert", type=Path, default=DEFAULT_V19_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    v19 = _read_json(args.v19_cert, "V19 KcsA evidence certificate")
    cert = build_lock(v19)
    write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "lock_status": cert["lock_status"],
        "positive_pressure_evidence_found": cert["positive_pressure_evidence_found"],
        "membrane_pore_role_evidence_found": cert["membrane_pore_role_evidence_found"],
        "positive_folding_evidence_found": cert["positive_folding_evidence_found"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "membrane_md_executed": cert["membrane_md_executed"],
        "lock_failed_checks": cert["lock_failed_checks"],
        "certificate": str(args.out_dir / "v19_kcsa_membrane_pore_evidence_lock_certificate.json"),
        "report": str(args.out_dir / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCK_REPORT.md"),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
