#!/usr/bin/env python3
from __future__ import annotations

"""V16 pressure-evidence gap lock.

This stage freezes what is still missing before the V16 pressure-role targets can
become real folding/contact evidence tests. It does not run MD, does not tune
thresholds, and does not upgrade zero-MD role classification into folding
evidence.
"""

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_ZERO_MD_CERT = RUN_ROOT / "V16_ZERO_MD_ROLE_TRANSFER_READOUT" / "v16_zero_md_role_transfer_readout_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V16_PRESSURE_EVIDENCE_GAP_LOCK"
REQUIRED_TARGET_IDS = {"p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}

TARGET_GAPS: dict[str, dict[str, Any]] = {
    "p53_TAD_MDM2": {
        "role_class": "disorder_partner_induced_object",
        "forbidden_misclassification_avoided": "autonomous_TAD_compact_fold",
        "missing_for_evidence_test": [
            "isolated_TAD_disorder_or_clean_abstain_context",
            "partner_bound_complex_context",
            "interface_or_contact_evidence",
            "state_labels_isolated_vs_MDM2_bound",
            "external_couplings_if_available",
            "leakage_guard_autonomous_fold_vs_partner_induced_fold",
        ],
        "next_evidence_test": "partner_induced_interface_or_helix_readout",
        "v17_readiness": "ready_for_partner_complex_context_gap_known_not_evidence_claim",
        "evidence_claim_forbidden_until": [
            "isolated_TAD_abstain_context_present",
            "interface_contact_evidence_present",
            "state_label_guard_present",
        ],
    },
    "KcsA": {
        "role_class": "membrane_pore_oligomer_object",
        "forbidden_misclassification_avoided": "soluble_single_domain_compact_core",
        "missing_for_evidence_test": [
            "membrane_topology_annotation",
            "pore_selectivity_filter_annotation",
            "oligomer_or_tetramer_assembly_context",
            "chain_and_interface_identity",
            "transmembrane_helix_roles",
            "external_couplings_if_available",
            "leakage_guard_against_soluble_core_misread",
        ],
        "next_evidence_test": "membrane_pore_role_evidence_readout",
        "v17_readiness": "ready_for_membrane_pore_context_gap_known_not_membrane_MD_claim",
        "evidence_claim_forbidden_until": [
            "membrane_topology_annotation_present",
            "assembly_or_interface_context_present",
            "pore_filter_residue_context_present",
        ],
    },
    "XCL1_lymphotactin": {
        "role_class": "metamorphic_switch_object",
        "forbidden_misclassification_avoided": "single_fold_forcing",
        "missing_for_evidence_test": [
            "state_A_and_state_B_labels",
            "monomer_dimer_context",
            "condition_labels_if_available",
            "state_specific_structural_or_contact_evidence",
            "state_specific_external_couplings_or_constraints_if_available",
            "leakage_guard_preventing_mixed_state_fake_core",
        ],
        "next_evidence_test": "state_specific_role_separation_readout",
        "v17_readiness": "ready_for_two_state_context_gap_known_not_fold_switch_claim",
        "evidence_claim_forbidden_until": [
            "state_specific_labels_present",
            "state_specific_evidence_present",
            "mixed_state_leakage_guard_present",
        ],
    },
}


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # pragma: no cover
        raise SystemExit(f"invalid JSON in {label}: {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _target_rows_by_id(zero_md: dict[str, Any]) -> dict[str, dict[str, Any]]:
    rows = zero_md.get("target_rows") if isinstance(zero_md.get("target_rows"), list) else []
    return {str(row.get("target_id")): row for row in rows if isinstance(row, dict)}


def _role_passed_targets(zero_md: dict[str, Any]) -> list[str]:
    # v22 conservative name first; fallback accepts v21 certificates if users run gap lock before rerunning zero-MD.
    for key in ("role_classification_passed_targets", "pressure_role_transfer_passed_targets", "positive_role_targets"):
        value = zero_md.get(key)
        if isinstance(value, list):
            return [str(x) for x in value]
    return []


def _gap_row(target_id: str, zero_row: dict[str, Any] | None, role_passed: set[str]) -> dict[str, Any]:
    base = TARGET_GAPS[target_id]
    zero_row = zero_row or {}
    forbidden_violations = zero_row.get("forbidden_misclassification_violations")
    if not isinstance(forbidden_violations, list):
        forbidden_violations = []
    passed = target_id in role_passed and not forbidden_violations
    return {
        "target_id": target_id,
        "role_class": base["role_class"],
        "pressure_class": zero_row.get("pressure_class"),
        "role_transfer_status": zero_row.get("target_transfer_status", "not_evaluated"),
        "role_classification_passed": passed,
        "zero_md_classification_only": True,
        "positive_folding_evidence_found": False,
        "forbidden_misclassification_avoided": not forbidden_violations,
        "forbidden_misclassification_guard": base["forbidden_misclassification_avoided"],
        "forbidden_misclassification_violations": forbidden_violations,
        "missing_for_evidence_test": base["missing_for_evidence_test"],
        "evidence_claim_forbidden_until": base["evidence_claim_forbidden_until"],
        "next_evidence_test": base["next_evidence_test"],
        "v17_readiness": base["v17_readiness"],
        "claim_allowed": False,
        "new_md_executed": False,
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "notes": [
            "zero_MD_pressure_role_transfer_is_classification_not_folding_evidence",
            "do_not_upgrade_role_context_to_contact_or_fold_claim",
        ],
    }


def build_gap_lock(zero_md: dict[str, Any]) -> dict[str, Any]:
    rows_by_id = _target_rows_by_id(zero_md)
    target_ids = set(rows_by_id)
    role_passed = set(_role_passed_targets(zero_md))
    gap_rows = [_gap_row(tid, rows_by_id.get(tid), role_passed) for tid in sorted(REQUIRED_TARGET_IDS)]

    violations = zero_md.get("forbidden_misclassification_violations")
    if not isinstance(violations, dict):
        violations = {}

    checks = {
        "zero_md_readout_completed": zero_md.get("role_transfer_status") == "V16_ZERO_MD_ROLE_TRANSFER_READOUT_COMPLETED_CLAIM_DISABLED",
        "claim_disabled": zero_md.get("claim_allowed") is False,
        "no_new_md": zero_md.get("new_md_executed") is False,
        "fixed_residue_cutoff_retired": zero_md.get("fixed_residue_cutoff_used") is False,
        "native_metrics_not_used_for_selection": zero_md.get("native_metrics_used_for_selection") is False,
        "target_specific_threshold_tuning_forbidden": zero_md.get("target_specific_threshold_tuning_allowed") is False,
        "required_targets_present": REQUIRED_TARGET_IDS.issubset(target_ids),
        "role_classification_passed_all_targets": REQUIRED_TARGET_IDS.issubset(role_passed),
        "forbidden_misclassification_avoided": not violations and all(r["forbidden_misclassification_avoided"] for r in gap_rows),
        "positive_folding_evidence_targets_empty": zero_md.get("positive_folding_evidence_targets", []) == [],
    }
    failed = [key for key, ok in checks.items() if not ok]
    status = "V16_PRESSURE_EVIDENCE_GAP_LOCKED" if not failed else "V16_PRESSURE_EVIDENCE_GAP_LOCK_BLOCKED"

    return {
        "kind": "V16_PRESSURE_EVIDENCE_GAP_LOCK_v0",
        "run_mode": "gap_lock_only_no_new_simulation_no_threshold_tuning",
        "gap_lock_status": status,
        "gap_lock_failed_checks": failed,
        "gap_lock_checks": checks,
        "source_zero_md_status": zero_md.get("role_transfer_status"),
        "source_v16_lock_status": zero_md.get("source_v16_lock_status"),
        "role_classification_passed_targets": sorted(role_passed),
        "pressure_role_transfer_passed_targets": sorted(role_passed),
        "positive_folding_evidence_targets": [],
        "positive_contact_evidence_targets": [],
        "evidence_claim_allowed_targets": [],
        "claim_allowed": False,
        "new_md_executed": False,
        "data_preflight_executed": True,
        "zero_md_role_transfer_executed": True,
        "fixed_threshold_policy": "forbidden",
        "threshold_policy": "purpose_aware_auto_selection_only_after_evidence_manifest_lock",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "native_metrics_not_used_for_selection": True,
        "global_missing_layer": "evidence_layer_not_grammar_or_role_classification_layer",
        "next_panel": "V17_PRESSURE_EVIDENCE_PANEL",
        "locked_interpretation": (
            "V16 passed zero-MD pressure-role transfer as role classification, not folding evidence. "
            "The gap lock freezes exactly what is missing before each pressure target can become a real folding/contact evidence test. "
            "No MD, threshold tuning, native-metric selection, fixed residue cutoff, or claim upgrade is allowed here."
        ),
        "target_gap_rows": gap_rows,
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V16 Pressure Evidence Gap Lock",
        "",
        "This is a gap lock, not a new evidence claim.",
        "",
        f"Status: `{cert.get('gap_lock_status')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        f"Role-classification passed targets: `{cert.get('role_classification_passed_targets')}`",
        f"Positive folding-evidence targets: `{cert.get('positive_folding_evidence_targets')}`",
        "",
        "## Locked interpretation",
        "",
        str(cert.get("locked_interpretation")),
        "",
        "## Target gaps",
    ]
    for row in cert.get("target_gap_rows", []):
        lines.extend([
            f"### {row.get('target_id')}",
            f"- Role class: `{row.get('role_class')}`",
            f"- Role transfer status: `{row.get('role_transfer_status')}`",
            f"- Role classification passed: `{row.get('role_classification_passed')}`",
            f"- Positive folding evidence found: `{row.get('positive_folding_evidence_found')}`",
            f"- Forbidden misclassification avoided: `{row.get('forbidden_misclassification_avoided')}`",
            f"- Next evidence test: `{row.get('next_evidence_test')}`",
            "- Missing for evidence test:",
        ])
        for item in row.get("missing_for_evidence_test", []):
            lines.append(f"  - `{item}`")
        lines.extend(["- Evidence claim forbidden until:"])
        for item in row.get("evidence_claim_forbidden_until", []):
            lines.append(f"  - `{item}`")
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Lock V16 pressure evidence gaps after zero-MD role transfer.")
    parser.add_argument("--zero-md-cert", default=str(DEFAULT_ZERO_MD_CERT))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    args = parser.parse_args()

    zero_md = _read_json(Path(args.zero_md_cert), "V16 zero-MD role transfer certificate")
    cert = build_gap_lock(zero_md)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v16_pressure_evidence_gap_lock_certificate.json"
    report_path = out_dir / "V16_PRESSURE_EVIDENCE_GAP_LOCK_REPORT.md"
    cert_path.write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(report_path, cert)

    print(json.dumps({
        "kind": cert.get("kind"),
        "certificate": str(cert_path),
        "report": str(report_path),
        "gap_lock_status": cert.get("gap_lock_status"),
        "gap_lock_failed_checks": cert.get("gap_lock_failed_checks"),
        "role_classification_passed_targets": cert.get("role_classification_passed_targets"),
        "positive_folding_evidence_targets": cert.get("positive_folding_evidence_targets"),
        "next_panel": cert.get("next_panel"),
        "claim_allowed": False,
        "new_md_executed": False,
    }, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
