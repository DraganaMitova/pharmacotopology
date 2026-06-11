#!/usr/bin/env python3
from __future__ import annotations

"""V21 pressure-evidence panel summary.

This is a reporting/locking summary only. It consumes the already-produced V18,
V18b, V19, and V20 zero-MD pressure-evidence artifacts and summarizes the panel
as pressure-context evidence, not folding evidence. It must not run MD, must not
upgrade claims, and must keep all pressure classes separated.
"""

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V18_LOCK_CERT = RUN_ROOT / "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCK" / "v18_p53_partner_induced_evidence_lock_certificate.json"
DEFAULT_V18B_CERT = RUN_ROOT / "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT" / "v18b_p53_isolated_tad_abstain_context_certificate.json"
DEFAULT_V19_LOCK_CERT = RUN_ROOT / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCK" / "v19_kcsa_membrane_pore_evidence_lock_certificate.json"
DEFAULT_V20_CERT = RUN_ROOT / "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT" / "v20_xcl1_state_specific_evidence_readout_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY"


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _has_false(row: dict[str, Any], key: str) -> bool:
    return row.get(key) is False


def _has_true(row: dict[str, Any], key: str) -> bool:
    return row.get(key) is True


def _p53_row(v18_lock: dict[str, Any], v18b: dict[str, Any]) -> dict[str, Any]:
    passed = (
        v18_lock.get("lock_status") == "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCKED"
        and v18b.get("test_status") == "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT_PASSED_CLAIM_DISABLED"
        and _has_true(v18_lock, "positive_pressure_evidence_found")
        and _has_true(v18b, "bound_complex_partner_induced_evidence_preserved")
        and _has_false(v18_lock, "positive_folding_evidence_found")
        and _has_false(v18b, "positive_folding_evidence_found")
        and _has_false(v18b, "isolated_TAD_autonomous_core_selected")
        and _has_false(v18b, "isolated_TAD_fold_claim_made")
    )
    missing = sorted(set((v18_lock.get("missing_evidence") or []) + (v18b.get("missing_evidence") or [])))
    return {
        "target_id": "p53_TAD_MDM2",
        "pressure_class": "disorder_partner_induced_binding",
        "role_class": "disorder_partner_induced_object",
        "source_artifacts": ["V18_P53_PARTNER_INDUCED_EVIDENCE_LOCK", "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT"],
        "pressure_evidence_status": "passed" if passed else "blocked_or_incomplete",
        "positive_pressure_evidence_found": bool(passed),
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "evidence_summary": "isolated TAD clean-abstains while MDM2-bound TAD preserves partner-induced interface/helix evidence",
        "selected_core_or_clean_abstain": v18b.get("selected_core_or_clean_abstain"),
        "available_evidence": sorted(set((v18_lock.get("available_evidence") or []) + (v18b.get("available_evidence") or []))),
        "missing_evidence": missing,
        "forbidden_misclassification_violations": sorted(set((v18_lock.get("forbidden_misclassification_violations") or []) + (v18b.get("forbidden_misclassification_violations") or []))),
        "stronger_test_requires": ["isolated_TAD_disorder_prior_or_external_annotation", "external_couplings_if_available"],
        "md_readiness": "no_new_MD_recommended_before_missing_external_context",
    }


def _kcsa_row(v19_lock: dict[str, Any]) -> dict[str, Any]:
    passed = (
        v19_lock.get("lock_status") == "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCKED"
        and _has_true(v19_lock, "positive_pressure_evidence_found")
        and _has_true(v19_lock, "membrane_pore_role_evidence_found")
        and _has_true(v19_lock, "soluble_core_misclassification_avoided")
        and _has_false(v19_lock, "positive_folding_evidence_found")
        and _has_false(v19_lock, "new_md_executed")
        and _has_false(v19_lock, "membrane_md_executed")
        and _has_false(v19_lock, "whole_channel_fold_claim_made")
        and _has_false(v19_lock, "tetramer_claim_made")
    )
    return {
        "target_id": "KcsA",
        "pressure_class": "membrane_pore_oligomer_environment",
        "role_class": "membrane_pore_oligomer_object",
        "source_artifacts": ["V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCK"],
        "pressure_evidence_status": "passed" if passed else "blocked_or_incomplete",
        "positive_pressure_evidence_found": bool(passed),
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "evidence_summary": "pore/filter diagnostic, potassium-ion shell, TM helix scaffold, and chain/interface context detected without soluble-core misclassification",
        "selected_core_or_clean_abstain": "membrane_pore_roles_detected_without_soluble_core_misclassification" if passed else "blocked_or_incomplete",
        "available_evidence": v19_lock.get("available_evidence", []),
        "missing_evidence": v19_lock.get("missing_evidence", []),
        "forbidden_misclassification_violations": v19_lock.get("forbidden_misclassification_violations", []),
        "stronger_test_requires": ["biological_tetramer_assembly_annotation_for_tetramer_claim", "external_couplings_if_available"],
        "md_readiness": "no_membrane_MD_recommended_before_assembly_and_external_coupling_context",
    }


def _xcl1_row(v20: dict[str, Any]) -> dict[str, Any]:
    passed = (
        v20.get("test_status") == "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED"
        and _has_true(v20, "positive_pressure_evidence_found")
        and _has_true(v20, "state_specific_role_evidence_found")
        and _has_true(v20, "state_A_role_evidence_found")
        and _has_true(v20, "state_B_role_evidence_found")
        and _has_false(v20, "positive_folding_evidence_found")
        and _has_false(v20, "mixed_state_pollution")
        and _has_false(v20, "single_fold_claim_made")
        and _has_false(v20, "fold_switch_claim_made")
        and _has_false(v20, "new_md_executed")
    )
    return {
        "target_id": "XCL1_lymphotactin",
        "pressure_class": "metamorphic_fold_switching",
        "role_class": "metamorphic_switch_object",
        "source_artifacts": ["V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT"],
        "pressure_evidence_status": "passed" if passed else "blocked_or_incomplete",
        "positive_pressure_evidence_found": bool(passed),
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "evidence_summary": "state A and state B role evidence kept in separate buckets with no mixed-state fake core or forced single-fold claim",
        "selected_core_or_clean_abstain": v20.get("selected_core_or_clean_abstain"),
        "available_evidence": v20.get("available_evidence", []),
        "missing_evidence": v20.get("missing_evidence", []),
        "forbidden_misclassification_violations": v20.get("forbidden_misclassification_violations", []),
        "stronger_test_requires": ["condition_labels_if_available", "state_specific_external_couplings_or_constraints_if_available"],
        "md_readiness": "no_MD_recommended_before_state_conditions_and_external_constraints",
    }


def build_summary(v18_lock: dict[str, Any], v18b: dict[str, Any], v19_lock: dict[str, Any], v20: dict[str, Any]) -> dict[str, Any]:
    rows = [_p53_row(v18_lock, v18b), _kcsa_row(v19_lock), _xcl1_row(v20)]
    failed_checks: list[str] = []
    for row in rows:
        if row["pressure_evidence_status"] != "passed":
            failed_checks.append(f"{row['target_id']}_pressure_evidence_passed")
        if row.get("positive_folding_evidence_found") is not False:
            failed_checks.append(f"{row['target_id']}_positive_folding_evidence_forbidden")
        if row.get("claim_allowed") is not False:
            failed_checks.append(f"{row['target_id']}_claim_disabled")
        if row.get("new_md_executed") is not False:
            failed_checks.append(f"{row['target_id']}_no_new_md")
        if row.get("forbidden_misclassification_violations") not in ([], None):
            failed_checks.append(f"{row['target_id']}_no_forbidden_misclassification")

    pressure_targets = [row["target_id"] for row in rows if row["positive_pressure_evidence_found"]]
    status = "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_LOCKED" if not failed_checks else "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_BLOCKED"
    all_missing = sorted({item for row in rows for item in row.get("missing_evidence", [])})
    return {
        "kind": "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_v0",
        "run_mode": "summary_lock_only_no_new_simulation_no_threshold_tuning",
        "summary_status": status,
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "positive_pressure_evidence_targets": pressure_targets,
        "positive_folding_evidence_targets": [],
        "pressure_classes_with_evidence": [row["pressure_class"] for row in rows if row["positive_pressure_evidence_found"]],
        "pressure_evidence_count": len(pressure_targets),
        "target_rows": rows,
        "global_missing_evidence": all_missing,
        "md_readiness_decision": "no_new_MD_recommended_yet_external_evidence_and_annotations_first",
        "next_recommended_panel": "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL",
        "summary_failed_checks": failed_checks,
        "locked_interpretation": (
            "V21 summarizes three zero-MD pressure-context evidence wins: p53/MDM2 disorder-partner pressure, "
            "KcsA membrane/pore pressure, and XCL1 metamorphic state-specific pressure. This is not a universal folding claim; "
            "positive_folding_evidence_targets remains empty and claim_allowed remains false."
        ),
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V21 Pressure Evidence Panel Summary",
        "",
        f"Summary status: `{cert.get('summary_status')}`",
        f"Positive pressure evidence targets: `{cert.get('positive_pressure_evidence_targets')}`",
        f"Positive folding evidence targets: `{cert.get('positive_folding_evidence_targets')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        f"Fixed residue cutoff used: `{cert.get('fixed_residue_cutoff_used')}`",
        "",
        "## Locked interpretation",
        cert.get("locked_interpretation", ""),
        "",
        "## Target rows",
    ]
    for row in cert.get("target_rows", []):
        lines += [
            f"### {row.get('target_id')}",
            f"- Pressure class: `{row.get('pressure_class')}`",
            f"- Evidence status: `{row.get('pressure_evidence_status')}`",
            f"- Summary: {row.get('evidence_summary')}",
            f"- Missing evidence: `{row.get('missing_evidence')}`",
            f"- MD readiness: `{row.get('md_readiness')}`",
            "",
        ]
    lines += [
        "## Next panel",
        f"`{cert.get('next_recommended_panel')}`",
        "",
        "## Failed checks",
        f"`{cert.get('summary_failed_checks')}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v21_pressure_evidence_panel_summary_certificate.json",
        "report": out_dir / "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_REPORT.md",
        "table": out_dir / "v21_pressure_evidence_panel_table.csv",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)
    lines = ["target_id,pressure_class,pressure_evidence_status,positive_pressure_evidence_found,positive_folding_evidence_found,claim_allowed,new_md_executed,md_readiness"]
    for row in cert.get("target_rows", []):
        def esc(v: Any) -> str:
            return str(v).replace(',', ';')
        lines.append(
            ",".join([
                esc(row.get("target_id")), esc(row.get("pressure_class")), esc(row.get("pressure_evidence_status")),
                esc(row.get("positive_pressure_evidence_found")), esc(row.get("positive_folding_evidence_found")),
                esc(row.get("claim_allowed")), esc(row.get("new_md_executed")), esc(row.get("md_readiness")),
            ])
        )
    paths["table"].write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v18-lock-cert", type=Path, default=DEFAULT_V18_LOCK_CERT)
    parser.add_argument("--v18b-cert", type=Path, default=DEFAULT_V18B_CERT)
    parser.add_argument("--v19-lock-cert", type=Path, default=DEFAULT_V19_LOCK_CERT)
    parser.add_argument("--v20-cert", type=Path, default=DEFAULT_V20_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    cert = build_summary(
        _read_json(args.v18_lock_cert, "V18 p53 lock certificate"),
        _read_json(args.v18b_cert, "V18b p53 contrast certificate"),
        _read_json(args.v19_lock_cert, "V19 KcsA lock certificate"),
        _read_json(args.v20_cert, "V20 XCL1 certificate"),
    )
    write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "summary_status": cert["summary_status"],
        "positive_pressure_evidence_targets": cert["positive_pressure_evidence_targets"],
        "positive_folding_evidence_targets": cert["positive_folding_evidence_targets"],
        "pressure_evidence_count": cert["pressure_evidence_count"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "md_readiness_decision": cert["md_readiness_decision"],
        "next_recommended_panel": cert["next_recommended_panel"],
        "summary_failed_checks": cert["summary_failed_checks"],
        "certificate": str(args.out_dir / "v21_pressure_evidence_panel_summary_certificate.json"),
        "report": str(args.out_dir / "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_REPORT.md"),
        "table": str(args.out_dir / "v21_pressure_evidence_panel_table.csv"),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
