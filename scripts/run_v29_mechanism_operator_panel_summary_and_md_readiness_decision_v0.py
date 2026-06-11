#!/usr/bin/env python3
from __future__ import annotations

"""V29 mechanism operator panel summary + MD readiness decision.

This is a synthesis layer after V25/V28. It summarizes the mechanism-operator
alphabet across the six active targets and makes an explicit no/yes MD-readiness
decision. It does not run MD, does not tune thresholds, and does not upgrade any
pressure-context evidence into a folding claim.
"""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V25_CERT = RUN_ROOT / "V25_FAST_MECHANISM_EVIDENCE_SPRINT" / "v25_fast_mechanism_evidence_sprint_certificate.json"
DEFAULT_V28_CERT = RUN_ROOT / "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST" / "v28_xcl1_state_condition_evidence_contrast_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_DECISION"

CORE_OPERATORS = ["role_detect", "evidence_assign", "pollution_guard"]
PRESSURE_OPERATORS = ["context_detect", "interface_context", "state_separation", "clean_abstain", "topology_context"]
FORBIDDEN_FLAGS = [
    "claim_allowed",
    "new_md_executed",
    "membrane_md_executed",
    "fixed_residue_cutoff_used",
    "native_metrics_used_for_selection",
]


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _operator_targets(v25: dict[str, Any]) -> list[dict[str, Any]]:
    table = v25.get("unified_mechanism_operator_table") if isinstance(v25.get("unified_mechanism_operator_table"), dict) else {}
    rows = table.get("targets") if isinstance(table.get("targets"), list) else []
    return [dict(row) for row in rows if isinstance(row, dict)]


def _target_row(rows: list[dict[str, Any]], target: str) -> dict[str, Any]:
    for row in rows:
        if row.get("target") == target:
            return row
    return {}


def _row_md_readiness(row: dict[str, Any]) -> dict[str, Any]:
    target = str(row.get("target"))
    missing = [str(x) for x in _as_list(row.get("missing_next_evidence"))]
    operators = [str(x) for x in _as_list(row.get("mechanism_operator_used"))]
    claim_allowed = row.get("claim_allowed") is True

    # MD is not ready unless external constraints/couplings or the target-specific
    # missing mechanism evidence have been locked. Here we deliberately prefer a
    # blocked/readiness explanation over a premature simulation recommendation.
    blocker_map = {
        "4AKE": ["interdomain_hinge_or_closure_evidence", "external_validation_if_available"],
        "1UBQ": ["cross_target_validation", "DCA_background_enrichment"],
        "1CLL": ["interdomain_hinge_evidence_if_claim_needed"],
        "p53_TAD_MDM2": ["external_couplings_if_available", "raw_local_external_annotation_files_if_needed"],
        "KcsA": ["external_couplings_if_available", "pore_filter_coupling_support", "expanded_biological_assembly_coordinate_or_author_defined_assembly_model_for_tetramer_claim"],
        "XCL1_lymphotactin": ["state_specific_external_couplings_or_constraints_if_available", "condition_labels_if_available"],
    }
    blockers = sorted(set(blocker_map.get(target, []) + missing))
    readiness = "not_ready_for_MD_external_constraints_or_target_specific_evidence_missing"
    if not blockers:
        readiness = "candidate_ready_for_deeper_non_claim_readout_not_MD_by_default"
    return {
        "target": target,
        "pressure_type": row.get("pressure_type"),
        "role_detected": row.get("role_detected") is True,
        "selected_signal": row.get("selected_signal"),
        "operators": operators,
        "missing_next_evidence": missing,
        "md_readiness": readiness,
        "md_blockers": blockers,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "claim_allowed": bool(claim_allowed),
    }


def build_v29(v25: dict[str, Any], v28: dict[str, Any]) -> dict[str, Any]:
    failed: list[str] = []

    if v25.get("sprint_status") != "V25_FAST_MECHANISM_EVIDENCE_SPRINT_LOCKED":
        failed.append("V25_fast_mechanism_sprint_locked")
    if v28.get("test_status") != "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST_PASSED_CLAIM_DISABLED":
        failed.append("V28_xcl1_state_condition_contrast_passed")

    # Boundary checks: no claim upgrade and no simulation in the summary layer.
    for label, cert in [("V25", v25), ("V28", v28)]:
        if cert.get("claim_allowed") is not False:
            failed.append(f"{label}_claim_disabled")
        if cert.get("new_md_executed") is not False:
            failed.append(f"{label}_no_new_md")
        if cert.get("positive_folding_evidence_found") not in (False, None):
            # V25 stores positive_folding_evidence_targets rather than the flag.
            failed.append(f"{label}_no_positive_folding_evidence")
        if cert.get("fixed_residue_cutoff_used") is not False and cert.get("fixed_residue_cutoff_used") is not None:
            failed.append(f"{label}_fixed_residue_cutoff_retired")
        if cert.get("native_metrics_used_for_selection") is not False and cert.get("native_metrics_used_for_selection") is not None:
            failed.append(f"{label}_native_metric_selection_forbidden")

    rows = _operator_targets(v25)
    if len(rows) < 6:
        failed.append("six_target_operator_table_present")

    operator_counts: Counter[str] = Counter()
    target_count = len(rows)
    for row in rows:
        operator_counts.update([str(x) for x in _as_list(row.get("mechanism_operator_used"))])

    repeated_core = [op for op, n in operator_counts.items() if n >= 3]
    universal_core_present = all(operator_counts.get(op, 0) == target_count for op in CORE_OPERATORS) if target_count else False
    if not universal_core_present:
        failed.append("core_operators_present_across_all_targets")

    xcl1_row = _target_row(rows, "XCL1_lymphotactin")
    xcl1_state_operator_preserved = bool(
        "state_separation" in _as_list(xcl1_row.get("mechanism_operator_used"))
        and v28.get("state_condition_contrast_preserved") is True
        and v28.get("mixed_state_pollution") is False
        and v28.get("single_fold_claim_made") is False
        and v28.get("fold_switch_claim_made") is False
    )
    if not xcl1_state_operator_preserved:
        failed.append("xcl1_state_condition_separation_operator_preserved")

    readiness_rows = [_row_md_readiness(row) for row in rows]
    md_ready_targets = [row["target"] for row in readiness_rows if row["md_readiness"].startswith("candidate_ready")]

    pressure_evidence_targets = [str(x) for x in _as_list(v25.get("positive_pressure_evidence_targets"))]
    folding_targets = [str(x) for x in _as_list(v25.get("positive_folding_evidence_targets"))]
    if folding_targets:
        failed.append("positive_folding_evidence_targets_empty")

    md_decision = {
        "kind": "V29_MD_READINESS_DECISION_v0",
        "decision_status": "V29_NO_NEW_MD_READY_EXTERNAL_CONSTRAINTS_REQUIRED",
        "md_ready_targets": md_ready_targets,
        "selected_next_panel": "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT",
        "selected_next_focus": "XCL1_lymphotactin_and_KcsA_constraints_before_any_MD",
        "reason": (
            "The operator grammar is coherent, but external/state-specific constraints remain missing. "
            "XCL1 needs state-specific external constraints; KcsA needs coupling/interface support and expanded assembly coordinates for any tetramer claim; "
            "p53 remains pressure-context evidence without coupling/MSA evidence."
        ),
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "membrane_MD_recommended": False,
        "parallel_MD_paths_allowed": False,
        "required_before_any_MD": sorted(set(
            blocker
            for row in readiness_rows
            for blocker in row.get("md_blockers", [])
            if blocker
        )),
        "claim_allowed": False,
        "positive_folding_evidence_targets": [],
    }

    summary = {
        "kind": "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_v0",
        "target_count": target_count,
        "targets": rows,
        "operator_counts": dict(sorted(operator_counts.items())),
        "repeated_mechanism_operators": sorted(repeated_core),
        "core_operator_policy": "role_detect_evidence_assign_pollution_guard_must_span_all_targets",
        "core_operators_across_all_targets": CORE_OPERATORS,
        "universal_core_operator_present": universal_core_present,
        "pressure_operator_family": PRESSURE_OPERATORS,
        "xcl1_state_condition_operator_preserved": xcl1_state_operator_preserved,
        "operator_interpretation": (
            "The shared mechanism is not a single threshold or a single fold. It is an operator grammar: "
            "role detection, context detection, evidence assignment, guard/pollution control, and target-specific pressure operators."
        ),
        "claim_allowed": False,
        "positive_folding_evidence_targets": [],
    }

    status = "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_LOCKED" if not failed else "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_BLOCKED"
    return {
        "kind": "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_DECISION_v0",
        "run_mode": "operator_panel_summary_and_md_readiness_no_simulation_no_threshold_tuning",
        "summary_status": status,
        "source_v25_status": v25.get("sprint_status"),
        "source_v28_status": v28.get("test_status"),
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "positive_pressure_evidence_targets": pressure_evidence_targets,
        "positive_folding_evidence_targets": [],
        "operator_panel_target_count": target_count,
        "repeated_mechanism_operators": summary["repeated_mechanism_operators"],
        "core_operators_across_all_targets": CORE_OPERATORS,
        "universal_core_operator_present": universal_core_present,
        "xcl1_state_condition_operator_preserved": xcl1_state_operator_preserved,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "md_readiness_decision_status": md_decision["decision_status"],
        "selected_next_panel": md_decision["selected_next_panel"],
        "summary_failed_checks": failed,
        "mechanism_operator_panel_summary": summary,
        "md_readiness_decision": md_decision,
        "md_readiness_rows": readiness_rows,
        "locked_interpretation": (
            "V29 summarizes the current protein-Esperanto operator panel across six targets and makes an explicit MD-readiness decision. "
            "The repeated mechanism is an operator grammar, not a universal folding claim. No MD is recommended until external/state-specific constraints and coupling evidence are acquired."
        ),
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V29 Mechanism Operator Panel Summary and MD Readiness Decision",
        "",
        f"Status: `{cert.get('summary_status')}`",
        f"Positive pressure evidence targets: `{cert.get('positive_pressure_evidence_targets')}`",
        f"Positive folding evidence targets: `{cert.get('positive_folding_evidence_targets')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        f"New MD recommended: `{cert.get('new_MD_recommended')}`",
        "",
        "## Locked interpretation",
        cert.get("locked_interpretation", ""),
        "",
        "## Repeated operators",
        ", ".join(cert.get("repeated_mechanism_operators", [])),
        "",
        "## MD readiness decision",
        f"Decision: `{cert.get('md_readiness_decision_status')}`",
        f"Selected next panel: `{cert.get('selected_next_panel')}`",
        f"Required before any MD: `{cert.get('md_readiness_decision', {}).get('required_before_any_MD')}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_rows_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "target", "pressure_type", "role_detected", "selected_signal", "operators", "missing_next_evidence", "md_readiness", "md_blockers", "new_MD_allowed", "new_MD_recommended", "claim_allowed",
        ])
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["operators"] = ";".join(str(x) for x in row.get("operators", []))
            out["missing_next_evidence"] = ";".join(str(x) for x in row.get("missing_next_evidence", []))
            out["md_blockers"] = ";".join(str(x) for x in row.get("md_blockers", []))
            writer.writerow(out)


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v29_mechanism_operator_panel_summary_and_md_readiness_certificate.json",
        "operator_panel": out_dir / "v29_mechanism_operator_panel_summary.json",
        "md_decision": out_dir / "v29_md_readiness_decision.json",
        "readiness_rows": out_dir / "v29_md_readiness_rows.json",
        "readiness_csv": out_dir / "v29_md_readiness_rows.csv",
        "report": out_dir / "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_DECISION_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["operator_panel"].write_text(json.dumps(cert["mechanism_operator_panel_summary"], indent=2, sort_keys=True), encoding="utf-8")
    paths["md_decision"].write_text(json.dumps(cert["md_readiness_decision"], indent=2, sort_keys=True), encoding="utf-8")
    paths["readiness_rows"].write_text(json.dumps(cert["md_readiness_rows"], indent=2, sort_keys=True), encoding="utf-8")
    _write_rows_csv(paths["readiness_csv"], cert["md_readiness_rows"])
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v25-cert", type=Path, default=DEFAULT_V25_CERT)
    parser.add_argument("--v28-cert", type=Path, default=DEFAULT_V28_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    v25 = _read_json(args.v25_cert, "V25 fast mechanism sprint certificate")
    v28 = _read_json(args.v28_cert, "V28 XCL1 state-condition contrast certificate")
    cert = build_v29(v25, v28)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "summary_status": cert["summary_status"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "new_MD_allowed": cert["new_MD_allowed"],
        "new_MD_recommended": cert["new_MD_recommended"],
        "positive_pressure_evidence_targets": cert["positive_pressure_evidence_targets"],
        "positive_folding_evidence_targets": cert["positive_folding_evidence_targets"],
        "operator_panel_target_count": cert["operator_panel_target_count"],
        "universal_core_operator_present": cert["universal_core_operator_present"],
        "xcl1_state_condition_operator_preserved": cert["xcl1_state_condition_operator_preserved"],
        "md_readiness_decision_status": cert["md_readiness_decision_status"],
        "selected_next_panel": cert["selected_next_panel"],
        "certificate": str(paths["certificate"]),
        "operator_panel": str(paths["operator_panel"]),
        "md_decision": str(paths["md_decision"]),
        "readiness_rows": str(paths["readiness_rows"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
