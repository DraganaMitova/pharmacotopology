#!/usr/bin/env python3
from __future__ import annotations

"""V33 constraint-backed operator readout.

V32 opens the door only after a real external source import selects a target.
V33 is still not a folding solver: it reads the V32-selected target's imported
external constraint files, performs a small no-MD operator readout, and keeps all
claim gates closed.

Scientific boundary:
- no MD, no membrane MD, no threshold tuning
- no native/coordinate metric is used for selection
- coordinate-derived source files may support operator-role readout only
- claim_allowed remains false
"""

import argparse
import csv
import json
import os
from pathlib import Path
from typing import Any, Iterable

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V32_CERT = RUN_ROOT / "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT" / "v32_external_constraint_source_import_preflight_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V33_CONSTRAINT_BACKED_OPERATOR_READOUT"

READY_V32_STATUS = "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED"
PASSED_STATUS = "V33_CONSTRAINT_BACKED_OPERATOR_READOUT_PASSED_CLAIM_DISABLED"
ABSTAIN_STATUS = "V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED"
BLOCKED_STATUS = "V33_BLOCKED_PRECONDITION_OR_PROVENANCE_FAILURE_CLAIM_DISABLED"


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT)).replace("\\", "/")
    except ValueError:
        return str(path).replace("\\", "/")


def _norm_rel(path_text: str) -> str:
    path_text = path_text.replace("\\", "/")
    p = Path(path_text)
    if p.is_absolute():
        return _rel(p)
    return path_text.lstrip("./")


def _load_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return [dict(row) for row in reader]


def _float_or_none(value: Any) -> float | None:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except Exception:
        return None


def _unique_nonempty(rows: Iterable[dict[str, str]], key: str) -> list[str]:
    return sorted({str(row.get(key, "")).strip() for row in rows if str(row.get(key, "")).strip()})


def _distance_summary(rows: list[dict[str, str]]) -> dict[str, Any]:
    distances = [_float_or_none(row.get("min_distance_angstrom")) for row in rows]
    values = sorted(d for d in distances if d is not None)
    if not values:
        return {"distance_values_present": False, "count": 0}
    mid = len(values) // 2
    median = values[mid] if len(values) % 2 else (values[mid - 1] + values[mid]) / 2.0
    return {
        "distance_values_present": True,
        "count": len(values),
        "min_distance_angstrom": round(values[0], 3),
        "median_distance_angstrom": round(median, 3),
        "max_distance_angstrom": round(values[-1], 3),
        "policy": "descriptive_only_no_threshold_selection",
    }


def _read_target_rows(v32: dict[str, Any], target: str) -> list[dict[str, Any]]:
    for target_row in _as_list(v32.get("target_rows")):
        if isinstance(target_row, dict) and target_row.get("target") == target:
            return [row for row in _as_list(target_row.get("valid_real_external_constraint_rows")) if isinstance(row, dict)]
    return []


def _source_rows_by_type(rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        out.setdefault(str(row.get("evidence_type", "")), []).append(row)
    return out


def _source_files_exist(rows: list[dict[str, Any]]) -> tuple[bool, list[str]]:
    missing: list[str] = []
    for row in rows:
        rel_path = _norm_rel(str(row.get("file_path", "")))
        if not rel_path or not (REPO_ROOT / rel_path).exists():
            missing.append(rel_path or "<missing_file_path>")
    return not missing, missing


def _kcsa_readout(source_rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_type = _source_rows_by_type(source_rows)
    pore_sources = by_type.get("pore_filter_coupling", [])
    interface_sources = by_type.get("assembly_interface_constraint", [])

    pore_rows: list[dict[str, str]] = []
    interface_rows: list[dict[str, str]] = []
    for row in pore_sources:
        pore_rows.extend(_load_csv(REPO_ROOT / _norm_rel(str(row.get("file_path", "")))))
    for row in interface_sources:
        interface_rows.extend(_load_csv(REPO_ROOT / _norm_rel(str(row.get("file_path", "")))))

    pore_chains = _unique_nonempty(pore_rows, "chain")
    pore_residue_numbers = _unique_nonempty(pore_rows, "residue_number") or _unique_nonempty(pore_rows, "residue_i")
    motif_labels = _unique_nonempty(pore_rows, "filter_motif")
    ion_names = _unique_nonempty(pore_rows, "ion_name")
    pore_constraint_classes = _unique_nonempty(pore_rows, "constraint_class")

    interface_chain_pairs = sorted({
        f"{str(row.get('chain_a', '')).strip()}-{str(row.get('chain_b', '')).strip()}"
        for row in interface_rows
        if str(row.get("chain_a", "")).strip() and str(row.get("chain_b", "")).strip()
    })
    interface_constraint_classes = _unique_nonempty(interface_rows, "constraint_class")

    operator_buckets: list[str] = []
    missing: list[str] = []
    if pore_rows:
        operator_buckets.append("pore_filter_operator")
    else:
        missing.append("pore_filter_external_constraint_rows")
    if interface_rows:
        operator_buckets.append("assembly_interface_operator")
    else:
        missing.append("assembly_interface_external_constraint_rows")
    if pore_rows and interface_rows:
        operator_buckets.append("membrane_pore_oligomer_constraint_context")

    readout = {
        "target": "KcsA",
        "readout_type": "constraint_backed_operator_readout_no_MD",
        "source_boundary": "external_coordinate_derived_source_import_readout_only_not_de_novo_folding_prediction",
        "pore_filter_source_file_count": len(pore_sources),
        "assembly_interface_source_file_count": len(interface_sources),
        "pore_filter_row_count": len(pore_rows),
        "assembly_interface_row_count": len(interface_rows),
        "pore_filter_chains": pore_chains,
        "pore_filter_residue_numbers_or_indices": pore_residue_numbers,
        "pore_filter_motif_labels": motif_labels,
        "pore_filter_ion_names": ion_names,
        "pore_filter_constraint_classes": pore_constraint_classes,
        "pore_filter_distance_summary": _distance_summary(pore_rows),
        "assembly_interface_chain_pair_count": len(interface_chain_pairs),
        "assembly_interface_chain_pairs_sample": interface_chain_pairs[:25],
        "assembly_interface_constraint_classes": interface_constraint_classes,
        "assembly_interface_distance_summary": _distance_summary(interface_rows),
        "operator_buckets_assigned": operator_buckets,
        "missing_operator_evidence": missing,
        "constraint_backed_operator_readout_found": not missing,
        "selection_threshold_used": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "claim_allowed": False,
        "new_MD_allowed": False,
    }
    return readout


def build_v33(v32: dict[str, Any]) -> dict[str, Any]:
    failed: list[str] = []
    if v32.get("preflight_status") != READY_V32_STATUS:
        failed.append("V32_status_not_ready_for_V33")
    for key in ["claim_allowed", "new_md_executed", "membrane_md_executed", "new_MD_allowed", "new_MD_recommended", "fixed_residue_cutoff_used", "native_metrics_used_for_selection"]:
        value = v32.get(key)
        if value not in (False, None):
            failed.append(f"V32_{key}_must_be_false")
    if v32.get("provenance_clean") is not True:
        failed.append("V32_provenance_not_clean")

    selected = v32.get("selected_V33_target")
    panel = v32.get("selected_V33_panel")
    if panel != "V33_CONSTRAINT_BACKED_OPERATOR_READOUT":
        failed.append("V32_selected_panel_not_V33_constraint_backed_operator_readout")
    if not selected:
        failed.append("V32_selected_V33_target_missing")

    source_rows = _read_target_rows(v32, str(selected)) if selected else []
    files_ok, missing_files = _source_files_exist(source_rows)
    if not files_ok:
        failed.append("selected_target_external_constraint_file_missing")

    readout: dict[str, Any] | None = None
    missing_operator_evidence: list[str] = []
    if not failed and selected == "KcsA":
        readout = _kcsa_readout(source_rows)
        missing_operator_evidence = _as_list(readout.get("missing_operator_evidence"))
    elif not failed:
        failed.append("selected_target_not_supported_by_V33_v0")

    if failed:
        status = BLOCKED_STATUS
    elif missing_operator_evidence:
        status = ABSTAIN_STATUS
    else:
        status = PASSED_STATUS

    available_evidence = []
    if readout and readout.get("pore_filter_row_count", 0):
        available_evidence.append("pore_filter_external_constraints")
    if readout and readout.get("assembly_interface_row_count", 0):
        available_evidence.append("assembly_interface_external_constraints")
    if readout and readout.get("constraint_backed_operator_readout_found"):
        available_evidence.append("constraint_backed_kcsa_operator_buckets")

    return {
        "kind": "V33_CONSTRAINT_BACKED_OPERATOR_READOUT_v0",
        "run_mode": "constraint_backed_operator_readout_no_MD_no_threshold_tuning_claim_disabled",
        "source_v32_status": v32.get("preflight_status"),
        "readout_status": status,
        "selected_V33_target": selected,
        "selected_V33_panel": panel,
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "provenance_clean": not failed and v32.get("provenance_clean") is True,
        "failed_checks": sorted(set(failed)),
        "missing_files": missing_files,
        "selected_target_import_rows": source_rows,
        "available_evidence": available_evidence,
        "missing_operator_evidence": missing_operator_evidence,
        "constraint_backed_operator_readout": readout,
        "operator_readout_found": bool(readout and readout.get("constraint_backed_operator_readout_found") and not failed),
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
        "locked_interpretation": (
            "V33 is a claim-disabled operator readout over V32-imported real external constraint files. "
            "For KcsA it can read pore/filter and assembly/interface buckets, but it is not a de novo folding prediction, not a whole-channel fold claim, and not a universal protein-folding solution."
        ),
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    readout = cert.get("constraint_backed_operator_readout") or {}
    lines = [
        "# V33 Constraint-Backed Operator Readout",
        "",
        f"Status: `{cert.get('readout_status')}`",
        f"Selected target: `{cert.get('selected_V33_target')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD allowed: `{cert.get('new_MD_allowed')}`",
        f"Positive folding evidence found: `{cert.get('positive_folding_evidence_found')}`",
        "",
        "## Locked interpretation",
        str(cert.get("locked_interpretation", "")),
        "",
        "## Readout summary",
        f"Pore/filter rows: `{readout.get('pore_filter_row_count')}`",
        f"Assembly/interface rows: `{readout.get('assembly_interface_row_count')}`",
        f"Operator buckets: `{readout.get('operator_buckets_assigned')}`",
        f"Missing operator evidence: `{cert.get('missing_operator_evidence')}`",
        "",
        "## Failed checks",
        json.dumps(cert.get("failed_checks", []), indent=2, sort_keys=True),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v33_constraint_backed_operator_readout_certificate.json"
    readout_path = out_dir / "v33_constraint_backed_operator_readout.json"
    report_path = out_dir / "V33_CONSTRAINT_BACKED_OPERATOR_READOUT_REPORT.md"
    decision_path = out_dir / "v33_next_decision.json"

    decision = {
        "kind": "V33_NEXT_DECISION_v0",
        "decision_status": cert.get("readout_status"),
        "selected_target": cert.get("selected_V33_target"),
        "next_action": (
            "run_V33_negative_controls_before_any_claim_or_MD"
            if cert.get("readout_status") == PASSED_STATUS
            else "fix_V32_source_import_or_operator_buckets_then_rerun_V33"
        ),
        "claim_allowed": False,
        "new_MD_allowed": False,
        "folding_problem_solved": False,
    }
    cert_with_artifacts = {
        **cert,
        "next_decision": decision,
        "artifacts": {
            "certificate": str(cert_path),
            "readout": str(readout_path),
            "report": str(report_path),
            "decision": str(decision_path),
        },
    }
    _write_json(cert_path, cert_with_artifacts)
    _write_json(readout_path, cert.get("constraint_backed_operator_readout"))
    _write_json(decision_path, decision)
    _write_report(report_path, cert_with_artifacts)
    return {"certificate": cert_path, "readout": readout_path, "report": report_path, "decision": decision_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V33 claim-disabled constraint-backed operator readout.")
    parser.add_argument("--v32-certificate", type=Path, default=DEFAULT_V32_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    v32 = _read_json(args.v32_certificate, "V32 external constraint source import preflight certificate")
    cert = build_v33(v32)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert.get("kind"),
        "readout_status": cert.get("readout_status"),
        "selected_V33_target": cert.get("selected_V33_target"),
        "operator_readout_found": cert.get("operator_readout_found"),
        "positive_folding_evidence_found": cert.get("positive_folding_evidence_found"),
        "folding_problem_solved": cert.get("folding_problem_solved"),
        "claim_allowed": cert.get("claim_allowed"),
        "new_MD_allowed": cert.get("new_MD_allowed"),
        "failed_checks": cert.get("failed_checks"),
        "missing_operator_evidence": cert.get("missing_operator_evidence"),
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
        "next_action": "run_V33_negative_controls_before_any_claim_or_MD" if cert.get("readout_status") == PASSED_STATUS else "fix_source_import_or_operator_buckets",
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
