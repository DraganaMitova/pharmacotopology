#!/usr/bin/env python3
from __future__ import annotations

"""V30 external constraint + coupling acquisition sprint.

This is a practical acquisition/readiness layer after V29. It scans local files
for target-specific external constraints/couplings, preserves the no-MD gate,
and selects the next acquisition panel. It does not synthesize couplings, does
not tune thresholds, and does not upgrade any pressure-context evidence into a
folding claim.
"""

import argparse
import csv
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V29_CERT = RUN_ROOT / "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_DECISION" / "v29_mechanism_operator_panel_summary_and_md_readiness_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT"

TARGET_ALIASES: dict[str, list[str]] = {
    "p53_TAD_MDM2": ["p53", "tp53", "mdm2", "1ycr", "dp00086", "disprot"],
    "KcsA": ["kcsa", "1k4c", "1bl8", "opm", "potassium", "tvgyg", "k_channel"],
    "XCL1_lymphotactin": ["xcl1", "lymphotactin", "2hdm", "2jp1", "2n54", "metamorphic"],
    "4AKE": ["4ake", "adenylate", "p69441"],
    "1UBQ": ["1ubq", "ubiquitin", "1d3z"],
    "1CLL": ["1cll", "calmodulin", "1cfc"],
}

CONSTRAINT_KEYWORDS = [
    "coupling", "couplings", "dca", "msa", "a3m", "sto", "stockholm", "alignment",
    "constraint", "constraints", "pfam", "evcoupling", "gremlin", "plmc",
]
ANNOTATION_KEYWORDS = [
    "annotation", "annotations", "condition", "state", "assembly", "biological",
    "opm", "disprot", "membrane", "topology", "tetramer", "oligomer",
]
SCAN_DIRS = ["data", "external_msa", "external_msa_free_predictors", "first_contact_clean_pharmacotopology_layer_run"]
SKIP_DIR_NAMES = {".git", ".pytest_cache", "__pycache__"}
MAX_SCAN_FILE_SIZE = 2_500_000


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _rel(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _looks_like(path_text: str, keywords: list[str]) -> bool:
    low = path_text.lower()
    return any(k in low for k in keywords)


def _file_mentions_target(path: Path, root: Path, target: str) -> bool:
    aliases = TARGET_ALIASES.get(target, [])
    rel = _rel(path, root).lower()
    if any(alias.lower() in rel for alias in aliases):
        return True
    # Carefully inspect small text/JSON/FASTA-ish files only; never parse big PDBs
    # or binary files for this acquisition layer.
    if path.suffix.lower() in {".json", ".csv", ".tsv", ".txt", ".md", ".fasta", ".fa", ".sto", ".a3m", ".aln"}:
        try:
            if path.stat().st_size <= MAX_SCAN_FILE_SIZE:
                text = path.read_text(encoding="utf-8", errors="ignore").lower()
                return any(alias.lower() in text for alias in aliases)
        except OSError:
            return False
    return False


def scan_local_constraint_files(root: Path) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {
        target: {
            "target": target,
            "candidate_coupling_or_constraint_files": [],
            "candidate_annotation_files": [],
            "raw_local_external_annotation_files": [],
        }
        for target in TARGET_ALIASES
    }
    candidates: list[Path] = []
    for dirname in SCAN_DIRS:
        base = root / dirname
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if any(part in SKIP_DIR_NAMES for part in path.parts):
                continue
            if not path.is_file():
                continue
            rel_text = _rel(path, root)
            if _looks_like(rel_text, CONSTRAINT_KEYWORDS + ANNOTATION_KEYWORDS):
                candidates.append(path)

    for target in TARGET_ALIASES:
        for path in candidates:
            if not _file_mentions_target(path, root, target):
                continue
            rel = _rel(path, root)
            if _looks_like(rel, CONSTRAINT_KEYWORDS):
                rows[target]["candidate_coupling_or_constraint_files"].append(rel)
            if _looks_like(rel, ANNOTATION_KEYWORDS):
                rows[target]["candidate_annotation_files"].append(rel)
            if _looks_like(rel, ANNOTATION_KEYWORDS) and not rel.lower().endswith((".pdb", ".cif", ".mmcif")):
                rows[target]["raw_local_external_annotation_files"].append(rel)
        for key in ["candidate_coupling_or_constraint_files", "candidate_annotation_files", "raw_local_external_annotation_files"]:
            rows[target][key] = sorted(set(rows[target][key]))
    return rows


def _md_row_map(v29: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {str(row.get("target")): dict(row) for row in _as_list(v29.get("md_readiness_rows")) if isinstance(row, dict)}


def _operator_row_map(v29: dict[str, Any]) -> dict[str, dict[str, Any]]:
    panel = v29.get("mechanism_operator_panel_summary") if isinstance(v29.get("mechanism_operator_panel_summary"), dict) else {}
    return {str(row.get("target")): dict(row) for row in _as_list(panel.get("targets")) if isinstance(row, dict)}


def _pressure_targets(v29: dict[str, Any]) -> list[str]:
    active = ["XCL1_lymphotactin", "KcsA", "p53_TAD_MDM2", "4AKE", "1UBQ", "1CLL"]
    seen = set()
    out: list[str] = []
    for item in active + [str(x) for x in _as_list(v29.get("positive_pressure_evidence_targets"))]:
        if item and item not in seen:
            seen.add(item)
            out.append(item)
    return out


def _target_priority(target: str) -> int:
    # V29 selected XCL1 and KcsA constraints before MD. Keep p53 third because
    # its external contrast is already locked and couplings are optional/future.
    return {
        "XCL1_lymphotactin": 1,
        "KcsA": 2,
        "p53_TAD_MDM2": 3,
        "4AKE": 4,
        "1UBQ": 5,
        "1CLL": 6,
    }.get(target, 99)


def _target_acquisition_row(target: str, scan: dict[str, Any], md_row: dict[str, Any], op_row: dict[str, Any]) -> dict[str, Any]:
    coupling_files = list(scan.get("candidate_coupling_or_constraint_files", []))
    annotation_files = list(scan.get("candidate_annotation_files", []))
    raw_annotation_files = list(scan.get("raw_local_external_annotation_files", []))
    blockers = [str(x) for x in _as_list(md_row.get("md_blockers"))]
    missing = sorted(set(blockers + [str(x) for x in _as_list(op_row.get("missing_next_evidence"))]))

    state_condition_context = target == "XCL1_lymphotactin"
    kcsa_membrane_context = target == "KcsA"
    p53_external_contrast = target == "p53_TAD_MDM2"

    # Coupling files are useful only if they are target-specific; do not treat
    # legacy generic files as evidence for V16/V18+ pressure targets.
    coupling_present = bool(coupling_files)
    annotation_present = bool(annotation_files) or state_condition_context or kcsa_membrane_context or p53_external_contrast

    acquisition_status = "constraint_or_coupling_missing_acquisition_needed"
    if coupling_present:
        acquisition_status = "local_constraint_or_coupling_files_present_ready_for_next_readout_preflight"
    elif annotation_present:
        acquisition_status = "annotation_context_locked_coupling_missing"

    if target == "XCL1_lymphotactin":
        needed_before_md = [
            "state_specific_external_couplings_or_constraints_if_available",
            "independent_condition_or_oligomerization_context_if_available",
            "mixed_state_leakage_guard_preserved",
        ]
        recommended_next_action = "acquire_or_import_state_specific_external_constraints_then_preflight_no_MD"
    elif target == "KcsA":
        needed_before_md = [
            "external_couplings_if_available",
            "pore_filter_coupling_support",
            "expanded_biological_assembly_coordinate_or_author_defined_assembly_model_for_tetramer_claim",
        ]
        recommended_next_action = "acquire_or_import_KcsA_external_couplings_and_assembly_interface_constraints_no_MD"
    elif target == "p53_TAD_MDM2":
        needed_before_md = [
            "external_couplings_if_available",
            "raw_local_external_annotation_files_if_needed",
        ]
        recommended_next_action = "optional_external_coupling_or_raw_annotation_import_no_MD"
    else:
        needed_before_md = sorted(set(missing))
        recommended_next_action = "keep_locked_until_target_specific_external_constraints_are_available"

    return {
        "target": target,
        "pressure_type": op_row.get("pressure_type"),
        "role_detected": op_row.get("role_detected") is True,
        "selected_signal": op_row.get("selected_signal"),
        "mechanism_operator_used": _as_list(op_row.get("mechanism_operator_used")),
        "candidate_coupling_or_constraint_files": coupling_files,
        "candidate_annotation_files": annotation_files,
        "raw_local_external_annotation_files": raw_annotation_files,
        "local_constraints_or_couplings_present": coupling_present,
        "local_annotation_context_present_or_locked": annotation_present,
        "missing_external_constraints_or_couplings": [x for x in needed_before_md if x not in coupling_files],
        "acquisition_status": acquisition_status,
        "recommended_next_action": recommended_next_action,
        "ready_for_MD": False,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "claim_allowed": False,
    }


def build_v30(v29: dict[str, Any], repo_root: Path) -> dict[str, Any]:
    failed: list[str] = []
    if v29.get("summary_status") != "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_LOCKED":
        failed.append("V29_operator_panel_summary_locked")
    if v29.get("md_readiness_decision_status") != "V29_NO_NEW_MD_READY_EXTERNAL_CONSTRAINTS_REQUIRED":
        failed.append("V29_no_new_md_ready_external_constraints_required")
    for key in ["claim_allowed", "new_md_executed", "fixed_residue_cutoff_used", "native_metrics_used_for_selection"]:
        value = v29.get(key)
        if value not in (False, None):
            failed.append(f"V29_{key}_must_be_false")

    scan_rows = scan_local_constraint_files(repo_root)
    md_rows = _md_row_map(v29)
    op_rows = _operator_row_map(v29)
    targets = _pressure_targets(v29)
    acquisition_rows = [
        _target_acquisition_row(target, scan_rows.get(target, {}), md_rows.get(target, {}), op_rows.get(target, {}))
        for target in targets
    ]
    acquisition_rows = sorted(acquisition_rows, key=lambda row: _target_priority(str(row.get("target"))))

    coupling_ready_targets = [str(row.get("target")) for row in acquisition_rows if row.get("local_constraints_or_couplings_present") is True]
    annotation_locked_targets = [str(row.get("target")) for row in acquisition_rows if row.get("local_annotation_context_present_or_locked") is True]
    md_ready_targets: list[str] = []

    # The next acquisition sprint can batch XCL1 and KcsA because this is not MD;
    # both are the V29-selected focus before any simulation.
    selected_v31_targets = ["XCL1_lymphotactin", "KcsA"]
    selected_next_panel = "V31_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_AND_PREFLIGHT_SPRINT"
    if coupling_ready_targets:
        selected_next_panel = "V31_CONSTRAINT_BACKED_OPERATOR_READOUT_PREFLIGHT"

    manifest = {
        "kind": "V30_EXTERNAL_CONSTRAINT_MANIFEST_v0",
        "purpose": "acquire_or_import_external_constraints_and_couplings_before_any_MD",
        "targets": acquisition_rows,
        "claim_allowed": False,
        "new_md_executed": False,
        "new_MD_allowed": False,
        "positive_folding_evidence_targets": [],
    }
    preflight = {
        "kind": "V30_LOCAL_CONSTRAINT_PREFLIGHT_v0",
        "targets_with_annotation_context_present_or_locked": annotation_locked_targets,
        "targets_with_local_constraints_or_couplings": coupling_ready_targets,
        "targets_ready_for_MD": md_ready_targets,
        "all_targets_have_explicit_status": all(bool(row.get("acquisition_status")) for row in acquisition_rows),
        "claim_allowed": False,
        "new_md_executed": False,
    }
    coupling_scan = {
        "kind": "V30_COUPLING_AND_CONSTRAINT_AVAILABILITY_SCAN_v0",
        "scan_policy": "target_specific_local_files_only_no_synthetic_couplings_no_web_fetch_in_readout",
        "targets": {
            str(row.get("target")): {
                "candidate_coupling_or_constraint_files": row.get("candidate_coupling_or_constraint_files", []),
                "candidate_annotation_files": row.get("candidate_annotation_files", []),
                "local_constraints_or_couplings_present": row.get("local_constraints_or_couplings_present", False),
            }
            for row in acquisition_rows
        },
        "claim_allowed": False,
        "new_md_executed": False,
    }
    decision = {
        "kind": "V30_NEXT_ACQUISITION_DECISION_v0",
        "decision_status": "V30_EXTERNAL_CONSTRAINT_ACQUISITION_LOCKED_NO_MD_READY",
        "selected_next_panel": selected_next_panel,
        "selected_V31_targets": selected_v31_targets,
        "primary_focus": "XCL1_state_specific_constraints_and_KcsA_coupling_interface_support",
        "reason": (
            "V29 shows coherent mechanism operators but no target is MD-ready. "
            "XCL1 needs state-specific external constraints/condition context; KcsA needs coupling/interface support and expanded assembly coordinates. "
            "The next sprint should import or build external constraint sources, not run MD."
        ),
        "coupling_ready_targets": coupling_ready_targets,
        "md_ready_targets": md_ready_targets,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "membrane_MD_recommended": False,
        "parallel_MD_paths_allowed": False,
        "claim_allowed": False,
        "required_before_any_MD": sorted(set(
            item
            for row in acquisition_rows
            for item in row.get("missing_external_constraints_or_couplings", [])
        )),
    }

    if md_ready_targets:
        failed.append("no_targets_should_be_md_ready_in_v30")
    if v29.get("positive_folding_evidence_targets") not in ([], None):
        failed.append("positive_folding_evidence_targets_empty")

    status = "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT_LOCKED" if not failed else "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT_BLOCKED"
    return {
        "kind": "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT_v0",
        "run_mode": "external_constraint_coupling_acquisition_no_simulation_no_threshold_tuning",
        "sprint_status": status,
        "source_v29_status": v29.get("summary_status"),
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
        "positive_pressure_evidence_targets": _as_list(v29.get("positive_pressure_evidence_targets")),
        "positive_folding_evidence_targets": [],
        "constraint_or_coupling_ready_targets": coupling_ready_targets,
        "annotation_context_locked_targets": annotation_locked_targets,
        "md_ready_targets": md_ready_targets,
        "selected_next_panel": selected_next_panel,
        "selected_V31_targets": selected_v31_targets,
        "sprint_failed_checks": failed,
        "external_constraint_manifest": manifest,
        "local_constraint_preflight": preflight,
        "coupling_availability_scan": coupling_scan,
        "next_acquisition_decision": decision,
        "target_rows": acquisition_rows,
        "locked_interpretation": (
            "V30 converts the V29 MD-readiness block into a practical acquisition sprint: local target-specific couplings/constraints are scanned, missing constraint layers are made explicit, and XCL1/KcsA are selected for the next external-constraint import/preflight step. No MD or folding claim is allowed."
        ),
    }


def _write_rows_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "target", "pressure_type", "local_constraints_or_couplings_present", "local_annotation_context_present_or_locked",
        "acquisition_status", "recommended_next_action", "missing_external_constraints_or_couplings", "ready_for_MD", "new_MD_allowed", "claim_allowed",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            out = {key: row.get(key) for key in fields}
            out["missing_external_constraints_or_couplings"] = ";".join(str(x) for x in row.get("missing_external_constraints_or_couplings", []))
            writer.writerow(out)


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V30 External Constraint and Coupling Acquisition Sprint",
        "",
        f"Status: `{cert.get('sprint_status')}`",
        f"Selected next panel: `{cert.get('selected_next_panel')}`",
        f"Selected V31 targets: `{cert.get('selected_V31_targets')}`",
        f"Constraint/coupling-ready targets: `{cert.get('constraint_or_coupling_ready_targets')}`",
        f"MD-ready targets: `{cert.get('md_ready_targets')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD allowed: `{cert.get('new_MD_allowed')}`",
        "",
        "## Locked interpretation",
        cert.get("locked_interpretation", ""),
        "",
        "## Target rows",
    ]
    for row in cert.get("target_rows", []):
        lines.extend([
            f"### {row.get('target')}",
            f"Acquisition status: `{row.get('acquisition_status')}`",
            f"Coupling/constraint files: `{row.get('candidate_coupling_or_constraint_files')}`",
            f"Missing: `{row.get('missing_external_constraints_or_couplings')}`",
            f"Recommended next action: `{row.get('recommended_next_action')}`",
            "",
        ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v30_external_constraint_and_coupling_acquisition_sprint_certificate.json",
        "manifest": out_dir / "v30_external_constraint_manifest.json",
        "preflight": out_dir / "v30_local_constraint_preflight.json",
        "coupling_scan": out_dir / "v30_coupling_and_constraint_availability_scan.json",
        "decision": out_dir / "v30_next_acquisition_decision.json",
        "rows": out_dir / "v30_external_constraint_target_rows.json",
        "rows_csv": out_dir / "v30_external_constraint_target_rows.csv",
        "report": out_dir / "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["manifest"].write_text(json.dumps(cert["external_constraint_manifest"], indent=2, sort_keys=True), encoding="utf-8")
    paths["preflight"].write_text(json.dumps(cert["local_constraint_preflight"], indent=2, sort_keys=True), encoding="utf-8")
    paths["coupling_scan"].write_text(json.dumps(cert["coupling_availability_scan"], indent=2, sort_keys=True), encoding="utf-8")
    paths["decision"].write_text(json.dumps(cert["next_acquisition_decision"], indent=2, sort_keys=True), encoding="utf-8")
    paths["rows"].write_text(json.dumps(cert["target_rows"], indent=2, sort_keys=True), encoding="utf-8")
    _write_rows_csv(paths["rows_csv"], cert["target_rows"])
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v29-cert", type=Path, default=DEFAULT_V29_CERT)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    v29 = _read_json(args.v29_cert, "V29 mechanism operator panel certificate")
    cert = build_v30(v29, args.repo_root)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "sprint_status": cert["sprint_status"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "new_MD_allowed": cert["new_MD_allowed"],
        "new_MD_recommended": cert["new_MD_recommended"],
        "constraint_or_coupling_ready_targets": cert["constraint_or_coupling_ready_targets"],
        "md_ready_targets": cert["md_ready_targets"],
        "selected_next_panel": cert["selected_next_panel"],
        "selected_V31_targets": cert["selected_V31_targets"],
        "certificate": str(paths["certificate"]),
        "manifest": str(paths["manifest"]),
        "preflight": str(paths["preflight"]),
        "coupling_scan": str(paths["coupling_scan"]),
        "decision": str(paths["decision"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
