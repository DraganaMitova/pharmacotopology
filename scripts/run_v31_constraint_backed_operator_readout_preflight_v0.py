#!/usr/bin/env python3
from __future__ import annotations

"""V31 constraint-backed operator readout preflight.

This is a strict provenance gate after V30. V30 intentionally scanned broadly and
may list local generated reports next to possible external constraints. V31
separates those file classes before any V32 constraint-backed operator readout.

The key rule is simple: internal generated artifacts can support audit/history,
but they are never allowed as external constraint or coupling evidence.
"""

import argparse
import csv
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V30_CERT = RUN_ROOT / "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT" / "v30_external_constraint_and_coupling_acquisition_sprint_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V31_CONSTRAINT_BACKED_OPERATOR_READOUT_PREFLIGHT"

SELECTED_SCOPE_DEFAULT = ["XCL1_lymphotactin", "KcsA"]
TARGET_PRIORITY = {
    "XCL1_lymphotactin": 1,
    "KcsA": 2,
    "p53_TAD_MDM2": 3,
    "4AKE": 4,
    "1UBQ": 5,
    "1CLL": 6,
}

CONSTRAINT_KEYWORDS = [
    "coupling", "couplings", "constraint", "constraints", "dca", "plmc",
    "evcoupling", "gremlin", "contact", "contacts", "pair", "pairs",
]
ALIGNMENT_KEYWORDS = ["msa", "alignment", "stockholm", "pfam", "a3m", "sto"]
ANNOTATION_KEYWORDS = [
    "annotation", "condition", "state", "assembly", "biological", "opm",
    "disprot", "membrane", "topology", "tetramer", "oligomer", "pore",
]
STRUCTURE_SUFFIXES = {".pdb", ".cif", ".mmcif"}
TEXT_SUFFIXES = {".json", ".csv", ".tsv", ".txt", ".md", ".fasta", ".fa", ".sto", ".a3m", ".aln"}


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _low(path: str) -> str:
    return path.replace("\\", "/").lower()


def _has_any(text: str, words: list[str]) -> bool:
    return any(word in text for word in words)


def _suffix(path_text: str) -> str:
    p = Path(path_text)
    if p.suffix.lower() == ".gz" and p.name.lower().endswith(".pdb.gz"):
        return ".pdb"
    return p.suffix.lower()


def classify_candidate_file(path_text: str) -> dict[str, Any]:
    """Classify one candidate path by provenance and allowed use.

    The classification is path/provenance based by design. It does not infer
    truth from file names beyond conservative source-class separation.
    """
    low = _low(path_text)
    suffix = _suffix(path_text)

    if low.startswith("first_contact_clean_pharmacotopology_layer_run/") or low.startswith("data/locked_runtime_certificates/"):
        return {
            "path": path_text,
            "evidence_source_class": "generated_internal_report",
            "allowed_use": "allowed_for_audit_only",
            "counts_as_real_external_constraint": False,
            "counts_as_annotation_context": False,
            "exclusion_reason": "runtime_or_exported_internal_artifact_not_external_evidence",
        }

    if low.startswith("external_msa/") or (suffix in {".sto", ".a3m", ".aln"} and _has_any(low, ALIGNMENT_KEYWORDS)):
        return {
            "path": path_text,
            "evidence_source_class": "real_external_alignment_source",
            "allowed_use": "allowed_for_constraint_derivation_preflight_only",
            "counts_as_real_external_constraint": False,
            "counts_as_annotation_context": False,
            "exclusion_reason": None,
        }

    if low.startswith("external_msa_free_predictors/"):
        # Useful for audit/model comparison, but not external coupling evidence.
        return {
            "path": path_text,
            "evidence_source_class": "model_prediction_or_msa_free_artifact",
            "allowed_use": "allowed_for_audit_only",
            "counts_as_real_external_constraint": False,
            "counts_as_annotation_context": False,
            "exclusion_reason": "model_prediction_not_external_constraint_or_coupling",
        }

    if suffix in STRUCTURE_SUFFIXES:
        if low.startswith("data/independent_contact_sources/") or low.startswith("data/v16_pressure_targets/"):
            return {
                "path": path_text,
                "evidence_source_class": "external_structure_source",
                "allowed_use": "allowed_for_structure_context_or_validation_only",
                "counts_as_real_external_constraint": False,
                "counts_as_annotation_context": True,
                "exclusion_reason": None,
            }
        return {
            "path": path_text,
            "evidence_source_class": "structure_file_unknown_or_generated",
            "allowed_use": "allowed_for_audit_only",
            "counts_as_real_external_constraint": False,
            "counts_as_annotation_context": False,
            "exclusion_reason": "structure_file_is_not_constraint_or_coupling_evidence",
        }

    if low.startswith("data/") and _has_any(low, CONSTRAINT_KEYWORDS) and suffix in TEXT_SUFFIXES:
        return {
            "path": path_text,
            "evidence_source_class": "real_external_constraint_or_coupling",
            "allowed_use": "allowed_for_constraint_backed_operator_readout",
            "counts_as_real_external_constraint": True,
            "counts_as_annotation_context": False,
            "exclusion_reason": None,
        }

    if (low.startswith("data/") or low.startswith("external_annotations/") or low.startswith("annotations/")) and _has_any(low, ANNOTATION_KEYWORDS) and suffix in TEXT_SUFFIXES:
        return {
            "path": path_text,
            "evidence_source_class": "annotation_only_external_context",
            "allowed_use": "allowed_for_role_context_only",
            "counts_as_real_external_constraint": False,
            "counts_as_annotation_context": True,
            "exclusion_reason": None,
        }

    if _has_any(low, CONSTRAINT_KEYWORDS) and suffix in TEXT_SUFFIXES:
        return {
            "path": path_text,
            "evidence_source_class": "unverified_constraint_like_file",
            "allowed_use": "excluded",
            "counts_as_real_external_constraint": False,
            "counts_as_annotation_context": False,
            "exclusion_reason": "constraint_like_path_not_under_accepted_external_source_root",
        }

    if _has_any(low, ANNOTATION_KEYWORDS) and suffix in TEXT_SUFFIXES:
        return {
            "path": path_text,
            "evidence_source_class": "unverified_annotation_like_file",
            "allowed_use": "excluded",
            "counts_as_real_external_constraint": False,
            "counts_as_annotation_context": False,
            "exclusion_reason": "annotation_like_path_not_under_accepted_external_source_root",
        }

    return {
        "path": path_text,
        "evidence_source_class": "unusable_or_unclassified",
        "allowed_use": "excluded",
        "counts_as_real_external_constraint": False,
        "counts_as_annotation_context": False,
        "exclusion_reason": "not_recognized_as_external_constraint_coupling_or_annotation_context",
    }


def _target_priority(target: str) -> int:
    return TARGET_PRIORITY.get(target, 99)


def _class_counts(classified: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for item in classified:
        key = str(item.get("evidence_source_class"))
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def _allowed_counts(classified: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for item in classified:
        key = str(item.get("allowed_use"))
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def _target_v31_row(v30_row: dict[str, Any], selected_scope: list[str]) -> dict[str, Any]:
    target = str(v30_row.get("target"))
    raw_candidates = []
    for key in ["candidate_coupling_or_constraint_files", "candidate_annotation_files", "raw_local_external_annotation_files"]:
        for item in _as_list(v30_row.get(key)):
            if isinstance(item, str):
                raw_candidates.append(item)
    candidates = sorted(set(raw_candidates))
    classified = [classify_candidate_file(path) for path in candidates]
    real_external = [item["path"] for item in classified if item.get("counts_as_real_external_constraint") is True]
    annotation_only = [item["path"] for item in classified if item.get("allowed_use") in {"allowed_for_role_context_only", "allowed_for_structure_context_or_validation_only"}]
    generated_internal = [item["path"] for item in classified if item.get("evidence_source_class") == "generated_internal_report"]
    excluded = [item["path"] for item in classified if item.get("allowed_use") == "excluded"]
    derivation_only = [item["path"] for item in classified if item.get("allowed_use") == "allowed_for_constraint_derivation_preflight_only"]

    in_selected_v31_scope = target in selected_scope
    constraint_backed_allowed = bool(real_external)
    status = "real_external_constraint_available_for_v32_preflight" if constraint_backed_allowed else "no_real_external_constraint_for_constraint_backed_readout"
    if not candidates:
        status = "no_candidate_files_listed_by_v30"

    return {
        "target": target,
        "in_selected_V31_scope": in_selected_v31_scope,
        "pressure_type": v30_row.get("pressure_type"),
        "selected_signal": v30_row.get("selected_signal"),
        "v30_acquisition_status": v30_row.get("acquisition_status"),
        "candidate_file_count": len(candidates),
        "classified_candidate_files": classified,
        "evidence_source_class_counts": _class_counts(classified),
        "allowed_use_counts": _allowed_counts(classified),
        "real_external_constraint_or_coupling_files": real_external,
        "annotation_context_files": annotation_only,
        "constraint_derivation_only_files": derivation_only,
        "generated_internal_report_files": generated_internal,
        "excluded_files": excluded,
        "constraint_backed_operator_readout_allowed": constraint_backed_allowed,
        "allowed_use": "allowed_for_constraint_backed_operator_readout" if constraint_backed_allowed else "not_allowed_for_constraint_backed_operator_readout",
        "preflight_status": status,
        "ready_for_MD": False,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "claim_allowed": False,
    }


def _select_v32_target(rows: list[dict[str, Any]], selected_scope: list[str]) -> str | None:
    eligible = [
        row for row in rows
        if row.get("in_selected_V31_scope") is True and row.get("constraint_backed_operator_readout_allowed") is True
    ]
    if not eligible:
        return None
    eligible = sorted(
        eligible,
        key=lambda row: (
            -len(_as_list(row.get("real_external_constraint_or_coupling_files"))),
            _target_priority(str(row.get("target"))),
        ),
    )
    return str(eligible[0].get("target"))


def build_v31(v30: dict[str, Any]) -> dict[str, Any]:
    failed: list[str] = []
    if v30.get("sprint_status") != "V30_EXTERNAL_CONSTRAINT_AND_COUPLING_ACQUISITION_SPRINT_LOCKED":
        failed.append("V30_external_constraint_acquisition_locked")
    for key in ["claim_allowed", "new_md_executed", "new_MD_allowed", "new_MD_recommended", "native_metrics_used_for_selection"]:
        value = v30.get(key)
        if value not in (False, None):
            failed.append(f"V30_{key}_must_be_false")

    selected_scope = [str(x) for x in _as_list(v30.get("selected_V31_targets"))] or list(SELECTED_SCOPE_DEFAULT)
    v30_rows = [row for row in _as_list(v30.get("target_rows")) if isinstance(row, dict)]
    target_rows = sorted([_target_v31_row(row, selected_scope) for row in v30_rows], key=lambda row: _target_priority(str(row.get("target"))))

    real_targets_all = [str(row.get("target")) for row in target_rows if row.get("constraint_backed_operator_readout_allowed") is True]
    real_targets_in_scope = [str(row.get("target")) for row in target_rows if row.get("in_selected_V31_scope") is True and row.get("constraint_backed_operator_readout_allowed") is True]
    selected_v32_target = _select_v32_target(target_rows, selected_scope)

    generated_misused = []
    annotation_misused = []
    for row in target_rows:
        for item in _as_list(row.get("classified_candidate_files")):
            if not isinstance(item, dict):
                continue
            if item.get("evidence_source_class") == "generated_internal_report" and item.get("counts_as_real_external_constraint") is True:
                generated_misused.append(str(item.get("path")))
            if item.get("allowed_use") in {"allowed_for_role_context_only", "allowed_for_structure_context_or_validation_only"} and item.get("counts_as_real_external_constraint") is True:
                annotation_misused.append(str(item.get("path")))

    provenance_clean = not generated_misused and not annotation_misused
    if not provenance_clean:
        failed.append("generated_or_annotation_files_misused_as_external_constraints")

    if failed:
        preflight_status = "V31_CONSTRAINT_BACKED_OPERATOR_READOUT_PREFLIGHT_BLOCKED"
    elif selected_v32_target:
        preflight_status = "V31_CONSTRAINT_BACKED_OPERATOR_READOUT_PREFLIGHT_PASSED_CLAIM_DISABLED"
    else:
        preflight_status = "V31_CLEAN_ABSTAIN_NO_REAL_EXTERNAL_CONSTRAINTS_FOR_SELECTED_TARGETS"

    selected_v32_panel = "V32_FIRST_CONSTRAINT_BACKED_OPERATOR_READOUT" if selected_v32_target else None
    next_action = (
        "run_V32_constraint_backed_operator_readout_for_selected_target_no_MD"
        if selected_v32_target
        else "acquire_or_import_real_external_constraints_for_XCL1_or_KcsA_before_constraint_backed_readout"
    )

    allowed_use_matrix = {
        "real_external_constraint_or_coupling": "allowed_for_constraint_backed_operator_readout",
        "real_external_alignment_source": "allowed_for_constraint_derivation_preflight_only",
        "annotation_only_external_context": "allowed_for_role_context_only",
        "external_structure_source": "allowed_for_structure_context_or_validation_only",
        "generated_internal_report": "allowed_for_audit_only",
        "model_prediction_or_msa_free_artifact": "allowed_for_audit_only",
        "unverified_constraint_like_file": "excluded",
        "unverified_annotation_like_file": "excluded",
        "unusable_or_unclassified": "excluded",
    }

    decision = {
        "kind": "V31_NEXT_OPERATOR_READOUT_DECISION_v0",
        "decision_status": preflight_status,
        "selected_V32_target": selected_v32_target,
        "selected_V32_panel": selected_v32_panel,
        "selected_V31_targets": selected_scope,
        "real_external_constraint_targets_all": real_targets_all,
        "real_external_constraint_targets_in_selected_scope": real_targets_in_scope,
        "next_action": next_action,
        "reason": (
            "V31 separates real external constraints/couplings from annotation context and generated internal reports. "
            "Only selected-scope targets with real external constraint/coupling files can enter V32 as constraint-backed readouts."
        ),
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "claim_allowed": False,
    }

    return {
        "kind": "V31_CONSTRAINT_BACKED_OPERATOR_READOUT_PREFLIGHT_v0",
        "run_mode": "provenance_classification_no_simulation_no_threshold_tuning",
        "preflight_status": preflight_status,
        "source_v30_status": v30.get("sprint_status"),
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
        "selected_V31_targets": selected_scope,
        "real_external_constraint_targets_all": real_targets_all,
        "real_external_constraint_targets_in_selected_scope": real_targets_in_scope,
        "selected_V32_target": selected_v32_target,
        "selected_V32_panel": selected_v32_panel,
        "md_ready_targets": [],
        "preflight_failed_checks": failed,
        "provenance_clean": provenance_clean,
        "generated_internal_reports_excluded_from_evidence_claims": True,
        "annotation_only_context_excluded_from_constraint_claims": True,
        "allowed_use_matrix": allowed_use_matrix,
        "target_rows": target_rows,
        "next_operator_readout_decision": decision,
        "locked_interpretation": (
            "V31 prevents self-confirmation leakage by classifying every V30 candidate file by provenance and allowed use. "
            "Generated runtime reports are audit-only; annotation files are role-context-only; only real external constraint/coupling files may support a V32 constraint-backed operator readout. No MD or folding claim is allowed."
        ),
    }


def _write_rows_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "target", "in_selected_V31_scope", "candidate_file_count", "constraint_backed_operator_readout_allowed",
        "preflight_status", "real_external_constraint_or_coupling_files", "annotation_context_files",
        "generated_internal_report_files", "excluded_files", "ready_for_MD", "new_MD_allowed", "claim_allowed",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            out = {key: row.get(key) for key in fields}
            for key in ["real_external_constraint_or_coupling_files", "annotation_context_files", "generated_internal_report_files", "excluded_files"]:
                out[key] = ";".join(str(x) for x in row.get(key, []))
            writer.writerow(out)


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V31 Constraint-Backed Operator Readout Preflight",
        "",
        f"Status: `{cert.get('preflight_status')}`",
        f"Selected V31 targets: `{cert.get('selected_V31_targets')}`",
        f"Selected V32 target: `{cert.get('selected_V32_target')}`",
        f"Real external constraint targets in scope: `{cert.get('real_external_constraint_targets_in_selected_scope')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD allowed: `{cert.get('new_MD_allowed')}`",
        "",
        "## Allowed-use policy",
    ]
    for source_class, allowed in cert.get("allowed_use_matrix", {}).items():
        lines.append(f"- `{source_class}` -> `{allowed}`")
    lines.extend([
        "",
        "## Locked interpretation",
        str(cert.get("locked_interpretation", "")),
        "",
        "## Target rows",
    ])
    for row in cert.get("target_rows", []):
        lines.extend([
            f"### {row.get('target')}",
            f"In selected V31 scope: `{row.get('in_selected_V31_scope')}`",
            f"Preflight status: `{row.get('preflight_status')}`",
            f"Real external constraints: `{row.get('real_external_constraint_or_coupling_files')}`",
            f"Generated internal reports: `{len(row.get('generated_internal_report_files', []))}`",
            f"Allowed-use counts: `{row.get('allowed_use_counts')}`",
            "",
        ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v31_constraint_backed_operator_readout_preflight_certificate.json",
        "classification": out_dir / "v31_candidate_file_provenance_classification.json",
        "allowed_use_matrix": out_dir / "v31_allowed_use_matrix.json",
        "decision": out_dir / "v31_next_operator_readout_decision.json",
        "rows": out_dir / "v31_target_preflight_rows.json",
        "rows_csv": out_dir / "v31_target_preflight_rows.csv",
        "report": out_dir / "V31_CONSTRAINT_BACKED_OPERATOR_READOUT_PREFLIGHT_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["classification"].write_text(json.dumps({"kind": "V31_CANDIDATE_FILE_PROVENANCE_CLASSIFICATION_v0", "targets": cert["target_rows"], "claim_allowed": False}, indent=2, sort_keys=True), encoding="utf-8")
    paths["allowed_use_matrix"].write_text(json.dumps({"kind": "V31_ALLOWED_USE_MATRIX_v0", "allowed_use_matrix": cert["allowed_use_matrix"], "claim_allowed": False}, indent=2, sort_keys=True), encoding="utf-8")
    paths["decision"].write_text(json.dumps(cert["next_operator_readout_decision"], indent=2, sort_keys=True), encoding="utf-8")
    paths["rows"].write_text(json.dumps(cert["target_rows"], indent=2, sort_keys=True), encoding="utf-8")
    _write_rows_csv(paths["rows_csv"], cert["target_rows"])
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v30-cert", type=Path, default=DEFAULT_V30_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    v30 = _read_json(args.v30_cert, "V30 external constraint acquisition certificate")
    cert = build_v31(v30)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "preflight_status": cert["preflight_status"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "new_MD_allowed": cert["new_MD_allowed"],
        "new_MD_recommended": cert["new_MD_recommended"],
        "selected_V31_targets": cert["selected_V31_targets"],
        "real_external_constraint_targets_in_selected_scope": cert["real_external_constraint_targets_in_selected_scope"],
        "selected_V32_target": cert["selected_V32_target"],
        "selected_V32_panel": cert["selected_V32_panel"],
        "provenance_clean": cert["provenance_clean"],
        "certificate": str(paths["certificate"]),
        "classification": str(paths["classification"]),
        "decision": str(paths["decision"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
