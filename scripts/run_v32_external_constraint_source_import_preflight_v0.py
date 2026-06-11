#!/usr/bin/env python3
from __future__ import annotations

"""V32 external constraint source import preflight.

V31 correctly clean-abstained for XCL1/KcsA because the selected-scope
candidate files were generated internal reports, not real external constraints.
V32 is therefore not a constraint-backed readout. It is an import/acquisition
preflight for real external constraint sources.

It never fabricates couplings, never promotes annotations to constraints, never
runs MD, and keeps claim_allowed=false.
"""

import argparse
import csv
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V31_CERT = RUN_ROOT / "V31_CONSTRAINT_BACKED_OPERATOR_READOUT_PREFLIGHT" / "v31_constraint_backed_operator_readout_preflight_certificate.json"
DEFAULT_IMPORT_MANIFEST = REPO_ROOT / "data" / "external_constraints" / "v32_external_constraint_source_import_manifest.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT"

SELECTED_SCOPE_DEFAULT = ["XCL1_lymphotactin", "KcsA"]
TARGET_PRIORITY = {"XCL1_lymphotactin": 1, "KcsA": 2, "p53_TAD_MDM2": 3, "4AKE": 4, "1UBQ": 5, "1CLL": 6}

REAL_CONSTRAINT_TYPES = {
    "external_coupling",
    "contact_constraint",
    "state_specific_constraint",
    "assembly_interface_constraint",
    "pore_filter_coupling",
}
ROLE_CONTEXT_TYPES = {"annotation_context", "state_label", "condition_label", "structure_context", "assembly_context"}
DERIVATION_ONLY_TYPES = {"alignment_source", "msa_source"}

ALLOWED_ROOTS = (
    "data/external_constraints/",
    "data/external_annotations/",
    "data/independent_contact_sources/",
    "data/v16_pressure_targets/",
    "external_msa/",
)

TEMPLATE_ROWS = [
    {
        "target": "XCL1_lymphotactin",
        "file_path": "data/external_constraints/XCL1_lymphotactin/state_A/<real_external_constraint_file>.csv",
        "evidence_type": "state_specific_constraint",
        "state_or_context": "state_A_chemokine_like",
        "source_name": "REPLACE_WITH_DATABASE_OR_PAPER_OR_PIPELINE_NAME",
        "source_url_or_citation": "REPLACE_WITH_STABLE_SOURCE_OR_CITATION",
        "source_date_or_version": "REPLACE_WITH_VERSION_OR_ACCESS_DATE",
        "provenance_notes": "Required: real external source, not first_contact runtime report.",
    },
    {
        "target": "XCL1_lymphotactin",
        "file_path": "data/external_constraints/XCL1_lymphotactin/state_B/<real_external_constraint_file>.csv",
        "evidence_type": "state_specific_constraint",
        "state_or_context": "state_B_dimer_beta_sandwich_like",
        "source_name": "REPLACE_WITH_DATABASE_OR_PAPER_OR_PIPELINE_NAME",
        "source_url_or_citation": "REPLACE_WITH_STABLE_SOURCE_OR_CITATION",
        "source_date_or_version": "REPLACE_WITH_VERSION_OR_ACCESS_DATE",
        "provenance_notes": "Both XCL1 states are required before state-switch constraint-backed readout.",
    },
    {
        "target": "KcsA",
        "file_path": "data/external_constraints/KcsA/pore_filter/<real_external_coupling_file>.csv",
        "evidence_type": "pore_filter_coupling",
        "state_or_context": "TVGYG_selectivity_filter_or_pore_filter",
        "source_name": "REPLACE_WITH_DATABASE_OR_PAPER_OR_PIPELINE_NAME",
        "source_url_or_citation": "REPLACE_WITH_STABLE_SOURCE_OR_CITATION",
        "source_date_or_version": "REPLACE_WITH_VERSION_OR_ACCESS_DATE",
        "provenance_notes": "Must support pore/filter/operator context, not only internal report text.",
    },
    {
        "target": "KcsA",
        "file_path": "data/external_constraints/KcsA/assembly_interface/<real_external_constraint_file>.csv",
        "evidence_type": "assembly_interface_constraint",
        "state_or_context": "tetramer_or_chain_interface_context",
        "source_name": "REPLACE_WITH_DATABASE_OR_PAPER_OR_PIPELINE_NAME",
        "source_url_or_citation": "REPLACE_WITH_STABLE_SOURCE_OR_CITATION",
        "source_date_or_version": "REPLACE_WITH_VERSION_OR_ACCESS_DATE",
        "provenance_notes": "Required before tetramer/assembly-backed claims. Still no whole-channel fold claim.",
    },
]


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


def _target_priority(target: str) -> int:
    return TARGET_PRIORITY.get(target, 99)


def _target_import_status(target: str, valid_real: list[dict[str, Any]], role_context: list[dict[str, Any]], derivation: list[dict[str, Any]]) -> tuple[str, bool, list[str]]:
    missing: list[str] = []
    if target == "XCL1_lymphotactin":
        contexts = " ".join(str(row.get("state_or_context", "")).lower() for row in valid_real)
        has_a = any(token in contexts for token in ["state_a", "chemokine", "monomer"])
        has_b = any(token in contexts for token in ["state_b", "dimer", "beta_sandwich", "beta-sandwich"])
        if not has_a:
            missing.append("state_A_real_external_constraint")
        if not has_b:
            missing.append("state_B_real_external_constraint")
        ready = not missing
        return ("target_ready_for_V33_constraint_backed_readout" if ready else "partial_or_missing_state_specific_external_constraints"), ready, missing
    if target == "KcsA":
        contexts = " ".join((str(row.get("state_or_context", "")) + " " + str(row.get("evidence_type", ""))).lower() for row in valid_real)
        has_filter_or_coupling = any(token in contexts for token in ["pore", "filter", "tvgyg", "selectivity", "coupling", "contact"])
        has_interface_or_assembly = any(token in contexts for token in ["interface", "assembly", "tetramer", "chain"])
        if not has_filter_or_coupling:
            missing.append("pore_filter_or_external_coupling_constraint")
        if not has_interface_or_assembly:
            missing.append("assembly_or_chain_interface_constraint")
        ready = not missing
        return ("target_ready_for_V33_constraint_backed_readout" if ready else "partial_or_missing_kcsa_coupling_interface_constraints"), ready, missing
    ready = bool(valid_real)
    return ("target_ready_for_V33_constraint_backed_readout" if ready else "no_real_external_constraints_imported"), ready, [] if ready else ["real_external_constraint_or_coupling"]


def _load_import_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {
            "kind": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_MANIFEST_v0",
            "claim_allowed": False,
            "manifest_present": False,
            "rows": [],
            "note": "No import manifest exists yet. Use the generated template under data/external_constraints/.",
        }
    data = _read_json(path, "V32 external constraint source import manifest")
    rows = _as_list(data.get("rows"))
    return {
        **data,
        "manifest_present": True,
        "rows": rows,
    }


def _row_path_exists(rel_path: str) -> bool:
    return (REPO_ROOT / rel_path).exists()


def classify_import_row(row: dict[str, Any], selected_scope: list[str]) -> dict[str, Any]:
    target = str(row.get("target", "")).strip()
    file_path = _norm_rel(str(row.get("file_path", "")).strip())
    evidence_type = str(row.get("evidence_type", "")).strip()
    source_name = str(row.get("source_name", "")).strip()
    source_ref = str(row.get("source_url_or_citation", "")).strip()
    state = str(row.get("state_or_context", "")).strip()
    errors: list[str] = []

    if not target:
        errors.append("target_missing")
    if target and target not in selected_scope:
        errors.append("target_not_in_selected_V31_scope")
    if not file_path:
        errors.append("file_path_missing")
    if file_path.startswith("first_contact_clean_pharmacotopology_layer_run/") or file_path.startswith("data/locked_runtime_certificates/"):
        errors.append("internal_runtime_or_exported_certificate_not_external_source")
    if file_path and not file_path.startswith(ALLOWED_ROOTS):
        errors.append("file_path_not_under_allowed_external_source_root")
    if file_path and not _row_path_exists(file_path):
        errors.append("file_path_does_not_exist")
    if not evidence_type:
        errors.append("evidence_type_missing")
    if not source_name or source_name.startswith("REPLACE_WITH"):
        errors.append("source_name_missing_or_template_placeholder")
    if not source_ref or source_ref.startswith("REPLACE_WITH"):
        errors.append("source_url_or_citation_missing_or_template_placeholder")

    if evidence_type in REAL_CONSTRAINT_TYPES and not errors:
        evidence_source_class = "real_external_constraint_or_coupling"
        allowed_use = "allowed_for_constraint_backed_operator_readout"
        counts_as_real_external_constraint = True
    elif evidence_type in DERIVATION_ONLY_TYPES and not errors:
        evidence_source_class = "real_external_alignment_source"
        allowed_use = "allowed_for_constraint_derivation_preflight_only"
        counts_as_real_external_constraint = False
    elif evidence_type in ROLE_CONTEXT_TYPES and not errors:
        evidence_source_class = "annotation_only_external_context"
        allowed_use = "allowed_for_role_context_only"
        counts_as_real_external_constraint = False
    else:
        evidence_source_class = "invalid_or_unusable_import_row"
        allowed_use = "excluded"
        counts_as_real_external_constraint = False
        if evidence_type and evidence_type not in REAL_CONSTRAINT_TYPES | DERIVATION_ONLY_TYPES | ROLE_CONTEXT_TYPES:
            errors.append("unsupported_evidence_type")

    return {
        "target": target,
        "file_path": file_path,
        "evidence_type": evidence_type,
        "state_or_context": state,
        "source_name": source_name,
        "source_url_or_citation": source_ref,
        "source_date_or_version": row.get("source_date_or_version"),
        "provenance_notes": row.get("provenance_notes"),
        "in_selected_V31_scope": target in selected_scope,
        "file_exists": _row_path_exists(file_path) if file_path else False,
        "evidence_source_class": evidence_source_class,
        "allowed_use": allowed_use,
        "counts_as_real_external_constraint": counts_as_real_external_constraint,
        "validation_errors": sorted(set(errors)),
        "valid_import_row": not errors,
        "claim_allowed": False,
    }


def _auto_select_v33(rows: list[dict[str, Any]]) -> str | None:
    eligible = [row for row in rows if row.get("target_ready_for_V33") is True]
    if not eligible:
        return None
    eligible = sorted(eligible, key=lambda row: (-len(_as_list(row.get("valid_real_external_constraint_rows"))), _target_priority(str(row.get("target")))))
    return str(eligible[0].get("target"))


def build_v32(v31: dict[str, Any], import_manifest: dict[str, Any]) -> dict[str, Any]:
    failed: list[str] = []
    if v31.get("preflight_status") not in {
        "V31_CLEAN_ABSTAIN_NO_REAL_EXTERNAL_CONSTRAINTS_FOR_SELECTED_TARGETS",
        "V31_CONSTRAINT_BACKED_OPERATOR_READOUT_PREFLIGHT_PASSED_CLAIM_DISABLED",
    }:
        failed.append("V31_preflight_status_acceptable")
    for key in ["claim_allowed", "new_md_executed", "new_MD_allowed", "new_MD_recommended", "native_metrics_used_for_selection"]:
        value = v31.get(key)
        if value not in (False, None):
            failed.append(f"V31_{key}_must_be_false")

    selected_scope = [str(x) for x in _as_list(v31.get("selected_V31_targets"))] or list(SELECTED_SCOPE_DEFAULT)
    rows_in = [row for row in _as_list(import_manifest.get("rows")) if isinstance(row, dict)]
    classified_import_rows = [classify_import_row(row, selected_scope) for row in rows_in]

    target_rows: list[dict[str, Any]] = []
    for target in sorted(selected_scope, key=_target_priority):
        target_classified = [row for row in classified_import_rows if row.get("target") == target]
        valid_real = [row for row in target_classified if row.get("counts_as_real_external_constraint") is True and row.get("valid_import_row") is True]
        role_context = [row for row in target_classified if row.get("allowed_use") == "allowed_for_role_context_only" and row.get("valid_import_row") is True]
        derivation = [row for row in target_classified if row.get("allowed_use") == "allowed_for_constraint_derivation_preflight_only" and row.get("valid_import_row") is True]
        invalid = [row for row in target_classified if row.get("valid_import_row") is False]
        status, ready, missing = _target_import_status(target, valid_real, role_context, derivation)
        if not target_classified:
            status = "no_import_rows_for_target"
            ready = False
            missing = ["real_external_constraint_import_rows"]
        target_rows.append({
            "target": target,
            "import_row_count": len(target_classified),
            "valid_real_external_constraint_rows": valid_real,
            "role_context_rows": role_context,
            "constraint_derivation_only_rows": derivation,
            "invalid_or_excluded_rows": invalid,
            "target_import_status": status,
            "target_ready_for_V33": ready,
            "missing_for_V33": missing,
            "ready_for_MD": False,
            "new_MD_allowed": False,
            "new_MD_recommended": False,
            "claim_allowed": False,
        })

    invalid_misuse = [row for row in classified_import_rows if "internal_runtime_or_exported_certificate_not_external_source" in _as_list(row.get("validation_errors"))]
    unsupported_or_missing = [row for row in classified_import_rows if row.get("valid_import_row") is False]
    if invalid_misuse:
        failed.append("internal_generated_artifact_supplied_as_external_source")

    selected_v33_target = _auto_select_v33(target_rows)
    if failed:
        status = "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_BLOCKED"
    elif selected_v33_target:
        status = "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED"
    else:
        status = "V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED"

    selected_v33_panel = "V33_CONSTRAINT_BACKED_OPERATOR_READOUT" if selected_v33_target else None
    next_action = (
        "run_V33_constraint_backed_operator_readout_for_auto_selected_target_no_MD"
        if selected_v33_target else
        "place_real_external_constraint_files_under_data/external_constraints_and_fill_v32_import_manifest_then_rerun_V32"
    )

    import_requirements = {
        "XCL1_lymphotactin": [
            "state_A real external constraint/coupling with state_or_context marking chemokine/monomer context",
            "state_B real external constraint/coupling with state_or_context marking dimer/beta-sandwich context",
            "no mixed-state pooling",
        ],
        "KcsA": [
            "pore/filter/selectivity or external coupling/contact support",
            "assembly/interface/tetramer/chain-interface context constraint",
            "no whole-channel/tetramer fold claim from annotation alone",
        ],
    }

    decision = {
        "kind": "V32_NEXT_CONSTRAINT_BACKED_OPERATOR_READOUT_DECISION_v0",
        "decision_status": status,
        "selected_V33_target": selected_v33_target,
        "selected_V33_panel": selected_v33_panel,
        "next_action": next_action,
        "reason": (
            "V32 imports only user-provided or locally present real external constraint sources. "
            "It does not treat V30/V31 runtime reports as evidence and does not run MD."
        ),
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "claim_allowed": False,
    }

    return {
        "kind": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_v0",
        "run_mode": "external_constraint_source_import_preflight_no_simulation_no_threshold_tuning",
        "source_v31_status": v31.get("preflight_status"),
        "import_manifest_present": bool(import_manifest.get("manifest_present")),
        "import_manifest_path": "data/external_constraints/v32_external_constraint_source_import_manifest.json",
        "selected_V31_targets": selected_scope,
        "preflight_status": status,
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
        "provenance_clean": not invalid_misuse,
        "preflight_failed_checks": failed,
        "classified_import_rows": classified_import_rows,
        "invalid_or_excluded_import_rows": unsupported_or_missing,
        "target_rows": target_rows,
        "selected_V33_target": selected_v33_target,
        "selected_V33_panel": selected_v33_panel,
        "md_ready_targets": [],
        "import_requirements": import_requirements,
        "next_constraint_backed_operator_readout_decision": decision,
        "locked_interpretation": (
            "V32 is an acquisition/import preflight, not a folding win. It requires real external constraints for XCL1 or KcsA before any constraint-backed operator readout. Generated internal reports remain audit-only; annotations remain context-only; no MD and no claim are allowed."
        ),
    }


def _write_template(path: Path) -> None:
    if path.exists():
        return
    template = {
        "kind": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_MANIFEST_v0",
        "claim_allowed": False,
        "instructions": [
            "Replace template rows with real external source rows.",
            "Do not use files under first_contact_clean_pharmacotopology_layer_run as evidence.",
            "Do not use data/locked_runtime_certificates as evidence.",
            "Real external constraint files should live under data/external_constraints/<target>/...",
            "Annotation-only files are allowed for role context only, not constraint-backed evidence claims.",
        ],
        "rows": TEMPLATE_ROWS,
    }
    _write_json(path, template)


def _write_rows_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "target", "import_row_count", "target_import_status", "target_ready_for_V33", "missing_for_V33",
        "valid_real_external_constraint_count", "role_context_count", "derivation_only_count", "invalid_or_excluded_count",
        "ready_for_MD", "new_MD_allowed", "claim_allowed",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "target": row.get("target"),
                "import_row_count": row.get("import_row_count"),
                "target_import_status": row.get("target_import_status"),
                "target_ready_for_V33": row.get("target_ready_for_V33"),
                "missing_for_V33": ";".join(str(x) for x in row.get("missing_for_V33", [])),
                "valid_real_external_constraint_count": len(row.get("valid_real_external_constraint_rows", [])),
                "role_context_count": len(row.get("role_context_rows", [])),
                "derivation_only_count": len(row.get("constraint_derivation_only_rows", [])),
                "invalid_or_excluded_count": len(row.get("invalid_or_excluded_rows", [])),
                "ready_for_MD": row.get("ready_for_MD"),
                "new_MD_allowed": row.get("new_MD_allowed"),
                "claim_allowed": row.get("claim_allowed"),
            })


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V32 External Constraint Source Import Preflight",
        "",
        f"Status: `{cert.get('preflight_status')}`",
        f"Selected V31 targets: `{cert.get('selected_V31_targets')}`",
        f"Selected V33 target: `{cert.get('selected_V33_target')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD allowed: `{cert.get('new_MD_allowed')}`",
        "",
        "## Locked interpretation",
        str(cert.get("locked_interpretation", "")),
        "",
        "## Target rows",
    ]
    for row in cert.get("target_rows", []):
        lines.extend([
            f"### {row.get('target')}",
            f"Status: `{row.get('target_import_status')}`",
            f"Ready for V33: `{row.get('target_ready_for_V33')}`",
            f"Missing: `{row.get('missing_for_V33')}`",
            f"Valid real constraints: `{len(row.get('valid_real_external_constraint_rows', []))}`",
            f"Invalid/excluded rows: `{len(row.get('invalid_or_excluded_rows', []))}`",
            "",
        ])
    lines.extend(["## Next decision", json.dumps(cert.get("next_constraint_backed_operator_readout_decision", {}), indent=2, sort_keys=True)])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any], import_manifest_path: Path) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_template(import_manifest_path)
    paths = {
        "certificate": out_dir / "v32_external_constraint_source_import_preflight_certificate.json",
        "import_template": import_manifest_path,
        "classified_import_rows": out_dir / "v32_classified_import_rows.json",
        "target_rows": out_dir / "v32_target_import_preflight_rows.json",
        "target_rows_csv": out_dir / "v32_target_import_preflight_rows.csv",
        "import_requirements": out_dir / "v32_external_constraint_import_requirements.json",
        "decision": out_dir / "v32_next_constraint_backed_operator_readout_decision.json",
        "report": out_dir / "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    _write_json(paths["classified_import_rows"], {"kind": "V32_CLASSIFIED_IMPORT_ROWS_v0", "rows": cert["classified_import_rows"], "claim_allowed": False})
    _write_json(paths["target_rows"], cert["target_rows"])
    _write_rows_csv(paths["target_rows_csv"], cert["target_rows"])
    _write_json(paths["import_requirements"], {"kind": "V32_EXTERNAL_CONSTRAINT_IMPORT_REQUIREMENTS_v0", "requirements": cert["import_requirements"], "claim_allowed": False})
    _write_json(paths["decision"], cert["next_constraint_backed_operator_readout_decision"])
    _write_json(paths["certificate"], cert)
    _write_report(paths["report"], cert)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v31-cert", type=Path, default=DEFAULT_V31_CERT)
    parser.add_argument("--import-manifest", type=Path, default=DEFAULT_IMPORT_MANIFEST)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    v31 = _read_json(args.v31_cert, "V31 constraint-backed preflight certificate")
    _write_template(args.import_manifest)
    import_manifest = _load_import_manifest(args.import_manifest)
    cert = build_v32(v31, import_manifest)
    paths = write_outputs(args.out_dir, cert, args.import_manifest)
    print(json.dumps({
        "kind": cert["kind"],
        "preflight_status": cert["preflight_status"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "new_MD_allowed": cert["new_MD_allowed"],
        "new_MD_recommended": cert["new_MD_recommended"],
        "selected_V31_targets": cert["selected_V31_targets"],
        "selected_V33_target": cert["selected_V33_target"],
        "selected_V33_panel": cert["selected_V33_panel"],
        "provenance_clean": cert["provenance_clean"],
        "import_manifest_path": cert["import_manifest_path"],
        "certificate": str(paths["certificate"]),
        "decision": str(paths["decision"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
