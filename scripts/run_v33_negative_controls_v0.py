#!/usr/bin/env python3
from __future__ import annotations

"""V33 negative controls for claim-disabled KcsA operator readout.

This script is intentionally narrow. It does not run MD, does not tune
thresholds, and does not make folding claims. It tests whether the V32/V33
source-import/readout gates refuse obvious leakage, incomplete evidence, and
wrong-target/operator-bucket misuse before any later claim or simulation work.
"""

import argparse
import copy
import importlib.util
import json
import os
from pathlib import Path
from typing import Any, Callable

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V31_CERT = RUN_ROOT / "V31_CONSTRAINT_BACKED_OPERATOR_READOUT_PREFLIGHT" / "v31_constraint_backed_operator_readout_preflight_certificate.json"
DEFAULT_V32_CERT = RUN_ROOT / "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT" / "v32_external_constraint_source_import_preflight_certificate.json"
DEFAULT_V33_CERT = RUN_ROOT / "V33_CONSTRAINT_BACKED_OPERATOR_READOUT" / "v33_constraint_backed_operator_readout_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V33_NEGATIVE_CONTROLS"

PASSED = "V33_NEGATIVE_CONTROLS_PASSED_CLAIM_DISABLED"
FAILED = "V33_NEGATIVE_CONTROLS_FAILED_CLAIM_DISABLED"


def _load_script(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if not spec or not spec.loader:
        raise SystemExit(f"cannot load script module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def _minimal_v31() -> dict[str, Any]:
    return {
        "preflight_status": "V31_CLEAN_ABSTAIN_NO_REAL_EXTERNAL_CONSTRAINTS_FOR_SELECTED_TARGETS",
        "selected_V31_targets": ["XCL1_lymphotactin", "KcsA"],
        "claim_allowed": False,
        "new_md_executed": False,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "native_metrics_used_for_selection": False,
    }


def _target_row(v32_cert: dict[str, Any], target: str = "KcsA") -> dict[str, Any]:
    for row in _as_list(v32_cert.get("target_rows")):
        if isinstance(row, dict) and row.get("target") == target:
            return row
    raise SystemExit(f"V32 certificate does not contain target row: {target}")


def _valid_kcsa_rows(v32_cert: dict[str, Any]) -> list[dict[str, Any]]:
    rows = _as_list(_target_row(v32_cert, "KcsA").get("valid_real_external_constraint_rows"))
    rows = [copy.deepcopy(row) for row in rows if isinstance(row, dict)]
    if not rows:
        raise SystemExit("V32 certificate has no valid KcsA rows; run KcsA source import and V32 first.")
    return rows


def _mutate_kcsa_rows(v32_cert: dict[str, Any], keep: Callable[[dict[str, Any]], bool]) -> dict[str, Any]:
    v32 = copy.deepcopy(v32_cert)
    for target_row in _as_list(v32.get("target_rows")):
        if isinstance(target_row, dict) and target_row.get("target") == "KcsA":
            target_row["valid_real_external_constraint_rows"] = [row for row in _as_list(target_row.get("valid_real_external_constraint_rows")) if isinstance(row, dict) and keep(row)]
            break
    return v32


def _manual_v32_for_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "preflight_status": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_READY_FOR_V33_CLAIM_DISABLED",
        "selected_V31_targets": ["XCL1_lymphotactin", "KcsA"],
        "selected_V33_target": "KcsA",
        "selected_V33_panel": "V33_CONSTRAINT_BACKED_OPERATOR_READOUT",
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "provenance_clean": True,
        "target_rows": [
            {"target": "XCL1_lymphotactin", "valid_real_external_constraint_rows": []},
            {"target": "KcsA", "valid_real_external_constraint_rows": copy.deepcopy(rows)},
        ],
    }


def _make_bad_internal_file() -> str:
    bad_path = RUN_ROOT / "V33_NEGATIVE_CONTROLS" / "internal_poison_source.csv"
    bad_path.parent.mkdir(parents=True, exist_ok=True)
    bad_path.write_text("constraint_id,chain,residue_number,constraint_class\nBAD,A,75,pore_filter_coupling\n", encoding="utf-8")
    return _rel(bad_path)


def _make_annotation_file() -> str:
    annot_path = REPO_ROOT / "data" / "external_annotations" / "KcsA" / "annotation_only_context.csv"
    annot_path.parent.mkdir(parents=True, exist_ok=True)
    annot_path.write_text("target,annotation\nKcsA,membrane potassium channel annotation only\n", encoding="utf-8")
    return _rel(annot_path)


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V33 Negative Controls",
        "",
        f"Status: `{cert.get('negative_control_status')}`",
        f"Passed controls: `{cert.get('passed_control_count')}` / `{cert.get('control_count')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD allowed: `{cert.get('new_MD_allowed')}`",
        "",
        "## Controls",
    ]
    for row in _as_list(cert.get("controls")):
        lines.extend([
            f"### {row.get('control_id')}",
            f"Observed: `{row.get('observed_status')}`",
            f"Expected: `{row.get('expected_status')}`",
            f"Passed: `{row.get('passed')}`",
            f"Reason: {row.get('reason')}",
            "",
        ])
    lines.extend([
        "## Locked interpretation",
        str(cert.get("locked_interpretation", "")),
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_negative_controls(v31_cert: dict[str, Any], v32_cert: dict[str, Any], v33_cert: dict[str, Any]) -> dict[str, Any]:
    v32_mod = _load_script("v32_import_preflight_for_negative_controls", REPO_ROOT / "scripts" / "run_v32_external_constraint_source_import_preflight_v0.py")
    v33_mod = _load_script("v33_readout_for_negative_controls", REPO_ROOT / "scripts" / "run_v33_constraint_backed_operator_readout_v0.py")
    v32_mod.REPO_ROOT = REPO_ROOT
    v33_mod.REPO_ROOT = REPO_ROOT

    controls: list[dict[str, Any]] = []

    # Positive baseline sanity inside the negative-control battery.
    baseline = v33_mod.build_v33(copy.deepcopy(v32_cert))
    controls.append({
        "control_id": "baseline_real_kcsa_sources_pass_before_negative_controls",
        "observed_status": baseline.get("readout_status"),
        "expected_status": "V33_CONSTRAINT_BACKED_OPERATOR_READOUT_PASSED_CLAIM_DISABLED",
        "passed": baseline.get("readout_status") == "V33_CONSTRAINT_BACKED_OPERATOR_READOUT_PASSED_CLAIM_DISABLED"
            and baseline.get("operator_readout_found") is True
            and baseline.get("claim_allowed") is False
            and baseline.get("new_MD_allowed") is False,
        "reason": "Real KcsA pore/filter + assembly/interface imported source rows should pass only as claim-disabled operator readout.",
    })

    missing_pore_v32 = _mutate_kcsa_rows(v32_cert, lambda row: row.get("evidence_type") != "pore_filter_coupling")
    missing_pore = v33_mod.build_v33(missing_pore_v32)
    controls.append({
        "control_id": "missing_pore_filter_bucket_abstains",
        "observed_status": missing_pore.get("readout_status"),
        "expected_status": "V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED",
        "observed_missing_operator_evidence": missing_pore.get("missing_operator_evidence"),
        "passed": missing_pore.get("readout_status") == "V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED"
            and "pore_filter_external_constraint_rows" in _as_list(missing_pore.get("missing_operator_evidence"))
            and missing_pore.get("claim_allowed") is False,
        "reason": "KcsA cannot produce a pore/filter operator readout without the pore/filter source bucket.",
    })

    missing_iface_v32 = _mutate_kcsa_rows(v32_cert, lambda row: row.get("evidence_type") != "assembly_interface_constraint")
    missing_iface = v33_mod.build_v33(missing_iface_v32)
    controls.append({
        "control_id": "missing_assembly_interface_bucket_abstains",
        "observed_status": missing_iface.get("readout_status"),
        "expected_status": "V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED",
        "observed_missing_operator_evidence": missing_iface.get("missing_operator_evidence"),
        "passed": missing_iface.get("readout_status") == "V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED"
            and "assembly_interface_external_constraint_rows" in _as_list(missing_iface.get("missing_operator_evidence"))
            and missing_iface.get("claim_allowed") is False,
        "reason": "KcsA cannot produce a tetramer/interface operator readout without the assembly/interface source bucket.",
    })

    empty_buckets_v32 = _mutate_kcsa_rows(v32_cert, lambda row: False)
    empty_buckets = v33_mod.build_v33(empty_buckets_v32)
    controls.append({
        "control_id": "empty_operator_buckets_abstain",
        "observed_status": empty_buckets.get("readout_status"),
        "expected_status": "V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED",
        "observed_missing_operator_evidence": empty_buckets.get("missing_operator_evidence"),
        "passed": empty_buckets.get("readout_status") == "V33_CLEAN_ABSTAIN_NO_V32_SELECTED_TARGET_OR_CONSTRAINT_BUCKETS_CLAIM_DISABLED"
            and len(_as_list(empty_buckets.get("missing_operator_evidence"))) >= 2
            and empty_buckets.get("claim_allowed") is False,
        "reason": "No source buckets means no operator readout.",
    })

    wrong_target = copy.deepcopy(v32_cert)
    wrong_target["selected_V33_target"] = "XCL1_lymphotactin"
    wrong_target_readout = v33_mod.build_v33(wrong_target)
    controls.append({
        "control_id": "wrong_target_xcl1_blocked_by_v33_v0_scope",
        "observed_status": wrong_target_readout.get("readout_status"),
        "expected_status": "V33_BLOCKED_PRECONDITION_OR_PROVENANCE_FAILURE_CLAIM_DISABLED",
        "observed_failed_checks": wrong_target_readout.get("failed_checks"),
        "passed": wrong_target_readout.get("readout_status") == "V33_BLOCKED_PRECONDITION_OR_PROVENANCE_FAILURE_CLAIM_DISABLED"
            and "selected_target_not_supported_by_V33_v0" in _as_list(wrong_target_readout.get("failed_checks"))
            and wrong_target_readout.get("claim_allowed") is False,
        "reason": "V33 v0 is deliberately scoped to KcsA after V32 selection; XCL1 needs its own state-specific import/readout.",
    })

    # V32 source-import controls.
    real_rows = _valid_kcsa_rows(v32_cert)
    internal_path = _make_bad_internal_file()
    internal_manifest = {
        "kind": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_MANIFEST_v0",
        "claim_allowed": False,
        "manifest_present": True,
        "rows": [
            {**real_rows[0], "file_path": internal_path, "source_name": "Internal runtime poison", "source_url_or_citation": "not_external", "source_date_or_version": "negative_control"},
            real_rows[1],
        ],
    }
    internal_v32 = v32_mod.build_v32(v31_cert or _minimal_v31(), internal_manifest)
    controls.append({
        "control_id": "internal_runtime_source_poison_blocked_by_v32",
        "observed_status": internal_v32.get("preflight_status"),
        "expected_status": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_BLOCKED",
        "observed_failed_checks": internal_v32.get("preflight_failed_checks"),
        "passed": internal_v32.get("preflight_status") == "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_PREFLIGHT_BLOCKED"
            and internal_v32.get("provenance_clean") is False
            and "internal_generated_artifact_supplied_as_external_source" in _as_list(internal_v32.get("preflight_failed_checks"))
            and internal_v32.get("claim_allowed") is False,
        "reason": "Files under first_contact_clean_pharmacotopology_layer_run must never become external evidence.",
    })

    annotation_path = _make_annotation_file()
    annotation_manifest = {
        "kind": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_MANIFEST_v0",
        "claim_allowed": False,
        "manifest_present": True,
        "rows": [
            {
                "target": "KcsA",
                "file_path": annotation_path,
                "evidence_type": "annotation_context",
                "state_or_context": "membrane_pore_annotation_context_only",
                "source_name": "External annotation negative control",
                "source_url_or_citation": "local_external_annotation_context_only_negative_control",
                "source_date_or_version": "negative_control",
                "provenance_notes": "Annotation-only context; must not count as real external constraint.",
            }
        ],
    }
    annotation_v32 = v32_mod.build_v32(v31_cert or _minimal_v31(), annotation_manifest)
    controls.append({
        "control_id": "annotation_only_context_does_not_select_v33",
        "observed_status": annotation_v32.get("preflight_status"),
        "expected_status": "V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED",
        "observed_selected_V33_target": annotation_v32.get("selected_V33_target"),
        "passed": annotation_v32.get("preflight_status") == "V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED"
            and annotation_v32.get("selected_V33_target") is None
            and annotation_v32.get("claim_allowed") is False,
        "reason": "Annotation-only rows can give role context, but must not open a constraint-backed operator readout.",
    })

    only_pore_manifest = {
        "kind": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_MANIFEST_v0",
        "claim_allowed": False,
        "manifest_present": True,
        "rows": [row for row in real_rows if row.get("evidence_type") == "pore_filter_coupling"],
    }
    only_pore_v32 = v32_mod.build_v32(v31_cert or _minimal_v31(), only_pore_manifest)
    controls.append({
        "control_id": "v32_only_pore_bucket_not_ready_for_kcsa_v33",
        "observed_status": only_pore_v32.get("preflight_status"),
        "expected_status": "V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED",
        "observed_selected_V33_target": only_pore_v32.get("selected_V33_target"),
        "passed": only_pore_v32.get("preflight_status") == "V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED"
            and only_pore_v32.get("selected_V33_target") is None,
        "reason": "KcsA needs both pore/filter and assembly/interface source buckets before V33 selection.",
    })

    only_iface_manifest = {
        "kind": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_MANIFEST_v0",
        "claim_allowed": False,
        "manifest_present": True,
        "rows": [row for row in real_rows if row.get("evidence_type") == "assembly_interface_constraint"],
    }
    only_iface_v32 = v32_mod.build_v32(v31_cert or _minimal_v31(), only_iface_manifest)
    controls.append({
        "control_id": "v32_only_interface_bucket_not_ready_for_kcsa_v33",
        "observed_status": only_iface_v32.get("preflight_status"),
        "expected_status": "V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED",
        "observed_selected_V33_target": only_iface_v32.get("selected_V33_target"),
        "passed": only_iface_v32.get("preflight_status") == "V32_CLEAN_ABSTAIN_REAL_EXTERNAL_CONSTRAINT_IMPORT_REQUIRED"
            and only_iface_v32.get("selected_V33_target") is None,
        "reason": "Assembly/interface context alone must not become a KcsA pore/filter readout.",
    })

    claim_locked = all(row.get("passed") is True for row in controls)
    all_claim_md_closed = all(
        obj.get("claim_allowed") is False and obj.get("new_MD_allowed") is False
        for obj in [baseline, missing_pore, missing_iface, empty_buckets, wrong_target_readout, internal_v32, annotation_v32, only_pore_v32, only_iface_v32]
    )

    status = PASSED if claim_locked and all_claim_md_closed else FAILED
    return {
        "kind": "V33_NEGATIVE_CONTROLS_v0",
        "run_mode": "negative_controls_no_MD_no_threshold_tuning_claim_disabled",
        "negative_control_status": status,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row.get("passed") is True),
        "controls": controls,
        "claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
        "baseline_v33_status": v33_cert.get("readout_status"),
        "baseline_v33_operator_readout_found": v33_cert.get("operator_readout_found"),
        "next_action": "run_discriminative_or_randomized_constraint_controls_before_any_MD" if status == PASSED else "fix_failed_negative_controls_before_V34_or_any_claim",
        "locked_interpretation": (
            "These controls validate V32/V33 evidence gating only. Passing them means the system blocks internal/runtime evidence, annotation-only rows, missing buckets, and wrong-target misuse while keeping claim and MD gates closed. It does not prove de novo folding or universal protein-folding prediction."
        ),
    }


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v33_negative_controls_certificate.json"
    report_path = out_dir / "V33_NEGATIVE_CONTROLS_REPORT.md"
    decision_path = out_dir / "v33_negative_controls_next_decision.json"
    decision = {
        "kind": "V33_NEGATIVE_CONTROLS_NEXT_DECISION_v0",
        "decision_status": cert.get("negative_control_status"),
        "next_action": cert.get("next_action"),
        "claim_allowed": False,
        "new_MD_allowed": False,
        "folding_problem_solved": False,
    }
    cert = {
        **cert,
        "next_decision": decision,
        "artifacts": {
            "certificate": str(cert_path),
            "report": str(report_path),
            "decision": str(decision_path),
        },
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, decision)
    _write_report(report_path, cert)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V33 claim-disabled negative controls.")
    parser.add_argument("--v31-certificate", type=Path, default=DEFAULT_V31_CERT)
    parser.add_argument("--v32-certificate", type=Path, default=DEFAULT_V32_CERT)
    parser.add_argument("--v33-certificate", type=Path, default=DEFAULT_V33_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    v31 = _read_json(args.v31_certificate, "V31 certificate") if args.v31_certificate.exists() else _minimal_v31()
    v32 = _read_json(args.v32_certificate, "V32 certificate")
    v33 = _read_json(args.v33_certificate, "V33 certificate")
    cert = build_negative_controls(v31, v32, v33)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert.get("kind"),
        "negative_control_status": cert.get("negative_control_status"),
        "control_count": cert.get("control_count"),
        "passed_control_count": cert.get("passed_control_count"),
        "claim_allowed": cert.get("claim_allowed"),
        "new_MD_allowed": cert.get("new_MD_allowed"),
        "positive_folding_evidence_found": cert.get("positive_folding_evidence_found"),
        "folding_problem_solved": cert.get("folding_problem_solved"),
        "next_action": cert.get("next_action"),
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
