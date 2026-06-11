#!/usr/bin/env python3
from __future__ import annotations

"""V36 real evidence dossier gate.

The gate validates external source dossiers for three hard classes without
using native coordinates, native metrics, internal runtime reports, or MD.
"""

import argparse
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "external_evidence_dossiers"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V36_REAL_EVIDENCE_DOSSIER_GATE"

READY_STATUS = "V36_REAL_EVIDENCE_DOSSIERS_READY_CLAIM_DISABLED"
PARTIAL_STATUS = "V36_PARTIAL_EVIDENCE_DOSSIERS_CLEAN_ABSTAIN"
BLOCKED_LEAKAGE_STATUS = "V36_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"
BLOCKED_UNTRUSTED_STATUS = "V36_BLOCKED_PLACEHOLDER_OR_UNTRUSTED_SOURCES"

TARGET_ORDER = ["KcsA", "XCL1_lymphotactin", "alpha_synuclein_SNCA"]

TARGET_REQUIREMENTS = {
    "KcsA": [
        "sequence_or_family_identity",
        "ion_selectivity_context",
        "membrane_topology_context",
        "filter_or_signature_context",
    ],
    "XCL1_lymphotactin": [
        "state_A_function_context",
        "state_B_function_context",
        "metamorphic_two_state_context",
        "no_mixed_state_pooling_rule",
    ],
    "alpha_synuclein_SNCA": [
        "intrinsic_disorder_context",
        "ensemble_context",
        "disorder_to_order_context",
        "no_single_native_fold_rule",
    ],
}

ALLOWED_SOURCE_TYPES = {
    "UniProt sequence/features/function annotations",
    "DisProt disorder annotations",
    "Pfam/InterPro family/domain signatures",
    "external MSA/family alignments",
    "external sequence conservation signatures",
    "literature-derived state/function annotations",
    "NMR/SAXS/FRET/crosslinking metadata",
    "experimentally described motifs/state labels",
}

COORDINATE_SOURCE_TYPES = {
    "PDB coordinate-derived contact tables",
    "RCSB coordinate contacts",
    "AlphaFold predicted coordinates",
    "ESMFold predicted coordinates",
    "RoseTTAFold predicted coordinates",
}

PLACEHOLDER_TOKENS = {"placeholder", "todo", "tbd", "citation needed", "example.com", "dummy", "unknown"}
COORDINATE_TOKENS = {
    "1bl8",
    "rcsb",
    "coordinate-derived",
    "coordinate_derived",
    "pdb coordinate",
    "pdb_coordinate",
    ".pdb",
    "alphafold",
    "esmfold",
    "rosettafold",
    "rmsd",
    "tm-score",
    "native metric",
}


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


def _bool_true(value: Any) -> bool:
    if value is True:
        return True
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return False


def _bool_false(value: Any) -> bool:
    if value is False:
        return True
    if isinstance(value, str):
        return value.strip().lower() in {"false", "0", "no"}
    return False


def _contains_placeholder(value: Any) -> bool:
    text = str(value or "").strip().lower()
    if not text:
        return True
    return any(token in text for token in PLACEHOLDER_TOKENS)


def _row_text_for_leakage(row: dict[str, Any]) -> str:
    keys = ["file_path", "source_type", "source_name", "source_url_or_citation", "evidence_statement"]
    return " ".join(str(row.get(key, "")) for key in keys).lower()


def _contains_coordinate_token(row: dict[str, Any]) -> bool:
    text = _row_text_for_leakage(row)
    return any(token in text for token in COORDINATE_TOKENS)


def _contains_internal_runtime(row: dict[str, Any]) -> bool:
    text = _row_text_for_leakage(row)
    return "first_contact_clean_pharmacotopology_layer_run" in text or _bool_true(row.get("internal_runtime_source"))


def _stable_source(value: Any) -> bool:
    text = str(value or "").strip().lower()
    if len(text) < 12:
        return False
    stable_tokens = ["http", "doi", "pubmed", "pmc", "uniprot", "disprot", "pfam", "interpro", "proc natl acad", "j biol chem", "curr opin"]
    return any(token in text for token in stable_tokens)


def _target_grammar_failures(target: str, dossier: dict[str, Any]) -> list[str]:
    failures: list[str] = []
    rules = dossier.get("grammar_rules")
    if not isinstance(rules, dict):
        return ["grammar_rules_missing"]
    if target == "XCL1_lymphotactin" and rules.get("mixed_state_pooling_allowed") is not False:
        failures.append("mixed_state_pooling_forbidden")
    if target == "alpha_synuclein_SNCA" and rules.get("single_native_fold_model_allowed") is not False:
        failures.append("single_native_fold_grammar_forbidden_for_SNCA")
    if target == "alpha_synuclein_SNCA" and rules.get("ensemble_not_single_fold") is not True:
        failures.append("SNCA_ensemble_not_single_fold_rule_missing")
    if rules.get("native_metrics_before_selection_allowed") is True:
        failures.append("native_metrics_before_selection_allowed_by_dossier")
    return failures


def evaluate_target(target: str, package: dict[str, Any]) -> dict[str, Any]:
    manifest = package.get("source_manifest") if isinstance(package.get("source_manifest"), dict) else {}
    dossier = package.get("evidence_dossier") if isinstance(package.get("evidence_dossier"), dict) else {}
    rows = [row for row in _as_list(manifest.get("rows")) if isinstance(row, dict)]
    required = TARGET_REQUIREMENTS[target]

    failed_checks: list[str] = []
    placeholder_count = 0
    coordinate_count = 0
    internal_count = 0
    native_metric_count = 0
    trusted_external_count = 0
    source_types: set[str] = set()
    buckets: set[str] = set()

    required_fields = {
        "target",
        "source_id",
        "evidence_bucket",
        "source_type",
        "source_name",
        "source_url_or_citation",
        "source_date_or_version",
        "evidence_statement",
        "allowed_use",
        "provenance_notes",
        "coordinate_derived",
        "internal_runtime_source",
        "native_metrics_used_for_selection",
        "coordinate_truth_used_before_selection",
        "claim_allowed",
    }

    for index, row in enumerate(rows):
        row_failures: list[str] = []
        missing = sorted(field for field in required_fields if field not in row)
        if missing:
            row_failures.append("missing_required_fields:" + ",".join(missing))

        source_type = str(row.get("source_type", "")).strip()
        source_types.add(source_type)
        if row.get("target") != target:
            row_failures.append("target_mismatch")
        if source_type in COORDINATE_SOURCE_TYPES:
            row_failures.append("coordinate_source_type_supplied")
        elif source_type not in ALLOWED_SOURCE_TYPES:
            row_failures.append("unsupported_source_type")
        if _bool_true(row.get("coordinate_derived")) or _bool_true(row.get("coordinate_truth_used_before_selection")) or _contains_coordinate_token(row):
            row_failures.append("coordinate_or_prediction_leakage")
        if _contains_internal_runtime(row):
            row_failures.append("internal_runtime_source_supplied")
        if _bool_true(row.get("native_metrics_used_for_selection")):
            row_failures.append("native_metrics_used_for_selection")
        if not _bool_false(row.get("claim_allowed")):
            row_failures.append("claim_allowed_field_not_false")

        for key in ["source_name", "source_url_or_citation", "source_date_or_version"]:
            if _contains_placeholder(row.get(key)):
                row_failures.append(f"placeholder_or_missing_{key}")
        if not _stable_source(row.get("source_url_or_citation")):
            row_failures.append("source_url_or_citation_not_stable")

        if any(item in row_failures for item in ["coordinate_source_type_supplied", "coordinate_or_prediction_leakage"]):
            coordinate_count += 1
        if "internal_runtime_source_supplied" in row_failures:
            internal_count += 1
        if "native_metrics_used_for_selection" in row_failures:
            native_metric_count += 1
        if any(item.startswith("placeholder_or_missing") or item in {"unsupported_source_type", "source_url_or_citation_not_stable"} for item in row_failures):
            placeholder_count += 1

        if not row_failures:
            trusted_external_count += 1
            bucket = str(row.get("evidence_bucket", "")).strip()
            if bucket:
                buckets.add(bucket)

        failed_checks.extend(f"{target}:row_{index}:{failure}" for failure in sorted(set(row_failures)))

    grammar_failures = _target_grammar_failures(target, dossier)
    failed_checks.extend(f"{target}:{failure}" for failure in grammar_failures)

    missing_buckets = [bucket for bucket in required if bucket not in buckets]
    if coordinate_count or internal_count or native_metric_count or grammar_failures:
        target_status = "blocked"
    elif placeholder_count:
        target_status = "blocked"
    elif missing_buckets:
        target_status = "partial"
        failed_checks.extend(f"{target}:missing_bucket:{bucket}" for bucket in missing_buckets)
    else:
        target_status = "ready"

    return {
        "target": target,
        "target_status": target_status,
        "source_row_count": len(rows),
        "trusted_external_source_count": trusted_external_count,
        "coordinate_derived_source_count": coordinate_count,
        "internal_runtime_source_count": internal_count,
        "placeholder_source_count": placeholder_count,
        "native_metric_source_count": native_metric_count,
        "source_types_used": sorted(source_types),
        "required_buckets": required,
        "present_buckets": sorted(buckets),
        "missing_buckets": missing_buckets,
        "failed_checks": sorted(set(failed_checks)),
    }


def load_packages(data_root: Path = DATA_ROOT) -> dict[str, dict[str, Any]]:
    packages: dict[str, dict[str, Any]] = {}
    for target in TARGET_ORDER:
        target_dir = data_root / target
        packages[target] = {
            "source_manifest": _read_json(target_dir / "source_manifest.json", f"{target} source manifest"),
            "evidence_dossier": _read_json(target_dir / "evidence_dossier.json", f"{target} evidence dossier"),
        }
    return packages


def build_v36(packages: dict[str, dict[str, Any]]) -> dict[str, Any]:
    target_results = [evaluate_target(target, packages.get(target, {})) for target in TARGET_ORDER]
    targets_ready = [row["target"] for row in target_results if row["target_status"] == "ready"]
    targets_partial = [row["target"] for row in target_results if row["target_status"] == "partial"]
    targets_blocked = [row["target"] for row in target_results if row["target_status"] == "blocked"]

    coordinate_count = sum(int(row["coordinate_derived_source_count"]) for row in target_results)
    internal_count = sum(int(row["internal_runtime_source_count"]) for row in target_results)
    native_metric_count = sum(int(row["native_metric_source_count"]) for row in target_results)
    placeholder_count = sum(int(row["placeholder_source_count"]) for row in target_results)
    external_count = sum(int(row["trusted_external_source_count"]) for row in target_results)
    source_counts_by_target = {row["target"]: row["trusted_external_source_count"] for row in target_results}
    source_types_by_target = {row["target"]: row["source_types_used"] for row in target_results}

    failed_checks: list[str] = []
    for row in target_results:
        failed_checks.extend(row["failed_checks"])

    if coordinate_count or internal_count or native_metric_count:
        status = BLOCKED_LEAKAGE_STATUS
    elif placeholder_count or targets_blocked:
        status = BLOCKED_UNTRUSTED_STATUS
    elif len(targets_ready) == len(TARGET_ORDER):
        status = READY_STATUS
    else:
        status = PARTIAL_STATUS

    controls = run_v36_controls()
    if not all(control.get("passed") is True for control in controls):
        failed_checks.append("V36_required_controls_failed")

    if status == READY_STATUS:
        next_action = "use_V36_dossiers_to_design_target_specific_operator_grammar_no_claim_no_MD"
    elif status == PARTIAL_STATUS:
        next_action = "complete_missing_real_external_dossier_buckets_before_operator_grammar"
    elif status == BLOCKED_LEAKAGE_STATUS:
        next_action = "remove_coordinate_internal_or_native_metric_leakage_before_V36"
    else:
        next_action = "replace_placeholder_or_invalid_target_grammar_sources_before_V36"

    return {
        "kind": "V36_REAL_EVIDENCE_DOSSIER_GATE_v0",
        "run_mode": "real_external_evidence_dossier_gate_no_coordinates_no_MD_claim_disabled",
        "control_status": status,
        "target_count": len(TARGET_ORDER),
        "ready_target_count": len(targets_ready),
        "partial_target_count": len(targets_partial),
        "blocked_target_count": len(targets_blocked),
        "targets_ready": targets_ready,
        "targets_partial": targets_partial,
        "targets_blocked": targets_blocked,
        "external_source_count": external_count,
        "source_counts_by_target": source_counts_by_target,
        "source_types_by_target": source_types_by_target,
        "coordinate_derived_source_count": coordinate_count,
        "internal_runtime_source_count": internal_count,
        "placeholder_source_count": placeholder_count,
        "native_metric_source_count": native_metric_count,
        "target_results": target_results,
        "controls": controls,
        "control_count": len(controls),
        "passed_control_count": sum(1 for control in controls if control.get("passed") is True),
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
        "next_action": next_action,
        "failed_checks": sorted(set(failed_checks)),
        "locked_interpretation": (
            "V36 readiness means the repo now has real external, non-coordinate dossiers for KcsA, XCL1, and SNCA operator grammar. "
            "It is not a folding prediction win, does not run MD, does not use native coordinates for selection, and does not solve protein folding."
        ),
    }


def _minimal_row(target: str, bucket: str, source_id: str = "source") -> dict[str, Any]:
    return {
        "target": target,
        "source_id": source_id,
        "evidence_bucket": bucket,
        "source_type": "literature-derived state/function annotations",
        "source_name": f"{target} trusted external source {bucket}",
        "source_url_or_citation": "PubMed stable literature citation",
        "source_date_or_version": "accessed_2026-06-11",
        "source_boundary": "external_noncoordinate_scientific_annotation_or_state_evidence",
        "allowed_use": "V36_claim_disabled_target_operator_grammar_dossier_only",
        "provenance_notes": "External evidence only; no coordinate contacts, no internal runtime reports, no MD.",
        "evidence_statement": f"{target} evidence for {bucket}",
        "grammar_role": bucket,
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "claim_allowed": False,
    }


def _ready_fixture_packages() -> dict[str, dict[str, Any]]:
    packages: dict[str, dict[str, Any]] = {}
    for target in TARGET_ORDER:
        rows = [_minimal_row(target, bucket, f"{target}_{idx}") for idx, bucket in enumerate(TARGET_REQUIREMENTS[target], start=1)]
        rules: dict[str, Any] = {"native_metrics_before_selection_allowed": False}
        if target == "XCL1_lymphotactin":
            rules["mixed_state_pooling_allowed"] = False
        if target == "alpha_synuclein_SNCA":
            rules["single_native_fold_model_allowed"] = False
            rules["ensemble_not_single_fold"] = True
        packages[target] = {
            "source_manifest": {"target": target, "rows": rows},
            "evidence_dossier": {"target": target, "grammar_rules": rules},
        }
    return packages


def _control(
    control_id: str,
    observed: dict[str, Any],
    expected_status: str,
    expected_ready: list[str] | None,
    expected_partial_or_blocked: str,
    reason: str,
) -> dict[str, Any]:
    passed = observed.get("control_status") == expected_status
    if expected_ready is not None:
        passed = passed and sorted(observed.get("targets_ready", [])) == sorted(expected_ready)
    if expected_partial_or_blocked:
        checks = " ".join(str(item) for item in observed.get("failed_checks", []))
        passed = passed and expected_partial_or_blocked in checks
    return {
        "control_id": control_id,
        "expected_status": expected_status,
        "observed_status": observed.get("control_status"),
        "expected_ready": expected_ready,
        "observed_ready": observed.get("targets_ready"),
        "expected_failure_contains": expected_partial_or_blocked,
        "observed_failed_checks": observed.get("failed_checks"),
        "passed": passed,
        "reason": reason,
    }


def _evaluate_control_packages(packages: dict[str, dict[str, Any]]) -> dict[str, Any]:
    target_results = [evaluate_target(target, packages.get(target, {})) for target in TARGET_ORDER]
    targets_ready = [row["target"] for row in target_results if row["target_status"] == "ready"]
    targets_partial = [row["target"] for row in target_results if row["target_status"] == "partial"]
    targets_blocked = [row["target"] for row in target_results if row["target_status"] == "blocked"]
    coordinate_count = sum(int(row["coordinate_derived_source_count"]) for row in target_results)
    internal_count = sum(int(row["internal_runtime_source_count"]) for row in target_results)
    native_metric_count = sum(int(row["native_metric_source_count"]) for row in target_results)
    placeholder_count = sum(int(row["placeholder_source_count"]) for row in target_results)
    failed_checks: list[str] = []
    for row in target_results:
        failed_checks.extend(row["failed_checks"])
    if coordinate_count or internal_count or native_metric_count:
        status = BLOCKED_LEAKAGE_STATUS
    elif placeholder_count or targets_blocked:
        status = BLOCKED_UNTRUSTED_STATUS
    elif len(targets_ready) == len(TARGET_ORDER):
        status = READY_STATUS
    else:
        status = PARTIAL_STATUS
    return {
        "control_status": status,
        "targets_ready": targets_ready,
        "targets_partial": targets_partial,
        "targets_blocked": targets_blocked,
        "failed_checks": sorted(set(failed_checks)),
    }


def run_v36_controls() -> list[dict[str, Any]]:
    controls: list[dict[str, Any]] = []

    packages = _ready_fixture_packages()
    bad = packages["KcsA"]["source_manifest"]["rows"][0]
    bad["source_type"] = "RCSB coordinate contacts"
    bad["source_name"] = "KcsA 1BL8 coordinate contact CSV"
    bad["source_url_or_citation"] = "data/external_constraints/KcsA/pore_filter/kcsa_1bl8_pore_filter_external_contacts.csv"
    bad["coordinate_derived"] = True
    controls.append(_control(
        "kcsa_v33_v34_coordinate_csvs_blocked",
        _evaluate_control_packages(packages),
        BLOCKED_LEAKAGE_STATUS,
        None,
        "coordinate_or_prediction_leakage",
        "KcsA V33/V34 coordinate-derived contact CSVs must be blocked as V36 evidence.",
    ))

    packages = _ready_fixture_packages()
    packages["KcsA"]["source_manifest"]["rows"] = [_minimal_row("KcsA", "sequence_or_family_identity", "generic_channel")]
    packages["KcsA"]["source_manifest"]["rows"][0]["evidence_statement"] = "generic channel annotation"
    controls.append(_control(
        "kcsa_generic_channel_annotation_partial",
        _evaluate_control_packages(packages),
        PARTIAL_STATUS,
        ["XCL1_lymphotactin", "alpha_synuclein_SNCA"],
        "KcsA:missing_bucket:ion_selectivity_context",
        "A generic channel label is not enough for KcsA readiness.",
    ))

    packages = _ready_fixture_packages()
    packages["XCL1_lymphotactin"]["source_manifest"]["rows"] = [
        row for row in packages["XCL1_lymphotactin"]["source_manifest"]["rows"]
        if row["evidence_bucket"] != "state_B_function_context"
    ]
    controls.append(_control(
        "xcl1_one_state_only_partial",
        _evaluate_control_packages(packages),
        PARTIAL_STATUS,
        ["KcsA", "alpha_synuclein_SNCA"],
        "XCL1_lymphotactin:missing_bucket:state_B_function_context",
        "XCL1 needs both native-state function contexts.",
    ))

    packages = _ready_fixture_packages()
    packages["XCL1_lymphotactin"]["evidence_dossier"]["grammar_rules"]["mixed_state_pooling_allowed"] = True
    controls.append(_control(
        "xcl1_mixed_state_pooling_blocked",
        _evaluate_control_packages(packages),
        BLOCKED_UNTRUSTED_STATUS,
        None,
        "mixed_state_pooling_forbidden",
        "XCL1 state-A/state-B pooling is a blocked grammar error.",
    ))

    packages = _ready_fixture_packages()
    packages["alpha_synuclein_SNCA"]["evidence_dossier"]["grammar_rules"]["single_native_fold_model_allowed"] = True
    controls.append(_control(
        "snca_single_fold_grammar_blocked",
        _evaluate_control_packages(packages),
        BLOCKED_UNTRUSTED_STATUS,
        None,
        "single_native_fold_grammar_forbidden_for_SNCA",
        "SNCA must remain an IDP/ensemble problem, not a forced single-fold problem.",
    ))

    packages = _ready_fixture_packages()
    packages["alpha_synuclein_SNCA"]["source_manifest"]["rows"] = [
        row for row in packages["alpha_synuclein_SNCA"]["source_manifest"]["rows"]
        if row["evidence_bucket"] != "intrinsic_disorder_context"
    ]
    controls.append(_control(
        "snca_without_disorder_evidence_partial",
        _evaluate_control_packages(packages),
        PARTIAL_STATUS,
        ["KcsA", "XCL1_lymphotactin"],
        "alpha_synuclein_SNCA:missing_bucket:intrinsic_disorder_context",
        "SNCA cannot be ready if disorder evidence is removed.",
    ))

    packages = _ready_fixture_packages()
    packages["KcsA"]["source_manifest"]["rows"][0]["source_name"] = "TODO placeholder"
    packages["KcsA"]["source_manifest"]["rows"][0]["source_url_or_citation"] = "citation needed"
    controls.append(_control(
        "placeholder_citation_blocked",
        _evaluate_control_packages(packages),
        BLOCKED_UNTRUSTED_STATUS,
        None,
        "placeholder_or_missing_source_name",
        "Placeholder source names or citations must be blocked.",
    ))

    packages = _ready_fixture_packages()
    packages["KcsA"]["source_manifest"]["rows"][0]["file_path"] = "first_contact_clean_pharmacotopology_layer_run/V35/report.json"
    packages["KcsA"]["source_manifest"]["rows"][0]["internal_runtime_source"] = True
    controls.append(_control(
        "internal_runtime_source_blocked",
        _evaluate_control_packages(packages),
        BLOCKED_LEAKAGE_STATUS,
        None,
        "internal_runtime_source_supplied",
        "Runtime reports under first_contact_clean_pharmacotopology_layer_run are not evidence.",
    ))

    packages = _ready_fixture_packages()
    packages["KcsA"]["source_manifest"]["rows"][0]["native_metrics_used_for_selection"] = True
    controls.append(_control(
        "native_coordinate_metrics_before_selection_blocked",
        _evaluate_control_packages(packages),
        BLOCKED_LEAKAGE_STATUS,
        None,
        "native_metrics_used_for_selection",
        "Native coordinate metrics before selection must be blocked.",
    ))

    return controls


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V36 Real Evidence Dossier Gate",
        "",
        f"Status: `{cert.get('control_status')}`",
        f"Targets ready: `{cert.get('targets_ready')}`",
        f"Targets partial: `{cert.get('targets_partial')}`",
        f"Targets blocked: `{cert.get('targets_blocked')}`",
        f"External source count: `{cert.get('external_source_count')}`",
        f"Coordinate-derived source count: `{cert.get('coordinate_derived_source_count')}`",
        f"Internal runtime source count: `{cert.get('internal_runtime_source_count')}`",
        f"Placeholder source count: `{cert.get('placeholder_source_count')}`",
        f"Controls: `{cert.get('passed_control_count')}` / `{cert.get('control_count')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD allowed: `{cert.get('new_MD_allowed')}`",
        f"Folding solved: `{cert.get('folding_problem_solved')}`",
        f"Next action: `{cert.get('next_action')}`",
        "",
        "## Source Counts By Target",
    ]
    for target, count in cert.get("source_counts_by_target", {}).items():
        lines.append(f"- `{target}`: `{count}`")
    lines.extend(["", "## Source Types By Target"])
    for target, types in cert.get("source_types_by_target", {}).items():
        lines.append(f"- `{target}`: `{types}`")
    lines.extend(["", "## Failed Checks"])
    failures = cert.get("failed_checks") or []
    if failures:
        lines.extend(f"- `{failure}`" for failure in failures)
    else:
        lines.append("- None")
    lines.extend(["", "## Controls"])
    for row in cert.get("controls", []):
        lines.extend([
            f"### {row.get('control_id')}",
            f"Passed: `{row.get('passed')}`",
            f"Observed status: `{row.get('observed_status')}`",
            f"Expected status: `{row.get('expected_status')}`",
            f"Reason: {row.get('reason')}",
            "",
        ])
    lines.extend(["## Locked Interpretation", str(cert.get("locked_interpretation", ""))])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v36_real_evidence_dossier_gate_certificate.json"
    report_path = out_dir / "V36_REAL_EVIDENCE_DOSSIER_GATE_REPORT.md"
    decision_path = out_dir / "v36_real_evidence_dossier_gate_next_decision.json"
    target_results_path = out_dir / "v36_target_dossier_evaluations.json"
    decision = {
        "kind": "V36_REAL_EVIDENCE_DOSSIER_GATE_NEXT_DECISION_v0",
        "decision_status": cert.get("control_status"),
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
            "target_results": str(target_results_path),
        },
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, decision)
    _write_json(target_results_path, cert.get("target_results", []))
    _write_report(report_path, cert)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path, "target_results": target_results_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V36 real evidence dossier gate.")
    parser.add_argument("--data-root", type=Path, default=DATA_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    packages = load_packages(args.data_root)
    cert = build_v36(packages)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert.get("kind"),
        "control_status": cert.get("control_status"),
        "target_count": cert.get("target_count"),
        "targets_ready": cert.get("targets_ready"),
        "targets_partial": cert.get("targets_partial"),
        "targets_blocked": cert.get("targets_blocked"),
        "external_source_count": cert.get("external_source_count"),
        "source_counts_by_target": cert.get("source_counts_by_target"),
        "source_types_by_target": cert.get("source_types_by_target"),
        "coordinate_derived_source_count": cert.get("coordinate_derived_source_count"),
        "internal_runtime_source_count": cert.get("internal_runtime_source_count"),
        "placeholder_source_count": cert.get("placeholder_source_count"),
        "control_count": cert.get("control_count"),
        "passed_control_count": cert.get("passed_control_count"),
        "claim_allowed": cert.get("claim_allowed"),
        "new_MD_allowed": cert.get("new_MD_allowed"),
        "new_md_executed": cert.get("new_md_executed"),
        "positive_folding_evidence_found": cert.get("positive_folding_evidence_found"),
        "folding_problem_solved": cert.get("folding_problem_solved"),
        "next_action": cert.get("next_action"),
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
