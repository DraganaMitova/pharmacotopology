#!/usr/bin/env python3
from __future__ import annotations

"""V35 non-coordinate evolutionary / blind holdout gate.

V35 is deliberately strict. It does not import coordinate-derived contacts,
does not run MD, and does not make a folding claim. Its job is to decide
whether a real external, non-coordinate evolutionary source exists locally and
contains enough potassium-channel operator context to open the next
claim-disabled holdout step.
"""

import argparse
import json
import os
import re
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"

DEFAULT_V34_CERT = RUN_ROOT / "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS" / "v34_kcsa_discriminative_content_controls_certificate.json"
DEFAULT_MANIFEST = REPO_ROOT / "data" / "external_constraints" / "v35_noncoordinate_external_source_manifest.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V35_NONCOORDINATE_EVOLUTIONARY_HOLDOUT"

KCSA_NONCOORDINATE_ROOT = REPO_ROOT / "data" / "external_constraints" / "KcsA" / "noncoordinate_evolutionary"

CLEAN_ABSTAIN = "V35_CLEAN_ABSTAIN_NO_NONCOORDINATE_EXTERNAL_EVIDENCE"
BLOCKED_COORDINATE = "V35_BLOCKED_COORDINATE_DERIVED_SOURCE_SUPPLIED"
BLOCKED_INTERNAL = "V35_BLOCKED_INTERNAL_RUNTIME_SOURCE_SUPPLIED"
READY = "V35_NONCOORDINATE_EVOLUTIONARY_HOLDOUT_READY_CLAIM_DISABLED"
PASSED = "V35_NONCOORDINATE_EVOLUTIONARY_CONTENT_CONTROLS_PASSED_CLAIM_DISABLED"

V34_READY_STATUS = "V34_KCSA_DISCRIMINATIVE_CONTENT_CONTROLS_PASSED_CLAIM_DISABLED"

REQUIRED_ROW_FIELDS = [
    "target",
    "file_path",
    "evidence_type",
    "source_name",
    "source_url_or_citation",
    "source_date_or_version",
    "source_boundary",
    "allowed_use",
    "provenance_notes",
    "coordinate_derived",
    "native_coordinate_used_before_selection",
    "claim_allowed",
]

VALID_EVIDENCE_TYPES = {
    "sequence_conservation_filter_signature",
    "external_msa_family_context",
    "evolutionary_coupling_candidate",
    "potassium_channel_family_signature",
}

EXCLUDED_EVIDENCE_TYPES = {
    "pdb_coordinate_contact",
    "coordinate_derived_contact",
    "internal_runtime_report",
    "annotation_only_claim",
    "generated_report_summary",
}

COORDINATE_EVIDENCE_TYPES = {"pdb_coordinate_contact", "coordinate_derived_contact"}
INTERNAL_EVIDENCE_TYPES = {"internal_runtime_report", "generated_report_summary"}

PLACEHOLDER_TOKENS = {
    "placeholder",
    "todo",
    "tbd",
    "unknown",
    "citation needed",
    "example.com",
    "dummy",
}

COORDINATE_TOKENS = {
    "1bl8",
    "rcsb pdb",
    "pdb 1bl8",
    ".pdb",
    "raw_sources",
    "coordinate-derived",
    "coordinate derived",
    "coordinate_contact",
}


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        if path == DEFAULT_MANIFEST:
            return {"kind": "V35_NONCOORDINATE_EXTERNAL_SOURCE_MANIFEST_v0", "rows": []}
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


def _norm_path_text(value: Any) -> str:
    return str(value or "").replace("\\", "/").strip()


def _resolve_source_path(file_path: Any) -> Path:
    text = str(file_path or "").strip()
    path = Path(text)
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path.resolve(strict=False)


def _is_under(path: Path, root: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(root.resolve(strict=False))
        return True
    except ValueError:
        return False


def _contains_placeholder(value: Any) -> bool:
    text = str(value or "").strip().lower()
    if not text:
        return True
    return any(token in text for token in PLACEHOLDER_TOKENS)


def _combined_row_text(row: dict[str, Any]) -> str:
    return " ".join(str(row.get(key, "")) for key in REQUIRED_ROW_FIELDS).lower()


def _row_is_coordinate_derived(row: dict[str, Any], norm_path: str) -> bool:
    evidence_type = str(row.get("evidence_type", "")).strip()
    combined = f"{norm_path} {_combined_row_text(row)}".lower()
    if evidence_type in COORDINATE_EVIDENCE_TYPES:
        return True
    if _bool_true(row.get("coordinate_derived")):
        return True
    if _bool_true(row.get("native_coordinate_used_before_selection")):
        return True
    return any(token in combined for token in COORDINATE_TOKENS)


def _row_is_internal_runtime(row: dict[str, Any], norm_path: str, resolved_path: Path) -> bool:
    evidence_type = str(row.get("evidence_type", "")).strip()
    combined = f"{norm_path} {_combined_row_text(row)}".lower()
    return (
        evidence_type in INTERNAL_EVIDENCE_TYPES
        or "first_contact_clean_pharmacotopology_layer_run" in combined
        or _is_under(resolved_path, RUN_ROOT)
    )


def _read_source_text(
    norm_path: str,
    resolved_path: Path,
    content_by_path: dict[str, str] | None,
) -> tuple[str | None, str | None]:
    if content_by_path is not None:
        keys = [norm_path, str(resolved_path)]
        try:
            keys.append(str(resolved_path.relative_to(REPO_ROOT)))
        except ValueError:
            pass
        for key in keys:
            if key in content_by_path:
                return content_by_path[key], None
    if not resolved_path.exists():
        return None, "source_file_missing"
    if resolved_path.is_dir():
        return None, "source_path_is_directory"
    try:
        return resolved_path.read_text(encoding="utf-8"), None
    except UnicodeDecodeError:
        return None, "source_file_not_utf8_text"


def _content_coordinate_like(text: str) -> bool:
    upper = text.upper()
    if re.search(r"(?m)^(ATOM|HETATM)\s+\d+", upper):
        return True
    coordinate_terms = ["CARTESIAN", "ANGSTROM", "COORDINATE-DERIVED", "RCSB PDB", "PDB 1BL8", "1BL8"]
    return any(term in upper for term in coordinate_terms)


def validate_noncoordinate_content(text: str, evidence_type: str) -> dict[str, Any]:
    failures: list[str] = []
    upper = text.upper()
    compact = re.sub(r"[^A-Z0-9+]+", "", upper)

    has_filter_signature = any(token in compact for token in ["TVGYG", "TIGYG", "TXGYG", "TXXTXGYG"])
    has_potassium_specificity = (
        "K+" in upper
        or "POTASSIUM" in upper
        or "K-SELECTIVE" in upper
        or "K_SELECTIVE" in upper
        or re.search(r"\bION[_\s-]*SPECIFICITY\s*[:=,]\s*K\b", upper) is not None
    )
    has_family_specificity = any(
        token in upper
        for token in ["POTASSIUM CHANNEL", "K CHANNEL", "K-CHANNEL", "K_CHANNEL", "KCSA", "PFAM", "INTERPRO"]
    )
    has_evolutionary_basis = any(
        token in upper
        for token in ["MSA", "ALIGNMENT", "CONSERVED", "CONSERVATION", "HOMOLOG", "FAMILY", "PFAM", "INTERPRO", "COUPLING", "COEVOLUTION", "EVOLUTIONARY"]
    )
    has_coupling_specificity = evidence_type != "evolutionary_coupling_candidate" or any(
        token in upper for token in ["COUPLING", "COEVOLUTION", "DIRECT COUPLING", "EVOLUTIONARY COUPLING"]
    )

    if not text.strip():
        failures.append("source_content_empty")
    if _content_coordinate_like(text):
        failures.append("coordinate_like_content_detected")
    if not has_filter_signature:
        failures.append("missing_potassium_channel_filter_signature")
    if not has_potassium_specificity:
        failures.append("missing_potassium_ion_specificity")
    if not has_family_specificity:
        failures.append("missing_potassium_channel_family_specificity")
    if not has_evolutionary_basis:
        failures.append("missing_external_evolutionary_or_family_basis")
    if not has_coupling_specificity:
        failures.append("missing_evolutionary_coupling_specificity")

    return {
        "valid": not failures,
        "failures": sorted(set(failures)),
        "summary": {
            "has_filter_signature": has_filter_signature,
            "has_potassium_specificity": bool(has_potassium_specificity),
            "has_family_specificity": has_family_specificity,
            "has_evolutionary_basis": has_evolutionary_basis,
            "has_coupling_specificity": has_coupling_specificity,
            "coordinate_like_content_detected": _content_coordinate_like(text),
        },
    }


def _row_preflight_failures(row: dict[str, Any], norm_path: str, resolved_path: Path) -> list[str]:
    failures: list[str] = []
    missing = [field for field in REQUIRED_ROW_FIELDS if field not in row]
    if missing:
        failures.append("missing_required_manifest_fields:" + ",".join(missing))

    if row.get("target") != "KcsA":
        failures.append("target_not_KcsA")

    evidence_type = str(row.get("evidence_type", "")).strip()
    if evidence_type in EXCLUDED_EVIDENCE_TYPES:
        failures.append(f"excluded_evidence_type:{evidence_type}")
    elif evidence_type not in VALID_EVIDENCE_TYPES:
        failures.append(f"unsupported_evidence_type:{evidence_type}")

    if not _bool_false(row.get("coordinate_derived")):
        failures.append("coordinate_derived_field_not_false")
    if not _bool_false(row.get("native_coordinate_used_before_selection")):
        failures.append("native_coordinate_used_before_selection_field_not_false")
    if not _bool_false(row.get("claim_allowed")):
        failures.append("claim_allowed_field_not_false")

    for key in ["source_name", "source_url_or_citation", "source_date_or_version"]:
        if _contains_placeholder(row.get(key)):
            failures.append(f"placeholder_or_missing_{key}")

    if len(str(row.get("source_url_or_citation", "")).strip()) < 12:
        failures.append("source_url_or_citation_not_stable_enough")
    if not str(row.get("source_boundary", "")).strip():
        failures.append("source_boundary_missing")
    if not str(row.get("allowed_use", "")).strip():
        failures.append("allowed_use_missing")
    if not str(row.get("provenance_notes", "")).strip():
        failures.append("provenance_notes_missing")

    if not _is_under(resolved_path, KCSA_NONCOORDINATE_ROOT):
        failures.append("source_path_not_in_kcsa_noncoordinate_evolutionary_boundary")
    if not norm_path:
        failures.append("file_path_missing")

    return failures


def evaluate_manifest_rows(
    rows: list[Any],
    content_by_path: dict[str, str] | None = None,
) -> dict[str, Any]:
    failed_checks: list[str] = []
    evaluated_rows: list[dict[str, Any]] = []
    coordinate_count = 0
    internal_count = 0
    noncoordinate_count = 0
    operator_candidate_found = False
    selected_sources: list[str] = []

    for index, raw_row in enumerate(rows):
        if not isinstance(raw_row, dict):
            failed_checks.append(f"row_{index}:manifest_row_not_object")
            evaluated_rows.append({"row_index": index, "valid_noncoordinate_source": False, "failures": ["manifest_row_not_object"]})
            continue

        norm_path = _norm_path_text(raw_row.get("file_path"))
        resolved_path = _resolve_source_path(norm_path)
        row_failures = _row_preflight_failures(raw_row, norm_path, resolved_path)

        internal = _row_is_internal_runtime(raw_row, norm_path, resolved_path)
        coordinate = _row_is_coordinate_derived(raw_row, norm_path)
        if internal:
            internal_count += 1
            row_failures.append("internal_runtime_source_supplied")
        if coordinate:
            coordinate_count += 1
            row_failures.append("coordinate_derived_source_supplied")

        content_result: dict[str, Any] | None = None
        if not internal and not coordinate and not row_failures:
            source_text, read_failure = _read_source_text(norm_path, resolved_path, content_by_path)
            if read_failure:
                row_failures.append(read_failure)
            else:
                content_result = validate_noncoordinate_content(source_text or "", str(raw_row.get("evidence_type", "")))
                if "coordinate_like_content_detected" in content_result.get("failures", []):
                    coordinate_count += 1
                    coordinate = True
                    row_failures.append("coordinate_like_content_detected")
                elif content_result.get("valid") is True:
                    noncoordinate_count += 1
                    operator_candidate_found = True
                    selected_sources.append(norm_path)
                else:
                    row_failures.extend(str(item) for item in _as_list(content_result.get("failures")))

        failed_checks.extend(f"row_{index}:{failure}" for failure in sorted(set(row_failures)))
        evaluated_rows.append({
            "row_index": index,
            "target": raw_row.get("target"),
            "file_path": norm_path,
            "evidence_type": raw_row.get("evidence_type"),
            "source_name": raw_row.get("source_name"),
            "internal_runtime_source": internal,
            "coordinate_derived_source": coordinate,
            "valid_noncoordinate_source": not row_failures and content_result is not None and content_result.get("valid") is True,
            "content_summary": content_result.get("summary") if content_result else None,
            "failures": sorted(set(row_failures)),
        })

    if coordinate_count:
        source_status = BLOCKED_COORDINATE
    elif internal_count:
        source_status = BLOCKED_INTERNAL
    elif noncoordinate_count:
        source_status = READY
    else:
        source_status = CLEAN_ABSTAIN

    return {
        "source_status": source_status,
        "evidence_row_count": len(rows),
        "noncoordinate_external_source_count": noncoordinate_count,
        "coordinate_derived_source_count": coordinate_count,
        "internal_runtime_source_count": internal_count,
        "operator_candidate_found": operator_candidate_found,
        "selected_source_paths": selected_sources,
        "evaluated_rows": evaluated_rows,
        "failed_checks": sorted(set(failed_checks)),
    }


def _valid_control_row(file_path: str) -> dict[str, Any]:
    return {
        "target": "KcsA",
        "file_path": file_path,
        "evidence_type": "sequence_conservation_filter_signature",
        "source_name": "External potassium channel family conservation source",
        "source_url_or_citation": "UniProt/Pfam-style external family conservation source, accessioned local holdout fixture",
        "source_date_or_version": "fixture_2026-06-11",
        "source_boundary": "external_noncoordinate_sequence_family_evidence_no_pdb_coordinates",
        "allowed_use": "claim_disabled_V35_noncoordinate_holdout_gate_only",
        "provenance_notes": "Sequence-family conservation context only; no PDB coordinates, no native metrics, no MD.",
        "coordinate_derived": False,
        "native_coordinate_used_before_selection": False,
        "claim_allowed": False,
    }


def _valid_control_content() -> str:
    return "\n".join([
        "target,family_context,filter_signature,ion_specificity,evolutionary_basis,notes",
        "KcsA,KcsA-like potassium channel family,TVGYG,K+,external MSA conservation across homolog family,sequence-only evidence no coordinates",
        "",
    ])


def _control_result(
    control_id: str,
    observed: dict[str, Any],
    expected_status: str,
    expected_failure_any: list[str],
    reason: str,
) -> dict[str, Any]:
    failed = " ".join(str(item) for item in _as_list(observed.get("failed_checks")))
    expected_hit = not expected_failure_any or any(item in failed for item in expected_failure_any)
    passed = observed.get("source_status") == expected_status and expected_hit
    return {
        "control_id": control_id,
        "expected_status": expected_status,
        "observed_status": observed.get("source_status"),
        "expected_failure_any": expected_failure_any,
        "observed_failed_checks": observed.get("failed_checks"),
        "passed": passed,
        "reason": reason,
    }


def run_v35_controls() -> list[dict[str, Any]]:
    path = "data/external_constraints/KcsA/noncoordinate_evolutionary/control_signature.csv"
    good_row = _valid_control_row(path)
    good_content = _valid_control_content()

    controls: list[dict[str, Any]] = []

    controls.append(_control_result(
        "missing_noncoordinate_source_clean_abstains",
        evaluate_manifest_rows([]),
        CLEAN_ABSTAIN,
        [],
        "No local non-coordinate external evidence must produce a clean abstain, not a fake pass.",
    ))

    coordinate_row = {
        **good_row,
        "file_path": "data/external_constraints/KcsA/pore_filter/kcsa_1bl8_pore_filter_external_contacts.csv",
        "evidence_type": "coordinate_derived_contact",
        "source_name": "RCSB PDB 1BL8 coordinate-derived contact CSV",
        "source_url_or_citation": "RCSB PDB 1BL8 coordinate-derived contacts",
        "coordinate_derived": True,
    }
    controls.append(_control_result(
        "coordinate_derived_kcsa_csv_is_blocked",
        evaluate_manifest_rows([coordinate_row]),
        BLOCKED_COORDINATE,
        ["coordinate_derived_source_supplied"],
        "Coordinate-derived KcsA CSVs from V33/V34 must not open V35.",
    ))

    runtime_row = {
        **good_row,
        "file_path": "first_contact_clean_pharmacotopology_layer_run/V33_NEGATIVE_CONTROLS/internal_poison_source.csv",
        "evidence_type": "internal_runtime_report",
        "source_name": "Internal runtime report",
        "source_url_or_citation": "first_contact_clean_pharmacotopology_layer_run generated runtime report",
    }
    controls.append(_control_result(
        "internal_runtime_source_is_blocked",
        evaluate_manifest_rows([runtime_row]),
        BLOCKED_INTERNAL,
        ["internal_runtime_source_supplied"],
        "Generated runtime artifacts are audit outputs, not external evidence.",
    ))

    placeholder_row = {
        **good_row,
        "source_name": "TODO placeholder",
        "source_url_or_citation": "citation needed",
    }
    controls.append(_control_result(
        "placeholder_source_name_or_citation_is_invalid",
        evaluate_manifest_rows([placeholder_row], {path: good_content}),
        CLEAN_ABSTAIN,
        ["placeholder_or_missing_source_name", "placeholder_or_missing_source_url_or_citation"],
        "A row without stable provenance must not open V35 even if its content looks plausible.",
    ))

    annotation_row = {
        **good_row,
        "evidence_type": "annotation_only_claim",
    }
    controls.append(_control_result(
        "annotation_only_source_does_not_open_v35",
        evaluate_manifest_rows([annotation_row], {path: good_content}),
        CLEAN_ABSTAIN,
        ["excluded_evidence_type:annotation_only_claim"],
        "Annotation-only rows remain context-only and cannot become non-coordinate evidence.",
    ))

    wrong_target_row = {
        **good_row,
        "target": "XCL1_lymphotactin",
    }
    controls.append(_control_result(
        "wrong_target_in_kcsa_v35_path_does_not_open",
        evaluate_manifest_rows([wrong_target_row], {path: good_content}),
        CLEAN_ABSTAIN,
        ["target_not_KcsA"],
        "V35 v0 is scoped to KcsA; a wrong target in the KcsA V35 boundary must not pass.",
    ))

    renamed_signature = good_content.replace("TVGYG", "GYAVT")
    controls.append(_control_result(
        "renamed_potassium_filter_signature_fails_content_validation",
        evaluate_manifest_rows([good_row], {path: renamed_signature}),
        CLEAN_ABSTAIN,
        ["missing_potassium_channel_filter_signature"],
        "A shuffled or renamed filter signature must fail sequence-content validation.",
    ))

    sodium_content = good_content.replace("K+", "Na+").replace("potassium channel", "sodium channel")
    controls.append(_control_result(
        "ion_specificity_relabel_away_from_potassium_fails",
        evaluate_manifest_rows([good_row], {path: sodium_content}),
        CLEAN_ABSTAIN,
        ["missing_potassium_ion_specificity"],
        "KcsA-like non-coordinate context must not survive relabeling away from K+.",
    ))

    generic_content = "\n".join([
        "target,annotation,notes",
        "protein,channel protein,generic label only",
        "",
    ])
    controls.append(_control_result(
        "generic_channel_annotation_is_not_enough",
        evaluate_manifest_rows([good_row], {path: generic_content}),
        CLEAN_ABSTAIN,
        ["missing_potassium_channel_filter_signature", "missing_potassium_channel_family_specificity"],
        "A generic channel label is not enough without filter, ion, and evolutionary/family specificity.",
    ))

    return controls


def _v34_precondition_failures(v34: dict[str, Any] | None) -> list[str]:
    if v34 is None:
        return []
    failures: list[str] = []
    if v34.get("control_status") != V34_READY_STATUS:
        failures.append("V34_control_status_not_passed_for_V35")
    for key in ["claim_allowed", "new_MD_allowed", "positive_folding_evidence_found", "folding_problem_solved"]:
        if v34.get(key) is not False:
            failures.append(f"V34_{key}_must_be_false")
    if v34.get("next_action") not in {None, "move_to_non_coordinate_external_evolutionary_or_blind_holdout_tests"}:
        failures.append("V34_next_action_not_V35")
    return sorted(set(failures))


def build_v35(
    manifest: dict[str, Any],
    v34: dict[str, Any] | None = None,
    content_by_path: dict[str, str] | None = None,
) -> dict[str, Any]:
    rows = _as_list(manifest.get("rows"))
    evaluated = evaluate_manifest_rows(rows, content_by_path)
    controls = run_v35_controls()
    precondition_failures = _v34_precondition_failures(v34)

    all_controls_passed = all(row.get("passed") is True for row in controls)
    source_status = evaluated["source_status"]
    if precondition_failures and source_status == READY:
        control_status = READY
    elif source_status == READY and all_controls_passed:
        control_status = PASSED
    else:
        control_status = source_status

    failed_checks = list(evaluated.get("failed_checks", []))
    failed_checks.extend(precondition_failures)
    if not all_controls_passed:
        failed_checks.append("V35_required_negative_controls_failed")

    if precondition_failures:
        next_action = "fix_V34_preconditions_before_V35"
    elif control_status == PASSED:
        next_action = "move_to_blind_holdout_operator_test_no_claim_no_MD"
    elif control_status == BLOCKED_COORDINATE:
        next_action = "remove_coordinate_derived_source_and_rerun_V35"
    elif control_status == BLOCKED_INTERNAL:
        next_action = "remove_internal_runtime_source_and_rerun_V35"
    elif control_status == READY:
        next_action = "fix_V35_content_controls_before_any_claim_or_MD"
    else:
        next_action = "acquire_real_noncoordinate_external_evolutionary_source_before_any_claim_or_MD"

    return {
        "kind": "V35_NONCOORDINATE_EVOLUTIONARY_OR_BLIND_HOLDOUT_TESTS_v0",
        "run_mode": "noncoordinate_external_evolutionary_holdout_no_coordinates_no_MD_claim_disabled",
        "source_boundary": "external_noncoordinate_evolutionary_or_blind_holdout_evidence_only_no_pdb_coordinates_no_internal_runtime_reports",
        "selected_target": "KcsA",
        "source_status": source_status,
        "evidence_row_count": evaluated["evidence_row_count"],
        "noncoordinate_external_source_count": evaluated["noncoordinate_external_source_count"],
        "coordinate_derived_source_count": evaluated["coordinate_derived_source_count"],
        "internal_runtime_source_count": evaluated["internal_runtime_source_count"],
        "operator_candidate_found": evaluated["operator_candidate_found"],
        "coordinate_truth_used_before_selection": False,
        "native_coordinate_used_before_selection": False,
        "native_metrics_used_for_selection": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
        "control_status": control_status,
        "passed_control_count": sum(1 for row in controls if row.get("passed") is True),
        "control_count": len(controls),
        "controls": controls,
        "failed_checks": sorted(set(str(item) for item in failed_checks)),
        "precondition_failures": precondition_failures,
        "evaluated_rows": evaluated["evaluated_rows"],
        "selected_source_paths": evaluated["selected_source_paths"],
        "next_action": next_action,
        "locked_interpretation": (
            "V35 only opens on real local non-coordinate external evolutionary evidence with stable provenance and KcsA-like potassium-channel sequence/family content. "
            "Coordinate-derived contacts, internal runtime reports, annotation-only rows, placeholders, wrong targets, and generic channel labels remain blocked or clean-abstained. "
            "Passing V35 is not a folding claim, does not run MD, and does not solve protein folding."
        ),
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V35 Non-Coordinate Evolutionary Holdout",
        "",
        f"Status: `{cert.get('control_status')}`",
        f"Source status: `{cert.get('source_status')}`",
        f"Selected target: `{cert.get('selected_target')}`",
        f"Evidence rows: `{cert.get('evidence_row_count')}`",
        f"Non-coordinate external sources: `{cert.get('noncoordinate_external_source_count')}`",
        f"Coordinate-derived sources: `{cert.get('coordinate_derived_source_count')}`",
        f"Internal runtime sources: `{cert.get('internal_runtime_source_count')}`",
        f"Operator candidate found: `{cert.get('operator_candidate_found')}`",
        f"Controls: `{cert.get('passed_control_count')}` / `{cert.get('control_count')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD allowed: `{cert.get('new_MD_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        f"Positive folding evidence found: `{cert.get('positive_folding_evidence_found')}`",
        f"Folding solved: `{cert.get('folding_problem_solved')}`",
        f"Next action: `{cert.get('next_action')}`",
        "",
        "## Selected Sources",
    ]
    selected = _as_list(cert.get("selected_source_paths"))
    if selected:
        lines.extend(f"- `{path}`" for path in selected)
    else:
        lines.append("- None")

    lines.extend(["", "## Failed Checks"])
    failures = _as_list(cert.get("failed_checks"))
    if failures:
        lines.extend(f"- `{failure}`" for failure in failures)
    else:
        lines.append("- None")

    lines.extend(["", "## Controls"])
    for row in _as_list(cert.get("controls")):
        lines.extend([
            f"### {row.get('control_id')}",
            f"Passed: `{row.get('passed')}`",
            f"Observed status: `{row.get('observed_status')}`",
            f"Expected status: `{row.get('expected_status')}`",
            f"Observed failed checks: `{row.get('observed_failed_checks')}`",
            f"Reason: {row.get('reason')}",
            "",
        ])

    lines.extend([
        "## Locked Interpretation",
        str(cert.get("locked_interpretation", "")),
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v35_noncoordinate_evolutionary_holdout_certificate.json"
    report_path = out_dir / "V35_NONCOORDINATE_EVOLUTIONARY_HOLDOUT_REPORT.md"
    decision_path = out_dir / "v35_noncoordinate_evolutionary_holdout_next_decision.json"
    rows_path = out_dir / "v35_evaluated_source_rows.json"

    decision = {
        "kind": "V35_NONCOORDINATE_EVOLUTIONARY_HOLDOUT_NEXT_DECISION_v0",
        "decision_status": cert.get("control_status"),
        "source_status": cert.get("source_status"),
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
            "evaluated_rows": str(rows_path),
        },
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, decision)
    _write_json(rows_path, cert.get("evaluated_rows", []))
    _write_report(report_path, cert)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path, "evaluated_rows": rows_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V35 non-coordinate evolutionary holdout gate.")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--v34-certificate", type=Path, default=DEFAULT_V34_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    manifest = _read_json(args.manifest, "V35 non-coordinate external source manifest")
    v34 = _read_json(args.v34_certificate, "V34 certificate") if args.v34_certificate.exists() else None
    cert = build_v35(manifest, v34)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert.get("kind"),
        "control_status": cert.get("control_status"),
        "source_status": cert.get("source_status"),
        "evidence_row_count": cert.get("evidence_row_count"),
        "noncoordinate_external_source_count": cert.get("noncoordinate_external_source_count"),
        "coordinate_derived_source_count": cert.get("coordinate_derived_source_count"),
        "internal_runtime_source_count": cert.get("internal_runtime_source_count"),
        "operator_candidate_found": cert.get("operator_candidate_found"),
        "passed_control_count": cert.get("passed_control_count"),
        "control_count": cert.get("control_count"),
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
