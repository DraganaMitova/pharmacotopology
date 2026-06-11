#!/usr/bin/env python3
from __future__ import annotations

"""Run V38 blind mechanism generalization panel."""

import argparse
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
PANEL_ROOT = REPO_ROOT / "data" / "blind_mechanism_panels" / "V38"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V38_BLIND_MECHANISM_GENERALIZATION_PANEL"

PASSED = "V38_BLIND_MECHANISM_GENERALIZATION_PASSED_CLAIM_DISABLED"
PARTIAL = "V38_PARTIAL_GENERALIZATION_CLEAN_ABSTAIN"
BLOCKED_KEY = "V38_BLOCKED_TARGET_NAME_OR_ANSWER_KEY_LEAKAGE"
BLOCKED_LEAKAGE = "V38_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"
FAILED_FP = "V38_FAILED_FALSE_POSITIVE_HARD_GRAMMAR"

HARD_CLASSES = {
    "membrane_pore_filter_oligomeric_ion_selectivity",
    "metamorphic_two_state_fold_switch",
    "intrinsic_disorder_contextual_ensemble",
}

SUPPORTED_CLASSES = HARD_CLASSES | {
    "other_membrane_or_transport_context",
    "soluble_single_or_contextual_fold_not_metamorphic",
    "other_disorder_or_low_confidence_ensemble",
    "insufficient_evidence_clean_abstain",
}

TARGET_NAME_TOKENS = {"kcsa", "xcl1", "lymphotactin", "alpha-synuclein", "snca", "bacteriorhodopsin", "aqp1", "cxcl8", "ubiquitin", "lysozyme", "myoglobin"}
COORDINATE_TOKENS = {"1bl8", "coordinate-derived", "coordinate_derived", "pdb coordinate", "alphafold", "esmfold", "rosettafold", "rmsd", "tm-score", ".pdb"}
PLACEHOLDER_TOKENS = {"placeholder", "todo", "tbd", "citation needed", "example.com", "dummy", "unknown"}


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


def _text(dossier: dict[str, Any]) -> str:
    return " ".join(
        str(piece or "")
        for piece in [
            " ".join(str(x) for x in _as_list(dossier.get("evidence_buckets"))),
            " ".join(str(x) for x in _as_list(dossier.get("mechanism_evidence"))),
            json.dumps(dossier.get("grammar_rules", {}), sort_keys=True),
            str(dossier.get("masked_id", "")),
        ]
    ).lower()


def _provenance_text(dossier: dict[str, Any]) -> str:
    values: list[str] = []
    for row in _as_list(dossier.get("source_provenance")):
        if isinstance(row, dict):
            for key in ["source_id", "source_type", "source_name", "source_url_or_citation", "evidence_statement"]:
                values.append(str(row.get(key, "")))
    return " ".join(values).lower()


def _leakage(dossier: dict[str, Any], require_masked: bool = True) -> dict[str, Any]:
    text = _text(dossier)
    provenance = _provenance_text(dossier)
    all_text = f"{text} {provenance}"
    failed: list[str] = []
    coord = 0
    internal = 0
    placeholder = 0
    key = 0
    if require_masked and any(token in text for token in TARGET_NAME_TOKENS):
        key += 1
        failed.append("target_name_token_in_decision_fields")
    if any(key_name in dossier for key_name in ["answer_key", "expected_mechanism_class", "expected_mechanism_class_for_answer_key_only", "target"]):
        key += 1
        failed.append("answer_key_or_unmasked_target_field_in_assignment_input")
    if _bool_true(dossier.get("answer_key_used_for_assignment")):
        key += 1
        failed.append("answer_key_used_for_assignment")
    if _bool_true(dossier.get("coordinate_truth_used_before_selection")) or _bool_true(dossier.get("native_metrics_used_for_selection")):
        coord += 1
        failed.append("native_coordinate_or_metric_used_before_selection")
    if any(token in all_text for token in COORDINATE_TOKENS):
        coord += 1
        failed.append("coordinate_or_predicted_structure_leakage")
    if "first_contact_clean_pharmacotopology_layer_run" in all_text:
        internal += 1
        failed.append("internal_runtime_source_supplied")
    for row in _as_list(dossier.get("source_provenance")):
        if not isinstance(row, dict):
            continue
        if _bool_true(row.get("coordinate_derived")) or _bool_true(row.get("coordinate_truth_used_before_selection")) or _bool_true(row.get("native_metrics_used_for_selection")):
            coord += 1
            failed.append("source_row_coordinate_or_native_metric_leakage")
        if _bool_true(row.get("internal_runtime_source")):
            internal += 1
            failed.append("source_row_internal_runtime_source")
        for key_name in ["source_name", "source_url_or_citation", "source_date_or_version"]:
            if key_name in row:
                val = str(row.get(key_name) or "").lower()
                if not val or any(token in val for token in PLACEHOLDER_TOKENS):
                    placeholder += 1
                    failed.append(f"source_row_placeholder_{key_name}")
    return {
        "coordinate_derived_source_count": coord,
        "internal_runtime_source_count": internal,
        "placeholder_source_count": placeholder,
        "answer_key_or_target_name_leakage_count": key,
        "failed_checks": sorted(set(failed)),
    }


def assign_masked_dossier(dossier: dict[str, Any], require_masked: bool = True) -> dict[str, Any]:
    text = _text(dossier)
    buckets = set(str(x) for x in _as_list(dossier.get("evidence_buckets")))
    rules = dossier.get("grammar_rules") if isinstance(dossier.get("grammar_rules"), dict) else {}
    leak = _leakage(dossier, require_masked=require_masked)
    failed = list(leak["failed_checks"])

    if leak["coordinate_derived_source_count"] or leak["internal_runtime_source_count"]:
        return _assignment(dossier, "insufficient_evidence_clean_abstain", 0, "blocked_leakage", "coordinate_or_internal_leakage", leak, failed)
    if leak["answer_key_or_target_name_leakage_count"]:
        return _assignment(dossier, "insufficient_evidence_clean_abstain", 0, "blocked_key", "target_name_or_answer_key_leakage", leak, failed)
    if leak["placeholder_source_count"]:
        return _assignment(dossier, "insufficient_evidence_clean_abstain", 0, "blocked_key", "placeholder_or_untrusted_source", leak, failed)

    kcsa_ops = {"membrane_topology_context", "ion_selectivity_context", "filter_or_signature_context"} <= buckets and ("k+" in text or "potassium" in text) and ("filter" in text or "tvgyg" in text)
    xcl1_ops = {"state_A_function_context", "state_B_function_context", "metamorphic_two_state_context", "no_mixed_state_pooling_rule"} <= buckets and rules.get("mixed_state_pooling_allowed") is False
    snca_ops = {"intrinsic_disorder_context", "ensemble_context", "disorder_to_order_context", "no_single_native_fold_rule"} <= buckets and rules.get("single_native_fold_model_allowed") is False

    if kcsa_ops:
        return _assignment(dossier, "membrane_pore_filter_oligomeric_ion_selectivity", 4, "assigned", None, leak, failed)
    if xcl1_ops:
        return _assignment(dossier, "metamorphic_two_state_fold_switch", 4, "assigned", None, leak, failed)
    if snca_ops:
        return _assignment(dossier, "intrinsic_disorder_contextual_ensemble", 4, "assigned", None, leak, failed)
    if ("membrane_topology_context" in buckets or "transport_context" in buckets or "transport_or_pump_context" in buckets or "water_channel_context" in buckets) and not ("filter_or_signature_context" in buckets and "ion_selectivity_context" in buckets):
        return _assignment(dossier, "other_membrane_or_transport_context", 3, "assigned", None, leak, failed)
    if "intrinsic_disorder_context" in buckets or "ensemble_context" in buckets:
        return _assignment(dossier, "other_disorder_or_low_confidence_ensemble", 2, "assigned", None, leak, failed)
    if any(bucket in buckets for bucket in ["soluble_chemokine_context", "soluble_single_fold_context", "compact_fold_context", "soluble_enzyme_context", "soluble_heme_protein_context"]):
        return _assignment(dossier, "soluble_single_or_contextual_fold_not_metamorphic", 3, "assigned", None, leak, failed)
    return _assignment(dossier, "insufficient_evidence_clean_abstain", 0, "clean_abstain", "generic_or_insufficient_evidence", leak, failed)


def _assignment(dossier: dict[str, Any], mechanism: str, strength: int, status: str, abstain_reason: str | None, leak: dict[str, Any], failed: list[str]) -> dict[str, Any]:
    return {
        "masked_id": dossier.get("masked_id"),
        "assigned_mechanism_class": mechanism,
        "confidence_strength": strength,
        "assignment_status": status,
        "abstain_reason": abstain_reason,
        "coordinate_derived_source_count": leak["coordinate_derived_source_count"],
        "internal_runtime_source_count": leak["internal_runtime_source_count"],
        "placeholder_source_count": leak["placeholder_source_count"],
        "answer_key_or_target_name_leakage_count": leak["answer_key_or_target_name_leakage_count"],
        "failed_checks": sorted(set(failed)),
    }


def load_masked_panel(panel_root: Path = PANEL_ROOT) -> list[dict[str, Any]]:
    masked_dir = panel_root / "masked_inputs"
    return [_read_json(path, f"masked dossier {path.name}") for path in sorted(masked_dir.glob("TARGET_*.json"))]


def load_unmasked_panel(panel_root: Path = PANEL_ROOT) -> list[dict[str, Any]]:
    dossiers_dir = panel_root / "dossiers"
    return [_read_json(path, f"unmasked dossier {path.name}") for path in sorted(dossiers_dir.glob("TARGET_*.json"))]


def _unmasked_to_assignment_input(dossier: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V38_UNMASKED_ASSIGNMENT_VIEW_v0",
        "masked_id": dossier.get("panel_id"),
        "evidence_buckets": dossier.get("evidence_buckets", []),
        "mechanism_evidence": dossier.get("mechanism_evidence", []),
        "grammar_rules": dossier.get("grammar_rules", {}),
        "source_provenance": dossier.get("source_rows", []),
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
    }


def evaluate_panel(masked_dossiers: list[dict[str, Any]], answer_key: dict[str, Any], unmasked_dossiers: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    assignments = [assign_masked_dossier(row, require_masked=True) for row in masked_dossiers]
    by_id = {row["masked_id"]: row for row in assignments}
    known_ids = [mid for mid, meta in answer_key.items() if meta.get("known_positive") is True]
    decoy_ids = [mid for mid, meta in answer_key.items() if meta.get("known_positive") is not True]
    hard_positive_correct = [
        mid for mid in known_ids
        if by_id.get(mid, {}).get("assigned_mechanism_class") == answer_key[mid]["expected_mechanism_class"]
    ]
    false_neg = [
        mid for mid in known_ids
        if by_id.get(mid, {}).get("assigned_mechanism_class") != answer_key[mid]["expected_mechanism_class"]
    ]
    false_pos = [
        mid for mid in decoy_ids
        if by_id.get(mid, {}).get("assigned_mechanism_class") in HARD_CLASSES
    ]
    confusion: dict[str, dict[str, int]] = {}
    for mid, meta in answer_key.items():
        expected = meta["expected_mechanism_class"]
        observed = by_id.get(mid, {}).get("assigned_mechanism_class", "missing_assignment")
        confusion.setdefault(expected, {})
        confusion[expected][observed] = confusion[expected].get(observed, 0) + 1

    unmasked_consistency = True
    if unmasked_dossiers is not None:
        unmasked_assignments = {
            row["panel_id"]: assign_masked_dossier(_unmasked_to_assignment_input(row), require_masked=False)
            for row in unmasked_dossiers
        }
        for mid, assignment in by_id.items():
            if unmasked_assignments.get(mid, {}).get("assigned_mechanism_class") != assignment.get("assigned_mechanism_class"):
                unmasked_consistency = False
                break

    return {
        "assignments": assignments,
        "known_positive_count": len(known_ids),
        "decoy_count": len(decoy_ids),
        "assigned_count": sum(1 for row in assignments if row["assignment_status"] == "assigned"),
        "clean_abstain_count": sum(1 for row in assignments if row["assignment_status"] == "clean_abstain"),
        "correct_known_positive_count": len(hard_positive_correct),
        "false_positive_hard_grammar_count": len(false_pos),
        "false_negative_known_positive_count": len(false_neg),
        "false_positive_ids": false_pos,
        "false_negative_ids": false_neg,
        "mechanism_confusion_matrix": confusion,
        "unmasked_masked_consistency_passed": unmasked_consistency,
        "coordinate_derived_source_count": sum(row["coordinate_derived_source_count"] for row in assignments),
        "internal_runtime_source_count": sum(row["internal_runtime_source_count"] for row in assignments),
        "placeholder_source_count": sum(row["placeholder_source_count"] for row in assignments),
        "answer_key_or_target_name_leakage_count": sum(row["answer_key_or_target_name_leakage_count"] for row in assignments),
        "failed_checks": sorted(set(check for row in assignments for check in row["failed_checks"])),
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def run_controls(masked: list[dict[str, Any]], answer_key: dict[str, Any], unmasked: list[dict[str, Any]]) -> list[dict[str, Any]]:
    baseline = evaluate_panel(masked, answer_key, unmasked)
    controls = [
        _control("baseline_masked_panel_assignment", baseline["correct_known_positive_count"] == 3 and baseline["false_positive_hard_grammar_count"] == 0, "Known positives correct and decoys not promoted.", baseline),
        _control("unmasked_vs_masked_consistency_check", baseline["unmasked_masked_consistency_passed"], "Masked and unmasked assignment classes must match."),
    ]

    swapped = [dict(row) for row in masked]
    swapped[0], swapped[2] = dict(swapped[2]), dict(swapped[0])
    swapped[0]["masked_id"], swapped[2]["masked_id"] = "TARGET_001", "TARGET_003"
    swapped_eval = evaluate_panel(swapped, answer_key)
    controls.append(_control("swapped_known_positive_dossiers_detected", swapped_eval["false_negative_known_positive_count"] >= 2, "Swapping KcsA/SNCA evidence must fail known-positive assignment.", swapped_eval))

    target4 = next(row for row in masked if row["masked_id"] == "TARGET_004")
    controls.append(_control("membrane_decoy_not_kcsa_without_filter_ion", assign_masked_dossier(target4)["assigned_mechanism_class"] != "membrane_pore_filter_oligomeric_ion_selectivity", "Membrane decoy lacks potassium filter/ion evidence."))
    target6 = next(row for row in masked if row["masked_id"] == "TARGET_006")
    controls.append(_control("chemokine_decoy_not_xcl1_without_two_state", assign_masked_dossier(target6)["assigned_mechanism_class"] != "metamorphic_two_state_fold_switch", "Chemokine decoy lacks two-state evidence."))
    target8 = next(row for row in masked if row["masked_id"] == "TARGET_008")
    controls.append(_control("folded_decoy_not_snca_without_disorder", assign_masked_dossier(target8)["assigned_mechanism_class"] != "intrinsic_disorder_contextual_ensemble", "Folded decoy lacks disorder evidence."))

    generic = {"kind": "V38_MASKED_DOSSIER_v0", "masked_id": "TARGET_999", "evidence_buckets": ["generic_annotation"], "mechanism_evidence": ["generic protein annotation"], "grammar_rules": {}, "source_provenance": []}
    controls.append(_control("generic_only_annotations_clean_abstain", assign_masked_dossier(generic)["assigned_mechanism_class"] == "insufficient_evidence_clean_abstain", "Generic-only annotations must abstain."))
    coordinate = dict(masked[0]); coordinate["source_provenance"] = [{"source_name": "1BL8 coordinates", "coordinate_derived": True}]
    controls.append(_control("coordinate_derived_source_blocked", assign_masked_dossier(coordinate)["assignment_status"] == "blocked_leakage", "Coordinate-derived source must block."))
    internal = dict(masked[0]); internal["source_provenance"] = [{"source_url_or_citation": "first_contact_clean_pharmacotopology_layer_run/V37/report.json"}]
    controls.append(_control("internal_runtime_source_blocked", assign_masked_dossier(internal)["assignment_status"] == "blocked_leakage", "Runtime source must block."))
    placeholder = dict(masked[0]); placeholder["source_provenance"] = [{"source_name": "TODO placeholder", "source_url_or_citation": "citation needed"}]
    controls.append(_control("placeholder_source_blocked", assign_masked_dossier(placeholder)["assignment_status"] == "blocked_key", "Placeholder source must block."))
    name_only = {"kind": "V38_MASKED_DOSSIER_v0", "masked_id": "TARGET_010", "evidence_buckets": [], "mechanism_evidence": ["KcsA"], "grammar_rules": {}, "source_provenance": []}
    controls.append(_control("target_name_only_assignment_attempt_blocked", assign_masked_dossier(name_only)["assignment_status"] == "blocked_key", "Target-name-only signal must block."))
    key_leak = dict(masked[0]); key_leak["expected_mechanism_class"] = "membrane_pore_filter_oligomeric_ion_selectivity"
    controls.append(_control("answer_key_leakage_attempt_blocked", assign_masked_dossier(key_leak)["assignment_status"] == "blocked_key", "Answer-key leakage must block."))
    return controls


def build_v38(panel_root: Path = PANEL_ROOT) -> dict[str, Any]:
    masked = load_masked_panel(panel_root)
    unmasked = load_unmasked_panel(panel_root)
    answer_key = _read_json(panel_root / "answer_key.json", "V38 answer key")
    panel_eval = evaluate_panel(masked, answer_key, unmasked)
    controls = run_controls(masked, answer_key, unmasked)
    failed_checks = list(panel_eval["failed_checks"])
    if not all(row["passed"] for row in controls):
        failed_checks.append("V38_required_controls_failed")

    if panel_eval["answer_key_or_target_name_leakage_count"]:
        status = BLOCKED_KEY
    elif panel_eval["coordinate_derived_source_count"] or panel_eval["internal_runtime_source_count"]:
        status = BLOCKED_LEAKAGE
    elif panel_eval["false_positive_hard_grammar_count"]:
        status = FAILED_FP
    elif panel_eval["correct_known_positive_count"] < panel_eval["known_positive_count"] or not panel_eval["unmasked_masked_consistency_passed"]:
        status = PARTIAL
    elif not all(row["passed"] for row in controls):
        status = PARTIAL
    else:
        status = PASSED

    if status == PASSED:
        next_action = "use_blind_generalized_mechanism_grammar_for_class_specific_readout_design_no_claim_no_MD"
    elif status == FAILED_FP:
        next_action = "tighten_decoy_rejection_before_any_mechanism_generalization_claim"
    elif status == BLOCKED_KEY:
        next_action = "remove_target_name_or_answer_key_leakage_from_blind_assignment"
    elif status == BLOCKED_LEAKAGE:
        next_action = "remove_coordinate_or_internal_runtime_leakage_from_V38"
    else:
        next_action = "complete_or_repair_blind_panel_assignment_before_generalization"

    return {
        "kind": "V38_BLIND_MECHANISM_GENERALIZATION_PANEL_v0",
        "run_mode": "blind_masked_mechanism_generalization_panel_no_coordinates_no_MD_claim_disabled",
        "control_status": status,
        "panel_target_count": len(masked),
        "masked_panel_used": True,
        "answer_key_used_for_assignment": False,
        "known_positive_count": panel_eval["known_positive_count"],
        "decoy_count": panel_eval["decoy_count"],
        "assigned_count": panel_eval["assigned_count"],
        "clean_abstain_count": panel_eval["clean_abstain_count"],
        "correct_known_positive_count": panel_eval["correct_known_positive_count"],
        "false_positive_hard_grammar_count": panel_eval["false_positive_hard_grammar_count"],
        "false_negative_known_positive_count": panel_eval["false_negative_known_positive_count"],
        "mechanism_confusion_matrix": panel_eval["mechanism_confusion_matrix"],
        "target_name_masking_passed": panel_eval["answer_key_or_target_name_leakage_count"] == 0,
        "unmasked_masked_consistency_passed": panel_eval["unmasked_masked_consistency_passed"],
        "swapped_dossier_detection_passed": next(row for row in controls if row["control_id"] == "swapped_known_positive_dossiers_detected")["passed"],
        "coordinate_derived_source_count": panel_eval["coordinate_derived_source_count"],
        "internal_runtime_source_count": panel_eval["internal_runtime_source_count"],
        "placeholder_source_count": panel_eval["placeholder_source_count"],
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
        "control_count": len(controls),
        "passed_control_count": sum(1 for row in controls if row["passed"]),
        "controls": controls,
        "masked_assignment_results": panel_eval["assignments"],
        "answer_key_evaluation": [
            {
                "masked_id": mid,
                "target": meta["target"],
                "group": meta["group"],
                "expected_mechanism_class": meta["expected_mechanism_class"],
                "assigned_mechanism_class": next(row for row in panel_eval["assignments"] if row["masked_id"] == mid)["assigned_mechanism_class"],
            }
            for mid, meta in sorted(answer_key.items())
        ],
        "failed_checks": sorted(set(failed_checks)),
        "next_action": next_action,
        "locked_interpretation": (
            "V38 passing means the mechanism grammar survived masked target names and near-decoy pressure: "
            "the three known positives were recovered, and membrane/soluble/folded decoys were not promoted into the hard positive grammars. "
            "This is mechanism-class discrimination, not atom prediction or a folding-solved claim."
        ),
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V38 Blind Mechanism Generalization Panel",
        "",
        f"Status: `{cert.get('control_status')}`",
        f"Known positive accuracy: `{cert.get('correct_known_positive_count')}` / `{cert.get('known_positive_count')}`",
        f"False positive hard-grammar count: `{cert.get('false_positive_hard_grammar_count')}`",
        f"False negative known-positive count: `{cert.get('false_negative_known_positive_count')}`",
        f"Controls: `{cert.get('passed_control_count')}` / `{cert.get('control_count')}`",
        f"Answer key used for assignment: `{cert.get('answer_key_used_for_assignment')}`",
        f"Target names masked: `{cert.get('masked_panel_used')}`",
        "",
        "## Masked Target Table",
        "| Masked ID | Assigned class | Strength | Abstain reason |",
        "|---|---|---:|---|",
    ]
    for row in cert.get("masked_assignment_results", []):
        lines.append(f"| {row.get('masked_id')} | {row.get('assigned_mechanism_class')} | {row.get('confidence_strength')} | {row.get('abstain_reason') or ''} |")
    lines.extend(["", "## Answer-Key Evaluation", "| Masked ID | Target | Expected | Assigned |", "|---|---|---|---|"])
    for row in cert.get("answer_key_evaluation", []):
        lines.append(f"| {row.get('masked_id')} | {row.get('target')} | {row.get('expected_mechanism_class')} | {row.get('assigned_mechanism_class')} |")
    lines.extend(["", "## Controls"])
    for row in cert.get("controls", []):
        lines.append(f"- `{row.get('control_id')}`: `{'PASS' if row.get('passed') else 'FAIL'}`")
    lines.extend(["", "## Plain English Interpretation", str(cert.get("locked_interpretation", ""))])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, panel_root: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    results_dir = panel_root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v38_blind_mechanism_generalization_panel_certificate.json"
    report_path = out_dir / "V38_BLIND_MECHANISM_GENERALIZATION_PANEL_REPORT.md"
    decision_path = out_dir / "v38_blind_mechanism_generalization_panel_next_decision.json"
    results_path = results_dir / "v38_masked_assignment_results.json"
    evaluation_path = results_dir / "v38_answer_key_evaluation.json"
    decision = {
        "kind": "V38_BLIND_MECHANISM_GENERALIZATION_PANEL_NEXT_DECISION_v0",
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
            "masked_assignment_results": str(results_path),
            "answer_key_evaluation": str(evaluation_path),
        },
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, decision)
    _write_json(results_path, cert["masked_assignment_results"])
    _write_json(evaluation_path, cert["answer_key_evaluation"])
    _write_report(report_path, cert)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path, "masked_assignment_results": results_path, "answer_key_evaluation": evaluation_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V38 blind mechanism generalization panel.")
    parser.add_argument("--panel-root", type=Path, default=PANEL_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    cert = build_v38(args.panel_root)
    paths = write_outputs(args.out_dir, args.panel_root, cert)
    print(json.dumps({
        "kind": cert.get("kind"),
        "control_status": cert.get("control_status"),
        "known_positive_accuracy": f"{cert.get('correct_known_positive_count')}/{cert.get('known_positive_count')}",
        "false_positive_hard_grammar_count": cert.get("false_positive_hard_grammar_count"),
        "false_negative_known_positive_count": cert.get("false_negative_known_positive_count"),
        "control_count": cert.get("control_count"),
        "passed_control_count": cert.get("passed_control_count"),
        "coordinate_derived_source_count": cert.get("coordinate_derived_source_count"),
        "internal_runtime_source_count": cert.get("internal_runtime_source_count"),
        "masked_panel_used": cert.get("masked_panel_used"),
        "answer_key_used_for_assignment": cert.get("answer_key_used_for_assignment"),
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
