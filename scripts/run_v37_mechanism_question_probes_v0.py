#!/usr/bin/env python3
from __future__ import annotations

"""V37 mechanism question probes.

V37 asks whether the V36 evidence dossiers can assign the correct mechanism
grammar class without native coordinates, MD, or target-name-only shortcuts.
"""

import argparse
import copy
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "external_evidence_dossiers"
MAP_ROOT = REPO_ROOT / "data" / "mechanism_maps" / "V37"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V37_MECHANISM_QUESTION_PROBES"

INPUT_V36_STATUS = "V36_REAL_EVIDENCE_DOSSIERS_READY_CLAIM_DISABLED"

PASSED = "V37_MECHANISM_QUESTION_PROBES_PASSED_CLAIM_DISABLED"
PARTIAL = "V37_PARTIAL_MECHANISM_ASSIGNMENT_CLEAN_ABSTAIN"
BLOCKED_LEAKAGE = "V37_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"
BLOCKED_NAME_ONLY = "V37_BLOCKED_TARGET_NAME_ONLY_ASSIGNMENT"
BLOCKED_WRONG_GRAMMAR = "V37_BLOCKED_FORCED_WRONG_GRAMMAR"

TARGET_ORDER = ["KcsA", "XCL1_lymphotactin", "alpha_synuclein_SNCA"]

TARGET_SPECS: dict[str, dict[str, Any]] = {
    "KcsA": {
        "mechanism_class": "membrane_pore_filter_oligomeric_ion_selectivity",
        "required": [
            "membrane_topology_operator",
            "ion_selectivity_operator",
            "filter_signature_operator",
            "oligomer_or_interface_context_operator",
        ],
        "forbidden": [
            "soluble_single_compact_core_operator",
            "metamorphic_two_state_operator",
            "intrinsic_disorder_ensemble_operator",
            "de_novo_whole_fold_claim_operator",
        ],
        "questions": [
            "Does the non-coordinate dossier support potassium-channel identity?",
            "Does it support membrane/pore/filter context?",
            "Does it separate filter/ion-selectivity grammar from whole-channel folding claims?",
            "Does it avoid using coordinate contacts?",
        ],
    },
    "XCL1_lymphotactin": {
        "mechanism_class": "metamorphic_two_state_fold_switch",
        "required": [
            "state_A_chemokine_monomer_operator",
            "state_B_beta_sandwich_dimer_operator",
            "state_specific_function_operator",
            "no_mixed_state_pooling_operator",
        ],
        "forbidden": [
            "single_consensus_fold_operator",
            "membrane_pore_filter_operator",
            "intrinsic_disorder_single_ensemble_operator",
            "de_novo_whole_fold_claim_operator",
        ],
        "questions": [
            "Does the dossier support two biologically relevant states?",
            "Does it prevent mixed-state pooling?",
            "Does it keep state A and state B function evidence separate?",
            "Does it refuse a single-fold claim?",
        ],
    },
    "alpha_synuclein_SNCA": {
        "mechanism_class": "intrinsic_disorder_contextual_ensemble",
        "required": [
            "intrinsic_disorder_operator",
            "ensemble_not_single_fold_operator",
            "disorder_to_order_context_operator",
            "context_bound_structure_operator",
        ],
        "forbidden": [
            "compact_native_single_fold_operator",
            "membrane_channel_pore_operator",
            "metamorphic_two_native_folds_operator",
            "de_novo_whole_fold_claim_operator",
        ],
        "questions": [
            "Does the dossier support free-state disorder?",
            "Does it represent membrane-bound helix as context-dependent, not solved folding?",
            "Does it block single-native-fold grammar?",
            "Does it keep aggregation/amyloid context from becoming a native-fold claim?",
        ],
    },
}

MECHANISM_TO_TARGET = {spec["mechanism_class"]: target for target, spec in TARGET_SPECS.items()}

PLACEHOLDER_TOKENS = {"placeholder", "todo", "tbd", "citation needed", "example.com", "dummy", "unknown"}
COORDINATE_TOKENS = {
    "1bl8",
    "rcsb",
    "pdb coordinate",
    "pdb_coordinate",
    "coordinate-derived",
    "coordinate_derived",
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


def _text(dossier: dict[str, Any]) -> str:
    pieces = [
        dossier.get("problem_class"),
        dossier.get("grammar_focus"),
        dossier.get("locked_interpretation"),
        dossier.get("display_name"),
        dossier.get("target"),
        json.dumps(dossier.get("grammar_rules", {}), sort_keys=True),
    ]
    return " ".join(str(piece or "") for piece in pieces).lower()


def _buckets(dossier: dict[str, Any]) -> set[str]:
    status = dossier.get("bucket_status") if isinstance(dossier.get("bucket_status"), dict) else {}
    present = _as_list(status.get("present_buckets"))
    if not present:
        present = _as_list(dossier.get("required_buckets"))
    return {str(bucket) for bucket in present}


def _rules(dossier: dict[str, Any]) -> dict[str, Any]:
    rules = dossier.get("grammar_rules")
    return rules if isinstance(rules, dict) else {}


def _source_rows(dossier: dict[str, Any]) -> list[dict[str, Any]]:
    return [row for row in _as_list(dossier.get("source_rows")) if isinstance(row, dict)]


def _row_text(row: dict[str, Any]) -> str:
    return " ".join(str(value or "") for value in row.values()).lower()


def _contains_placeholder(value: Any) -> bool:
    text = str(value or "").lower().strip()
    if not text:
        return True
    return any(token in text for token in PLACEHOLDER_TOKENS)


def _source_leakage_counts(dossier: dict[str, Any]) -> dict[str, Any]:
    coordinate_count = 0
    internal_count = 0
    placeholder_count = 0
    native_metric_count = 0
    failed: list[str] = []

    if _bool_true(dossier.get("coordinate_truth_used_before_selection")):
        coordinate_count += 1
        failed.append("dossier_coordinate_truth_used_before_selection")
    if _bool_true(dossier.get("native_metrics_used_for_selection")):
        native_metric_count += 1
        failed.append("dossier_native_metrics_used_for_selection")
    if _bool_true(dossier.get("claim_allowed")):
        failed.append("dossier_claim_allowed_true")
    if _bool_true(dossier.get("positive_folding_evidence_found")) or _bool_true(dossier.get("folding_problem_solved")):
        failed.append("dossier_folding_claim_true")

    for index, row in enumerate(_source_rows(dossier)):
        row_text = _row_text(row)
        if _bool_true(row.get("coordinate_derived")) or _bool_true(row.get("coordinate_truth_used_before_selection")) or any(token in row_text for token in COORDINATE_TOKENS):
            coordinate_count += 1
            failed.append(f"source_row_{index}:coordinate_or_prediction_leakage")
        if _bool_true(row.get("internal_runtime_source")) or "first_contact_clean_pharmacotopology_layer_run" in row_text:
            internal_count += 1
            failed.append(f"source_row_{index}:internal_runtime_source")
        if _bool_true(row.get("native_metrics_used_for_selection")):
            native_metric_count += 1
            failed.append(f"source_row_{index}:native_metrics_used_for_selection")
        for key in ["source_name", "source_url_or_citation", "source_date_or_version"]:
            if key in row and _contains_placeholder(row.get(key)):
                placeholder_count += 1
                failed.append(f"source_row_{index}:placeholder_or_missing_{key}")

    return {
        "coordinate_derived_source_count": coordinate_count,
        "internal_runtime_source_count": internal_count,
        "placeholder_source_count": placeholder_count,
        "native_metric_source_count": native_metric_count,
        "failed_checks": sorted(set(failed)),
    }


def detect_required_operators(dossier: dict[str, Any], mechanism_class: str) -> list[str]:
    buckets = _buckets(dossier)
    rules = _rules(dossier)
    text = _text(dossier)
    found: list[str] = []

    if mechanism_class == "membrane_pore_filter_oligomeric_ion_selectivity":
        if "membrane_topology_context" in buckets:
            found.append("membrane_topology_operator")
        if "ion_selectivity_context" in buckets:
            found.append("ion_selectivity_operator")
        if "filter_or_signature_context" in buckets:
            found.append("filter_signature_operator")
        if "interface" in text or "oligomer" in text or "homotetramer" in text:
            found.append("oligomer_or_interface_context_operator")

    elif mechanism_class == "metamorphic_two_state_fold_switch":
        if "state_A_function_context" in buckets:
            found.append("state_A_chemokine_monomer_operator")
        if "state_B_function_context" in buckets:
            found.append("state_B_beta_sandwich_dimer_operator")
        if {"state_A_function_context", "state_B_function_context"} <= buckets:
            found.append("state_specific_function_operator")
        if "no_mixed_state_pooling_rule" in buckets and rules.get("mixed_state_pooling_allowed") is False:
            found.append("no_mixed_state_pooling_operator")

    elif mechanism_class == "intrinsic_disorder_contextual_ensemble":
        if "intrinsic_disorder_context" in buckets:
            found.append("intrinsic_disorder_operator")
        if "no_single_native_fold_rule" in buckets or rules.get("ensemble_not_single_fold") is True:
            found.append("ensemble_not_single_fold_operator")
        if "disorder_to_order_context" in buckets:
            found.append("disorder_to_order_context_operator")
        if "disorder_to_order_context" in buckets:
            found.append("context_bound_structure_operator")

    return sorted(set(found))


def detect_forbidden_operators(dossier: dict[str, Any], expected_target: str) -> list[str]:
    buckets = _buckets(dossier)
    rules = _rules(dossier)
    text = _text(dossier)
    forbidden: list[str] = []

    if _bool_true(dossier.get("claim_allowed")) or _bool_true(dossier.get("positive_folding_evidence_found")) or _bool_true(dossier.get("folding_problem_solved")):
        forbidden.append("de_novo_whole_fold_claim_operator")

    has_kcsa = bool({"ion_selectivity_context", "membrane_topology_context", "filter_or_signature_context"} & buckets) or "pore/filter" in text or "potassium" in text
    has_xcl1 = bool({"state_A_function_context", "state_B_function_context", "metamorphic_two_state_context"} & buckets) or "metamorphic" in text or "state-switch" in text
    has_snca = bool({"intrinsic_disorder_context", "ensemble_context", "disorder_to_order_context"} & buckets) or "intrinsic" in text and "disorder" in text

    if expected_target == "KcsA":
        if has_xcl1:
            forbidden.append("metamorphic_two_state_operator")
        if has_snca:
            forbidden.append("intrinsic_disorder_ensemble_operator")
        if "single compact" in text:
            forbidden.append("soluble_single_compact_core_operator")
    elif expected_target == "XCL1_lymphotactin":
        if has_kcsa:
            forbidden.append("membrane_pore_filter_operator")
        if has_snca:
            forbidden.append("intrinsic_disorder_single_ensemble_operator")
        if rules.get("mixed_state_pooling_allowed") is True or "single consensus" in text:
            forbidden.append("single_consensus_fold_operator")
    elif expected_target == "alpha_synuclein_SNCA":
        if has_kcsa:
            forbidden.append("membrane_channel_pore_operator")
        if has_xcl1:
            forbidden.append("metamorphic_two_native_folds_operator")
        if rules.get("single_native_fold_model_allowed") is True or "compact native single fold" in text:
            forbidden.append("compact_native_single_fold_operator")

    return sorted(set(forbidden))


def assign_mechanism_class(dossier: dict[str, Any]) -> dict[str, Any]:
    scores: dict[str, int] = {}
    operators_by_class: dict[str, list[str]] = {}
    text = _text(dossier)
    for spec in TARGET_SPECS.values():
        mechanism_class = spec["mechanism_class"]
        ops = detect_required_operators(dossier, mechanism_class)
        operators_by_class[mechanism_class] = ops
        score = len(ops)
        if mechanism_class.startswith("membrane") and ("pore/filter" in text or "potassium-channel" in text or "potassium channel" in text):
            score += 1
        if mechanism_class.startswith("metamorphic") and ("metamorphic" in text or "state-switch" in text):
            score += 1
        if mechanism_class.startswith("intrinsic") and ("idp" in text or "ensemble" in text or "disorder" in text):
            score += 1
        scores[mechanism_class] = score

    ordered = sorted(scores.items(), key=lambda item: item[1], reverse=True)
    if not ordered or ordered[0][1] < 3:
        assigned = None
    elif len(ordered) > 1 and ordered[0][1] == ordered[1][1]:
        assigned = None
    else:
        assigned = ordered[0][0]
    return {
        "assigned_mechanism_class": assigned,
        "assignment_scores": scores,
        "operators_by_class": operators_by_class,
    }


def answer_scientific_questions(expected_target: str, dossier: dict[str, Any], required_found: list[str], forbidden: list[str], leakage: dict[str, Any]) -> list[dict[str, Any]]:
    found = set(required_found)
    no_leakage = leakage["coordinate_derived_source_count"] == 0 and leakage["internal_runtime_source_count"] == 0
    rules = _rules(dossier)
    rows: list[dict[str, Any]] = []

    if expected_target == "KcsA":
        answers = [
            "ion_selectivity_operator" in found,
            {"membrane_topology_operator", "filter_signature_operator"} <= found,
            "de_novo_whole_fold_claim_operator" not in forbidden and _bool_false(dossier.get("folding_problem_solved")),
            no_leakage,
        ]
    elif expected_target == "XCL1_lymphotactin":
        answers = [
            {"state_A_chemokine_monomer_operator", "state_B_beta_sandwich_dimer_operator"} <= found,
            "no_mixed_state_pooling_operator" in found and rules.get("mixed_state_pooling_allowed") is False,
            "state_specific_function_operator" in found,
            "single_consensus_fold_operator" not in forbidden and _bool_false(dossier.get("folding_problem_solved")),
        ]
    else:
        answers = [
            "intrinsic_disorder_operator" in found,
            "context_bound_structure_operator" in found and _bool_false(dossier.get("folding_problem_solved")),
            "ensemble_not_single_fold_operator" in found and "compact_native_single_fold_operator" not in forbidden,
            "de_novo_whole_fold_claim_operator" not in forbidden and _bool_false(dossier.get("positive_folding_evidence_found")),
        ]

    for question, answer in zip(TARGET_SPECS[expected_target]["questions"], answers):
        rows.append({"question": question, "answer": bool(answer)})
    return rows


def evaluate_dossier(expected_target: str, dossier: dict[str, Any]) -> dict[str, Any]:
    expected_class = TARGET_SPECS[expected_target]["mechanism_class"]
    assignment = assign_mechanism_class(dossier)
    assigned = assignment["assigned_mechanism_class"]
    required_found = detect_required_operators(dossier, expected_class)
    required_missing = [op for op in TARGET_SPECS[expected_target]["required"] if op not in required_found]
    forbidden = detect_forbidden_operators(dossier, expected_target)
    leakage = _source_leakage_counts(dossier)
    failed_checks = [f"{expected_target}:{check}" for check in leakage["failed_checks"]]

    if assigned is None:
        failed_checks.append(f"{expected_target}:no_high_confidence_mechanism_assignment")
    elif assigned != expected_class:
        failed_checks.append(f"{expected_target}:mechanism_class_mismatch:{assigned}")
    for op in required_missing:
        failed_checks.append(f"{expected_target}:missing_required_operator:{op}")
    for op in forbidden:
        failed_checks.append(f"{expected_target}:forbidden_operator_triggered:{op}")

    if leakage["coordinate_derived_source_count"] or leakage["internal_runtime_source_count"] or leakage["native_metric_source_count"]:
        target_status = "blocked_leakage"
    elif leakage["placeholder_source_count"] or forbidden or (assigned is not None and assigned != expected_class):
        target_status = "blocked_wrong_grammar"
    elif assigned is None or required_missing:
        target_status = "partial"
    else:
        target_status = "assigned"

    answers = answer_scientific_questions(expected_target, dossier, required_found, forbidden, leakage)
    return {
        "target": expected_target,
        "dossier_target_label": dossier.get("target"),
        "expected_mechanism_class": expected_class,
        "assigned_mechanism_class": assigned,
        "target_status": target_status,
        "assignment_scores": assignment["assignment_scores"],
        "required_operator_buckets_present": required_found,
        "required_operator_buckets_missing": required_missing,
        "forbidden_operator_buckets_triggered": forbidden,
        "scientific_question_answers": answers,
        "coordinate_derived_source_count": leakage["coordinate_derived_source_count"],
        "internal_runtime_source_count": leakage["internal_runtime_source_count"],
        "placeholder_source_count": leakage["placeholder_source_count"],
        "native_metric_source_count": leakage["native_metric_source_count"],
        "failed_checks": sorted(set(failed_checks)),
    }


def load_dossiers(data_root: Path = DATA_ROOT) -> dict[str, dict[str, Any]]:
    return {
        target: _read_json(data_root / target / "evidence_dossier.json", f"{target} V36 evidence dossier")
        for target in TARGET_ORDER
    }


def _mask_dossier_names(packages: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    masked = copy.deepcopy(packages)
    aliases = {"KcsA": "TARGET_A", "XCL1_lymphotactin": "TARGET_B", "alpha_synuclein_SNCA": "TARGET_C"}
    for target, dossier in masked.items():
        alias = aliases[target]
        dossier["target"] = alias
        dossier["display_name"] = alias
        if isinstance(dossier.get("locked_interpretation"), str):
            dossier["locked_interpretation"] = dossier["locked_interpretation"].replace(target, alias)
    return masked


def target_name_masking_passes(packages: dict[str, dict[str, Any]]) -> bool:
    masked = _mask_dossier_names(packages)
    for target in TARGET_ORDER:
        result = evaluate_dossier(target, masked[target])
        if result["assigned_mechanism_class"] != TARGET_SPECS[target]["mechanism_class"]:
            return False
    return True


def swapped_dossier_detection_passes(packages: dict[str, dict[str, Any]]) -> bool:
    kcsa_as_snca = copy.deepcopy(packages["KcsA"])
    kcsa_as_snca["target"] = "alpha_synuclein_SNCA"
    snca_result = evaluate_dossier("alpha_synuclein_SNCA", kcsa_as_snca)
    snca_as_kcsa = copy.deepcopy(packages["alpha_synuclein_SNCA"])
    snca_as_kcsa["target"] = "KcsA"
    kcsa_result = evaluate_dossier("KcsA", snca_as_kcsa)
    return (
        snca_result["target_status"] == "blocked_wrong_grammar"
        and kcsa_result["target_status"] == "blocked_wrong_grammar"
        and snca_result["assigned_mechanism_class"] == TARGET_SPECS["KcsA"]["mechanism_class"]
        and kcsa_result["assigned_mechanism_class"] == TARGET_SPECS["alpha_synuclein_SNCA"]["mechanism_class"]
    )


def _aggregate_status(target_results: list[dict[str, Any]], masking_passed: bool, swapped_passed: bool, controls: list[dict[str, Any]]) -> str:
    if any(row["coordinate_derived_source_count"] or row["internal_runtime_source_count"] or row["native_metric_source_count"] for row in target_results):
        return BLOCKED_LEAKAGE
    if not masking_passed:
        return BLOCKED_NAME_ONLY
    if not swapped_passed:
        return BLOCKED_WRONG_GRAMMAR
    if any(row["target_status"] == "blocked_wrong_grammar" for row in target_results):
        return BLOCKED_WRONG_GRAMMAR
    if any(control.get("passed") is not True for control in controls):
        return BLOCKED_WRONG_GRAMMAR
    if any(row["target_status"] != "assigned" for row in target_results):
        return PARTIAL
    return PASSED


def build_v37(packages: dict[str, dict[str, Any]]) -> dict[str, Any]:
    target_results = [evaluate_dossier(target, packages[target]) for target in TARGET_ORDER]
    masking_passed = target_name_masking_passes(packages)
    swapped_passed = swapped_dossier_detection_passes(packages)
    controls = run_v37_controls(packages)
    status = _aggregate_status(target_results, masking_passed, swapped_passed, controls)

    assigned = [row for row in target_results if row["target_status"] == "assigned"]
    partial = [row for row in target_results if row["target_status"] == "partial"]
    blocked = [row for row in target_results if row["target_status"].startswith("blocked")]
    mechanism_classes_assigned = {row["target"]: row["assigned_mechanism_class"] for row in target_results}
    required_present = {row["target"]: row["required_operator_buckets_present"] for row in target_results}
    forbidden_triggered = {row["target"]: row["forbidden_operator_buckets_triggered"] for row in target_results}

    failed_checks: list[str] = []
    for row in target_results:
        failed_checks.extend(row["failed_checks"])
    if not masking_passed:
        failed_checks.append("target_name_masking_failed_assignment_depends_on_target_name")
    if not swapped_passed:
        failed_checks.append("swapped_dossier_detection_failed")
    if any(control.get("passed") is not True for control in controls):
        failed_checks.append("V37_required_controls_failed")

    if status == PASSED:
        next_action = "use_mechanism_question_maps_to_design_class_specific_readouts_no_claim_no_MD"
    elif status == PARTIAL:
        next_action = "complete_missing_operator_evidence_before_mechanism_map_claims"
    elif status == BLOCKED_LEAKAGE:
        next_action = "remove_coordinate_internal_or_native_metric_leakage_before_V37"
    elif status == BLOCKED_NAME_ONLY:
        next_action = "fix_assignment_logic_so_masked_target_names_still_work"
    else:
        next_action = "fix_forced_wrong_grammar_or_forbidden_operator_triggers_before_V37"

    return {
        "kind": "V37_MECHANISM_QUESTION_PROBES_v0",
        "run_mode": "mechanism_class_assignment_from_v36_noncoordinate_dossiers_no_coordinates_no_MD_claim_disabled",
        "input_v36_status": INPUT_V36_STATUS,
        "control_status": status,
        "target_count": len(TARGET_ORDER),
        "assigned_target_count": len(assigned),
        "partial_target_count": len(partial),
        "blocked_target_count": len(blocked),
        "mechanism_classes_assigned": mechanism_classes_assigned,
        "required_operator_buckets_present": required_present,
        "forbidden_operator_buckets_triggered": forbidden_triggered,
        "target_name_masking_passed": masking_passed,
        "swapped_dossier_detection_passed": swapped_passed,
        "coordinate_derived_source_count": sum(row["coordinate_derived_source_count"] for row in target_results),
        "internal_runtime_source_count": sum(row["internal_runtime_source_count"] for row in target_results),
        "placeholder_source_count": sum(row["placeholder_source_count"] for row in target_results),
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
        "control_count": len(controls),
        "passed_control_count": sum(1 for control in controls if control.get("passed") is True),
        "controls": controls,
        "target_results": target_results,
        "failed_checks": sorted(set(failed_checks)),
        "next_action": next_action,
        "locked_interpretation": (
            "V37 passing means the V36 non-coordinate dossiers can classify the folding-problem type: "
            "KcsA as membrane pore/filter ion-selectivity grammar, XCL1 as metamorphic state-switch grammar, "
            "and SNCA as intrinsic-disorder ensemble grammar. This is not structure prediction, not MD, and not a folding-solved claim."
        ),
    }


def _minimal_dossier(target: str, buckets: list[str], rules: dict[str, Any], focus: str, problem_class: str = "") -> dict[str, Any]:
    return {
        "kind": "V36_TARGET_EVIDENCE_DOSSIER_v0",
        "target": target,
        "display_name": target,
        "problem_class": problem_class,
        "grammar_focus": focus,
        "bucket_status": {"present_buckets": buckets, "required_buckets": buckets, "missing_buckets": [], "complete": True},
        "required_buckets": buckets,
        "grammar_rules": rules,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
    }


def _ready_fixture_dossiers() -> dict[str, dict[str, Any]]:
    return {
        "KcsA": _minimal_dossier(
            "KcsA",
            ["sequence_or_family_identity", "ion_selectivity_context", "membrane_topology_context", "filter_or_signature_context"],
            {"coordinate_contact_tables_allowed": False, "native_metrics_before_selection_allowed": False},
            "membrane channel pore/filter/interface grammar potassium ion selectivity",
            "membrane_pore_filter_oligomeric_interface_problem",
        ),
        "XCL1_lymphotactin": _minimal_dossier(
            "XCL1_lymphotactin",
            ["state_A_function_context", "state_B_function_context", "metamorphic_two_state_context", "no_mixed_state_pooling_rule"],
            {"mixed_state_pooling_allowed": False, "native_metrics_before_selection_allowed": False},
            "metamorphic state-switch grammar with state-specific function separation",
            "metamorphic_two_fold_state_switch_problem",
        ),
        "alpha_synuclein_SNCA": _minimal_dossier(
            "alpha_synuclein_SNCA",
            ["intrinsic_disorder_context", "ensemble_context", "disorder_to_order_context", "no_single_native_fold_rule"],
            {"single_native_fold_model_allowed": False, "ensemble_not_single_fold": True, "native_metrics_before_selection_allowed": False},
            "IDP ensemble disorder-to-order context-bound membrane-bound helix grammar",
            "intrinsically_disordered_ensemble_problem",
        ),
    }


def _control(control_id: str, observed: dict[str, Any], expected_status: str, expected_failure_contains: str, reason: str) -> dict[str, Any]:
    checks = " ".join(str(item) for item in observed.get("failed_checks", []))
    passed = observed.get("control_status") == expected_status and (not expected_failure_contains or expected_failure_contains in checks)
    return {
        "control_id": control_id,
        "expected_status": expected_status,
        "observed_status": observed.get("control_status"),
        "expected_failure_contains": expected_failure_contains,
        "observed_failed_checks": observed.get("failed_checks"),
        "passed": passed,
        "reason": reason,
    }


def _evaluate_control_packages(packages: dict[str, dict[str, Any]]) -> dict[str, Any]:
    target_results = [evaluate_dossier(target, packages[target]) for target in TARGET_ORDER]
    status = _aggregate_status(target_results, True, True, [])
    failed: list[str] = []
    for row in target_results:
        failed.extend(row["failed_checks"])
    return {"control_status": status, "target_results": target_results, "failed_checks": sorted(set(failed))}


def run_v37_controls(packages: dict[str, dict[str, Any]] | None = None) -> list[dict[str, Any]]:
    base = copy.deepcopy(packages) if packages is not None else _ready_fixture_dossiers()
    controls: list[dict[str, Any]] = []

    controls.append(_control(
        "baseline_v36_dossiers_assign_all_three_mechanism_classes",
        _evaluate_control_packages(copy.deepcopy(base)),
        PASSED,
        "",
        "Baseline V36 dossiers must assign KcsA, XCL1, and SNCA mechanism classes.",
    ))

    masked = _mask_dossier_names(copy.deepcopy(base))
    controls.append(_control(
        "target_name_masked_dossiers_still_assign_from_evidence_content",
        _evaluate_control_packages(masked),
        PASSED,
        "",
        "Masking target names must not break evidence-content assignment.",
    ))

    swapped = copy.deepcopy(base)
    swapped["KcsA"], swapped["alpha_synuclein_SNCA"] = copy.deepcopy(base["alpha_synuclein_SNCA"]), copy.deepcopy(base["KcsA"])
    swapped["KcsA"]["target"] = "KcsA"
    swapped["alpha_synuclein_SNCA"]["target"] = "alpha_synuclein_SNCA"
    controls.append(_control(
        "swapped_dossiers_are_detected_not_renamed_into_target_grammar",
        _evaluate_control_packages(swapped),
        BLOCKED_WRONG_GRAMMAR,
        "mechanism_class_mismatch",
        "Swapped evidence must be detected even when target names are changed.",
    ))

    kcsa_partial = copy.deepcopy(base)
    kcsa_partial["KcsA"]["bucket_status"]["present_buckets"] = ["sequence_or_family_identity", "membrane_topology_context"]
    kcsa_partial["KcsA"]["required_buckets"] = ["sequence_or_family_identity", "membrane_topology_context"]
    controls.append(_control(
        "kcsa_without_ion_filter_evidence_is_partial",
        _evaluate_control_packages(kcsa_partial),
        PARTIAL,
        "KcsA:missing_required_operator:ion_selectivity_operator",
        "KcsA must not be ready without ion/filter operators.",
    ))

    xcl1_partial = copy.deepcopy(base)
    xcl1_partial["XCL1_lymphotactin"]["bucket_status"]["present_buckets"] = [
        "state_A_function_context",
        "metamorphic_two_state_context",
        "no_mixed_state_pooling_rule",
    ]
    controls.append(_control(
        "xcl1_without_state_b_is_partial",
        _evaluate_control_packages(xcl1_partial),
        PARTIAL,
        "XCL1_lymphotactin:missing_required_operator:state_B_beta_sandwich_dimer_operator",
        "XCL1 needs state B evidence.",
    ))

    xcl1_mixed = copy.deepcopy(base)
    xcl1_mixed["XCL1_lymphotactin"]["grammar_rules"]["mixed_state_pooling_allowed"] = True
    controls.append(_control(
        "xcl1_forced_mixed_state_pooling_is_blocked",
        _evaluate_control_packages(xcl1_mixed),
        BLOCKED_WRONG_GRAMMAR,
        "forbidden_operator_triggered:single_consensus_fold_operator",
        "Forced mixed-state pooling must be blocked.",
    ))

    snca_partial = copy.deepcopy(base)
    snca_partial["alpha_synuclein_SNCA"]["bucket_status"]["present_buckets"] = [
        "ensemble_context",
        "disorder_to_order_context",
        "no_single_native_fold_rule",
    ]
    controls.append(_control(
        "snca_without_disorder_evidence_is_partial",
        _evaluate_control_packages(snca_partial),
        PARTIAL,
        "alpha_synuclein_SNCA:missing_required_operator:intrinsic_disorder_operator",
        "SNCA cannot be ready without disorder evidence.",
    ))

    snca_single = copy.deepcopy(base)
    snca_single["alpha_synuclein_SNCA"]["grammar_rules"]["single_native_fold_model_allowed"] = True
    controls.append(_control(
        "snca_forced_compact_single_fold_is_blocked",
        _evaluate_control_packages(snca_single),
        BLOCKED_WRONG_GRAMMAR,
        "forbidden_operator_triggered:compact_native_single_fold_operator",
        "SNCA compact single-fold grammar is forbidden.",
    ))

    generic = copy.deepcopy(base)
    generic["KcsA"] = _minimal_dossier("KcsA", ["sequence_or_family_identity"], {"native_metrics_before_selection_allowed": False}, "generic protein annotation")
    controls.append(_control(
        "generic_only_annotations_do_not_assign_high_confidence_grammar",
        _evaluate_control_packages(generic),
        PARTIAL,
        "KcsA:no_high_confidence_mechanism_assignment",
        "Generic annotations must not become high-confidence mechanism grammar.",
    ))

    coordinate = copy.deepcopy(base)
    coordinate["KcsA"]["source_rows"] = [{
        "source_name": "KcsA 1BL8 coordinate contact CSV",
        "source_url_or_citation": "data/external_constraints/KcsA/pore_filter/kcsa_1bl8_pore_filter_external_contacts.csv",
        "coordinate_derived": True,
    }]
    controls.append(_control(
        "coordinate_derived_source_supplied_to_v37_is_blocked",
        _evaluate_control_packages(coordinate),
        BLOCKED_LEAKAGE,
        "coordinate_or_prediction_leakage",
        "Coordinate-derived V33/V34-style evidence must be blocked.",
    ))

    internal = copy.deepcopy(base)
    internal["KcsA"]["source_rows"] = [{
        "source_name": "internal runtime report",
        "source_url_or_citation": "first_contact_clean_pharmacotopology_layer_run/V36/report.json",
        "internal_runtime_source": True,
    }]
    controls.append(_control(
        "internal_runtime_source_supplied_to_v37_is_blocked",
        _evaluate_control_packages(internal),
        BLOCKED_LEAKAGE,
        "internal_runtime_source",
        "Runtime reports must not become mechanism evidence.",
    ))

    placeholder = copy.deepcopy(base)
    placeholder["KcsA"]["source_rows"] = [{
        "source_name": "TODO placeholder",
        "source_url_or_citation": "citation needed",
        "source_date_or_version": "unknown",
    }]
    controls.append(_control(
        "placeholder_source_row_supplied_to_v37_is_blocked",
        _evaluate_control_packages(placeholder),
        BLOCKED_WRONG_GRAMMAR,
        "placeholder_or_missing_source_name",
        "Placeholder source rows must be blocked.",
    ))

    return controls


def _target_map(result: dict[str, Any]) -> dict[str, Any]:
    return {
        "kind": "V37_TARGET_MECHANISM_QUESTION_MAP_v0",
        "target": result["target"],
        "mechanism_class": result["assigned_mechanism_class"],
        "expected_mechanism_class": result["expected_mechanism_class"],
        "target_status": result["target_status"],
        "required_operator_buckets_present": result["required_operator_buckets_present"],
        "required_operator_buckets_missing": result["required_operator_buckets_missing"],
        "forbidden_operator_buckets_triggered": result["forbidden_operator_buckets_triggered"],
        "scientific_question_answers": result["scientific_question_answers"],
        "claim_allowed": False,
        "new_MD_allowed": False,
        "folding_problem_solved": False,
    }


def _write_answer_table(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V37 Mechanism Question Answer Table",
        "",
        "| Target | Mechanism class | Required operators | Forbidden operators | Scientific questions all true |",
        "|---|---|---|---|---|",
    ]
    for row in cert.get("target_results", []):
        all_true = all(answer.get("answer") is True for answer in row.get("scientific_question_answers", []))
        lines.append(
            "| {target} | {mechanism} | {required} | {forbidden} | {all_true} |".format(
                target=row.get("target"),
                mechanism=row.get("assigned_mechanism_class"),
                required=", ".join(row.get("required_operator_buckets_present", [])),
                forbidden=", ".join(row.get("forbidden_operator_buckets_triggered", [])) or "none",
                all_true=all_true,
            )
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V37 Mechanism Question Probes",
        "",
        f"Status: `{cert.get('control_status')}`",
        f"Input V36 status: `{cert.get('input_v36_status')}`",
        f"Assigned targets: `{cert.get('assigned_target_count')}` / `{cert.get('target_count')}`",
        f"Partial targets: `{cert.get('partial_target_count')}`",
        f"Blocked targets: `{cert.get('blocked_target_count')}`",
        f"Target-name masking passed: `{cert.get('target_name_masking_passed')}`",
        f"Swapped dossier detection passed: `{cert.get('swapped_dossier_detection_passed')}`",
        f"Coordinate-derived source count: `{cert.get('coordinate_derived_source_count')}`",
        f"Internal runtime source count: `{cert.get('internal_runtime_source_count')}`",
        f"Controls: `{cert.get('passed_control_count')}` / `{cert.get('control_count')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD allowed: `{cert.get('new_MD_allowed')}`",
        f"Folding solved: `{cert.get('folding_problem_solved')}`",
        f"Next action: `{cert.get('next_action')}`",
        "",
        "## Targets",
    ]
    for row in cert.get("target_results", []):
        lines.extend([
            f"### {row.get('target')}",
            f"Mechanism class: `{row.get('assigned_mechanism_class')}`",
            f"Required operators found: `{row.get('required_operator_buckets_present')}`",
            f"Forbidden operators triggered: `{row.get('forbidden_operator_buckets_triggered')}`",
            "",
            "Scientific questions:",
        ])
        for answer in row.get("scientific_question_answers", []):
            lines.append(f"- `{answer.get('answer')}` {answer.get('question')}")
        lines.append("")
    lines.extend(["## Controls"])
    for control in cert.get("controls", []):
        lines.append(f"- `{control.get('control_id')}`: `{'PASS' if control.get('passed') else 'FAIL'}`")
    lines.extend(["", "## Failed Checks"])
    failures = cert.get("failed_checks") or []
    lines.extend([f"- `{failure}`" for failure in failures] if failures else ["- None"])
    lines.extend(["", "## Locked Interpretation", str(cert.get("locked_interpretation", ""))])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, map_root: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    map_root.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v37_mechanism_question_probes_certificate.json"
    report_path = out_dir / "V37_MECHANISM_QUESTION_PROBES_REPORT.md"
    decision_path = out_dir / "v37_mechanism_question_probes_next_decision.json"
    target_results_path = out_dir / "v37_target_mechanism_probe_results.json"
    answer_table_path = map_root / "V37_MECHANISM_QUESTION_ANSWER_TABLE.md"
    map_paths = {
        "KcsA": map_root / "KcsA_mechanism_question_map.json",
        "XCL1_lymphotactin": map_root / "XCL1_lymphotactin_mechanism_question_map.json",
        "alpha_synuclein_SNCA": map_root / "alpha_synuclein_SNCA_mechanism_question_map.json",
    }
    decision = {
        "kind": "V37_MECHANISM_QUESTION_PROBES_NEXT_DECISION_v0",
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
            "answer_table": str(answer_table_path),
            "target_maps": {target: str(path) for target, path in map_paths.items()},
        },
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, decision)
    _write_json(target_results_path, cert.get("target_results", []))
    for row in cert.get("target_results", []):
        _write_json(map_paths[row["target"]], _target_map(row))
    _write_answer_table(answer_table_path, cert)
    _write_report(report_path, cert)
    paths = {"certificate": cert_path, "report": report_path, "decision": decision_path, "target_results": target_results_path, "answer_table": answer_table_path}
    paths.update({f"{target}_map": path for target, path in map_paths.items()})
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V37 mechanism question probes.")
    parser.add_argument("--data-root", type=Path, default=DATA_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--map-root", type=Path, default=MAP_ROOT)
    args = parser.parse_args()
    packages = load_dossiers(args.data_root)
    cert = build_v37(packages)
    paths = write_outputs(args.out_dir, args.map_root, cert)
    print(json.dumps({
        "kind": cert.get("kind"),
        "control_status": cert.get("control_status"),
        "input_v36_status": cert.get("input_v36_status"),
        "mechanism_classes_assigned": cert.get("mechanism_classes_assigned"),
        "assigned_target_count": cert.get("assigned_target_count"),
        "partial_target_count": cert.get("partial_target_count"),
        "blocked_target_count": cert.get("blocked_target_count"),
        "target_name_masking_passed": cert.get("target_name_masking_passed"),
        "swapped_dossier_detection_passed": cert.get("swapped_dossier_detection_passed"),
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
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
        "answer_table": str(paths["answer_table"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
