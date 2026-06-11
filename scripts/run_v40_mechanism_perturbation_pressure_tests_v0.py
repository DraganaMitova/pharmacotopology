#!/usr/bin/env python3
from __future__ import annotations

"""Run V40 mechanism perturbation pressure validation."""

import argparse
import copy
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "mechanism_perturbations" / "V40"
SOURCE_ROOT = DATA_ROOT / "sources"
PRED_ROOT = DATA_ROOT / "predictions"
VALIDATION_ROOT = DATA_ROOT / "validation"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V40_MECHANISM_PERTURBATION_PRESSURE_TESTS"

TARGET_ORDER = ["KcsA", "XCL1_lymphotactin", "alpha_synuclein_SNCA"]
REQUIRED_BUCKETS_BY_TARGET = {
    "KcsA": {
        "filter_signature_perturbation",
        "ion_selectivity_perturbation",
        "membrane_topology_context_perturbation",
        "oligomer_interface_context_perturbation",
    },
    "XCL1_lymphotactin": {
        "state_A_loss_or_weakening",
        "state_B_loss_or_weakening",
        "mixed_state_pooling_error",
        "state_function_decoupling",
        "two_state_balance_shift",
    },
    "alpha_synuclein_SNCA": {
        "disorder_evidence_loss",
        "ensemble_to_single_fold_forcing",
        "context_bound_helix_overpromotion",
        "aggregation_context_overpromotion",
        "disorder_to_order_context_shift",
    },
}

PASSED = "V40_MECHANISM_PERTURBATION_PRESSURE_PASSED_CLAIM_DISABLED"
PARTIAL = "V40_PARTIAL_PERTURBATION_SUPPORT_CLEAN_ABSTAIN"
BLOCKED_LEAKAGE = "V40_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"
BLOCKED_KEY = "V40_BLOCKED_ANSWER_KEY_LEAKAGE"
FAILED = "V40_FAILED_PERTURBATION_PRESSURE_TEST"

SUPPORTED = "perturbation_supported_by_holdout"
PARTIAL_LEVEL = "perturbation_partially_supported_clean_abstain"
CONTRADICTED = "perturbation_contradicted_or_failed"
BLOCKED = "blocked_for_leakage"

PLACEHOLDER_TOKENS = {"placeholder", "todo", "tbd", "citation needed", "example.com", "dummy", "unknown"}
COORDINATE_TOKENS = {
    "1bl8",
    "pdb coordinate",
    "coordinate-derived",
    "coordinate_derived",
    "alphafold",
    "esmfold",
    "rosettafold",
    "rmsd",
    "tm-score",
    ".pdb",
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


def _row_text(row: dict[str, Any]) -> str:
    fields = [
        "source_id",
        "source_role",
        "source_type",
        "source_name",
        "source_url_or_citation",
        "evidence_statement",
    ]
    return " ".join(str(row.get(field, "")) for field in fields).lower()


def _source_failures(rows: list[dict[str, Any]], prefix: str) -> dict[str, Any]:
    coord = 0
    internal = 0
    key = 0
    placeholder = 0
    failed: list[str] = []
    for index, row in enumerate(rows):
        text = _row_text(row)
        if (
            _bool_true(row.get("coordinate_derived"))
            or _bool_true(row.get("coordinate_truth_used_before_selection"))
            or _bool_true(row.get("native_metrics_used_for_selection"))
            or any(token in text for token in COORDINATE_TOKENS)
        ):
            coord += 1
            failed.append(f"{prefix}:row_{index}:coordinate_or_native_metric_leakage")
        if _bool_true(row.get("internal_runtime_source")) or "first_contact_clean_pharmacotopology_layer_run" in text:
            internal += 1
            failed.append(f"{prefix}:row_{index}:internal_runtime_source")
        if _bool_true(row.get("answer_key_source")) or "answer_key" in text or "blind_mechanism_panels/v38/answer_key" in text:
            key += 1
            failed.append(f"{prefix}:row_{index}:answer_key_leakage")
        for field in ["source_name", "source_url_or_citation", "source_date_or_version"]:
            value = str(row.get(field, "")).lower()
            if not value or any(token in value for token in PLACEHOLDER_TOKENS):
                placeholder += 1
                failed.append(f"{prefix}:row_{index}:placeholder_{field}")
    return {
        "coordinate_derived_source_count": coord,
        "internal_runtime_source_count": internal,
        "answer_key_source_count": key,
        "placeholder_source_count": placeholder,
        "failed_checks": sorted(set(failed)),
    }


def validate_target(packet: dict[str, Any], manifest: dict[str, Any]) -> dict[str, Any]:
    target = str(packet.get("target"))
    manifest_target = str(manifest.get("target"))
    perturbations = [row for row in _as_list(packet.get("perturbation_predictions")) if isinstance(row, dict)]
    prediction_sources = [row for row in _as_list(manifest.get("prediction_seed_sources")) if isinstance(row, dict)]
    validation_sources = [row for row in _as_list(manifest.get("validation_holdout_sources")) if isinstance(row, dict)]
    failed: list[str] = []

    if manifest_target != target:
        failed.append(f"{target}:manifest_target_mismatch:{manifest_target}")
    if _bool_true(packet.get("answer_key_used_for_prediction")) or _bool_true(manifest.get("answer_key_used_for_prediction")):
        failed.append(f"{target}:answer_key_used_for_prediction")
    for flag in ["claim_allowed", "positive_folding_evidence_found", "folding_problem_solved", "new_md_executed"]:
        if _bool_true(packet.get(flag)) or _bool_true(manifest.get(flag)):
            failed.append(f"{target}:forbidden_flag_enabled:{flag}")

    pred_fail = _source_failures(prediction_sources, f"{target}:prediction_source")
    holdout_fail = _source_failures(validation_sources, f"{target}:validation_holdout")
    failed.extend(pred_fail["failed_checks"])
    failed.extend(holdout_fail["failed_checks"])

    prediction_source_ids = {str(row.get("source_id")) for row in prediction_sources}
    validation_source_ids = {str(row.get("source_id")) for row in validation_sources}
    overlap = sorted(prediction_source_ids & validation_source_ids)
    if overlap:
        failed.append(f"{target}:prediction_validation_source_overlap:{','.join(overlap)}")
    for row in validation_sources:
        if row.get("independent_from_prediction_sources") is not True:
            failed.append(f"{target}:validation_holdout_not_independent:{row.get('source_id')}")

    expected_buckets = {
        str(row.get("expected_validation_bucket") or row.get("perturbation_bucket"))
        for row in perturbations
        if row.get("expected_validation_bucket") or row.get("perturbation_bucket")
    }
    required_buckets = REQUIRED_BUCKETS_BY_TARGET.get(target, expected_buckets)
    observed_buckets = {str(row.get("evidence_bucket")) for row in validation_sources if row.get("evidence_bucket")}
    supported = sorted(expected_buckets & observed_buckets)
    missing_holdout = sorted(expected_buckets - observed_buckets)
    missing_required = sorted(required_buckets - expected_buckets)
    if missing_required:
        failed.append(f"{target}:missing_required_perturbation_buckets:{','.join(missing_required)}")

    coord = pred_fail["coordinate_derived_source_count"] + holdout_fail["coordinate_derived_source_count"]
    internal = pred_fail["internal_runtime_source_count"] + holdout_fail["internal_runtime_source_count"]
    key = pred_fail["answer_key_source_count"] + holdout_fail["answer_key_source_count"] + (1 if f"{target}:answer_key_used_for_prediction" in failed else 0)
    placeholder = pred_fail["placeholder_source_count"] + holdout_fail["placeholder_source_count"]

    if coord or internal or key or placeholder:
        level = BLOCKED
    elif not perturbations or len(expected_buckets) < 3 or len(supported) == 0:
        level = CONTRADICTED
    elif missing_holdout or missing_required:
        level = PARTIAL_LEVEL
    else:
        level = SUPPORTED

    return {
        "target": target,
        "mechanism_class": packet.get("mechanism_class"),
        "scientist_question": packet.get("scientist_question"),
        "perturbation_count": len(perturbations),
        "prediction_source_count": len(prediction_sources),
        "validation_holdout_source_count": len(validation_sources),
        "supported_perturbation_buckets": supported,
        "missing_holdout_buckets": missing_holdout,
        "missing_required_perturbation_buckets": missing_required,
        "validation_level": level,
        "coordinate_derived_source_count": coord,
        "internal_runtime_source_count": internal,
        "answer_key_source_count": key,
        "placeholder_source_count": placeholder,
        "failed_checks": sorted(set(failed)),
    }


def load_inputs() -> dict[str, tuple[dict[str, Any], dict[str, Any]]]:
    return {
        target: (
            _read_json(PRED_ROOT / target / "perturbation_prediction_packet.json", f"{target} perturbation packet"),
            _read_json(SOURCE_ROOT / target / "perturbation_source_manifest.json", f"{target} source manifest"),
        )
        for target in TARGET_ORDER
    }


def _evaluate_inputs(inputs: dict[str, tuple[dict[str, Any], dict[str, Any]]]) -> list[dict[str, Any]]:
    return [validate_target(packet, manifest) for packet, manifest in inputs.values()]


def _mutate(inputs: dict[str, tuple[dict[str, Any], dict[str, Any]]], target: str, packet_fn=None, manifest_fn=None) -> dict[str, tuple[dict[str, Any], dict[str, Any]]]:
    out = copy.deepcopy(inputs)
    packet, manifest = out[target]
    if packet_fn:
        packet_fn(packet)
    if manifest_fn:
        manifest_fn(manifest)
    return out


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _observed(result: dict[str, Any]) -> dict[str, Any]:
    return {
        "target": result["target"],
        "validation_level": result["validation_level"],
        "supported_perturbation_buckets": result["supported_perturbation_buckets"],
        "missing_holdout_buckets": result["missing_holdout_buckets"],
        "missing_required_perturbation_buckets": result["missing_required_perturbation_buckets"],
        "coordinate_derived_source_count": result["coordinate_derived_source_count"],
        "internal_runtime_source_count": result["internal_runtime_source_count"],
        "answer_key_source_count": result["answer_key_source_count"],
        "placeholder_source_count": result["placeholder_source_count"],
    }


def _generic_prediction(bucket: str) -> dict[str, str]:
    return {
        "perturbation_id": "adversarial_control",
        "operator": "adversarial_control",
        "perturbation_bucket": bucket,
        "expected_validation_bucket": bucket,
        "predicted_pressure": "block",
        "expected_behavior": "adversarial control should not validate",
        "falsifiable_if": "adversarial control validates",
    }


def run_controls(inputs: dict[str, tuple[dict[str, Any], dict[str, Any]]]) -> list[dict[str, Any]]:
    baseline = _evaluate_inputs(inputs)
    coord_count = sum(row["coordinate_derived_source_count"] for row in baseline)
    controls = [
        _control("prediction_perturbations_use_no_coordinate_sources", coord_count == 0, "Prediction perturbations must use no coordinate-derived evidence.", {"coordinate_derived_source_count": coord_count}),
        _control("validation_holdouts_use_no_coordinate_sources", all(row["validation_level"] != BLOCKED for row in baseline), "Validation holdouts must use no coordinate-derived evidence.", {"validation_levels": {row["target"]: row["validation_level"] for row in baseline}}),
    ]

    internal = _mutate(
        inputs,
        "KcsA",
        manifest_fn=lambda manifest: manifest["validation_holdout_sources"].append({
            "source_id": "bad_internal_runtime",
            "source_name": "runtime report",
            "source_url_or_citation": "first_contact_clean_pharmacotopology_layer_run/V39/report.json",
            "source_date_or_version": "accessed_2026-06-11",
            "internal_runtime_source": True,
        }),
    )
    internal_result = validate_target(*internal["KcsA"])
    controls.append(_control("internal_runtime_sources_are_blocked", internal_result["validation_level"] == BLOCKED and internal_result["internal_runtime_source_count"] > 0, "Internal runtime sources must block validation.", _observed(internal_result)))

    key = _mutate(inputs, "KcsA", packet_fn=lambda packet: packet.update({"answer_key_used_for_prediction": True}))
    key_result = validate_target(*key["KcsA"])
    controls.append(_control("v38_answer_key_leakage_attempt_is_blocked", key_result["validation_level"] == BLOCKED and key_result["answer_key_source_count"] > 0, "V38 answer-key leakage must block validation.", _observed(key_result)))

    kcsa_generic = _mutate(
        inputs,
        "KcsA",
        packet_fn=lambda packet: packet.update({"perturbation_predictions": [_generic_prediction("generic_channel_annotation_only")]}),
    )
    kcsa_generic_result = validate_target(*kcsa_generic["KcsA"])
    controls.append(_control("kcsa_generic_channel_only_does_not_preserve_hard_grammar", kcsa_generic_result["validation_level"] != SUPPORTED, "Generic channel annotation alone must not preserve hard KcsA grammar.", _observed(kcsa_generic_result)))

    kcsa_removed = _mutate(
        inputs,
        "KcsA",
        packet_fn=lambda packet: packet.update({
            "perturbation_predictions": [
                row for row in packet["perturbation_predictions"]
                if row["perturbation_id"] not in {"KCSA_V40_P1_FILTER_SIGNATURE", "KCSA_V40_P2_ION_SELECTIVITY"}
            ]
        }),
    )
    kcsa_removed_result = validate_target(*kcsa_removed["KcsA"])
    controls.append(_control("kcsa_filter_ion_selectivity_removed_weakens_or_invalidates", kcsa_removed_result["validation_level"] in {PARTIAL_LEVEL, CONTRADICTED}, "Removing KcsA filter/ion-selectivity perturbations must weaken or invalidate hard grammar.", _observed(kcsa_removed_result)))

    xcl1_state_a = _mutate(
        inputs,
        "XCL1_lymphotactin",
        packet_fn=lambda packet: packet.update({
            "perturbation_predictions": [
                row for row in packet["perturbation_predictions"] if row["perturbation_id"] != "XCL1_V40_P1_STATE_A_LOSS"
            ]
        }),
    )
    xcl1_state_a_result = validate_target(*xcl1_state_a["XCL1_lymphotactin"])
    controls.append(_control("xcl1_state_a_removed_becomes_partial_or_invalid", xcl1_state_a_result["validation_level"] in {PARTIAL_LEVEL, CONTRADICTED}, "Removing XCL1 state A must become partial/invalid.", _observed(xcl1_state_a_result)))

    xcl1_state_b = _mutate(
        inputs,
        "XCL1_lymphotactin",
        packet_fn=lambda packet: packet.update({
            "perturbation_predictions": [
                row for row in packet["perturbation_predictions"] if row["perturbation_id"] != "XCL1_V40_P2_STATE_B_LOSS"
            ]
        }),
    )
    xcl1_state_b_result = validate_target(*xcl1_state_b["XCL1_lymphotactin"])
    controls.append(_control("xcl1_state_b_removed_becomes_partial_or_invalid", xcl1_state_b_result["validation_level"] in {PARTIAL_LEVEL, CONTRADICTED}, "Removing XCL1 state B must become partial/invalid.", _observed(xcl1_state_b_result)))

    xcl1_pool = _mutate(
        inputs,
        "XCL1_lymphotactin",
        packet_fn=lambda packet: packet["perturbation_predictions"].append(_generic_prediction("mixed_state_pooling_allowed")),
    )
    xcl1_pool_result = validate_target(*xcl1_pool["XCL1_lymphotactin"])
    controls.append(_control("xcl1_mixed_state_pooling_forced_is_blocked", xcl1_pool_result["validation_level"] != SUPPORTED and "mixed_state_pooling_allowed" in xcl1_pool_result["missing_holdout_buckets"], "Forced mixed-state pooling must be blocked.", _observed(xcl1_pool_result)))

    snca_disorder_removed = _mutate(
        inputs,
        "alpha_synuclein_SNCA",
        packet_fn=lambda packet: packet.update({
            "perturbation_predictions": [
                row for row in packet["perturbation_predictions"] if row["perturbation_id"] != "SNCA_V40_P1_DISORDER_LOSS"
            ]
        }),
    )
    snca_disorder_result = validate_target(*snca_disorder_removed["alpha_synuclein_SNCA"])
    controls.append(_control("snca_disorder_evidence_removed_becomes_partial_or_invalid", snca_disorder_result["validation_level"] in {PARTIAL_LEVEL, CONTRADICTED}, "Removing SNCA disorder evidence must become partial/invalid.", _observed(snca_disorder_result)))

    snca_single = _mutate(
        inputs,
        "alpha_synuclein_SNCA",
        packet_fn=lambda packet: packet["perturbation_predictions"].append(_generic_prediction("compact_single_fold_native_claim")),
    )
    snca_single_result = validate_target(*snca_single["alpha_synuclein_SNCA"])
    controls.append(_control("snca_compact_single_fold_forcing_is_blocked", snca_single_result["validation_level"] != SUPPORTED and "compact_single_fold_native_claim" in snca_single_result["missing_holdout_buckets"], "SNCA compact single-fold forcing must be blocked.", _observed(snca_single_result)))

    snca_aggregation = _mutate(
        inputs,
        "alpha_synuclein_SNCA",
        packet_fn=lambda packet: packet["perturbation_predictions"].append(_generic_prediction("aggregation_native_fold_claim")),
    )
    snca_aggregation_result = validate_target(*snca_aggregation["alpha_synuclein_SNCA"])
    controls.append(_control("snca_aggregation_context_overpromoted_to_native_fold_is_blocked", snca_aggregation_result["validation_level"] != SUPPORTED and "aggregation_native_fold_claim" in snca_aggregation_result["missing_holdout_buckets"], "SNCA aggregation context must not become a native-fold claim.", _observed(snca_aggregation_result)))

    swapped = _mutate(inputs, "KcsA", manifest_fn=lambda manifest: manifest.update(inputs["alpha_synuclein_SNCA"][1]))
    swapped_result = validate_target(*swapped["KcsA"])
    controls.append(_control("swapped_perturbation_holdouts_fail_validation", swapped_result["validation_level"] != SUPPORTED, "Swapped perturbation holdouts must be detected or fail validation.", _observed(swapped_result)))

    placeholder = _mutate(
        inputs,
        "KcsA",
        manifest_fn=lambda manifest: manifest["validation_holdout_sources"].append({
            "source_id": "placeholder",
            "source_name": "TODO placeholder",
            "source_url_or_citation": "citation needed",
            "source_date_or_version": "accessed_2026-06-11",
        }),
    )
    placeholder_result = validate_target(*placeholder["KcsA"])
    controls.append(_control("placeholder_citations_source_rows_are_blocked", placeholder_result["validation_level"] == BLOCKED and placeholder_result["placeholder_source_count"] > 0, "Placeholder citations/source rows must block validation.", _observed(placeholder_result)))
    return controls


def build_v40() -> dict[str, Any]:
    inputs = load_inputs()
    results = _evaluate_inputs(inputs)
    controls = run_controls(inputs)
    supported = [row for row in results if row["validation_level"] == SUPPORTED]
    partial = [row for row in results if row["validation_level"] == PARTIAL_LEVEL]
    contradicted = [row for row in results if row["validation_level"] == CONTRADICTED]
    blocked = [row for row in results if row["validation_level"] == BLOCKED]
    coord = sum(row["coordinate_derived_source_count"] for row in results)
    internal = sum(row["internal_runtime_source_count"] for row in results)
    key = sum(row["answer_key_source_count"] for row in results)
    answer_key_used = any(
        _bool_true(packet.get("answer_key_used_for_prediction")) or _bool_true(manifest.get("answer_key_used_for_prediction"))
        for packet, manifest in inputs.values()
    )
    failed = sorted(set(check for row in results for check in row["failed_checks"]))
    if not all(control["passed"] for control in controls):
        failed.append("V40_required_controls_failed")

    enough_buckets = all(row["perturbation_count"] >= 3 for row in results)
    if key or answer_key_used:
        status = BLOCKED_KEY
    elif coord or internal or blocked:
        status = BLOCKED_LEAKAGE
    elif contradicted or not all(control["passed"] for control in controls) or not enough_buckets:
        status = FAILED
    elif len(supported) >= 2:
        status = PASSED
    else:
        status = PARTIAL

    if status == PASSED:
        next_action = "advance_to_unseen_target_perturbation_grammar_no_claim_no_MD"
    elif status == PARTIAL:
        next_action = "collect_more_independent_perturbation_holdouts_before_generalization"
    elif status == BLOCKED_KEY:
        next_action = "remove_answer_key_leakage_before_V40"
    elif status == BLOCKED_LEAKAGE:
        next_action = "remove_coordinate_or_internal_leakage_before_V40"
    else:
        next_action = "repair_failed_perturbation_pressure_controls_before_next_step"

    return {
        "kind": "V40_MECHANISM_PERTURBATION_PRESSURE_TESTS_v0",
        "run_mode": "mechanism_perturbation_pressure_no_coordinates_no_MD_claim_disabled",
        "control_status": status,
        "target_count": len(TARGET_ORDER),
        "perturbation_packet_count": len(results),
        "perturbation_count_by_target": {row["target"]: row["perturbation_count"] for row in results},
        "validated_perturbation_count_by_target": {row["target"]: len(row["supported_perturbation_buckets"]) for row in results},
        "validation_level_by_target": {row["target"]: row["validation_level"] for row in results},
        "supported_target_count": len(supported),
        "partial_target_count": len(partial),
        "contradicted_target_count": len(contradicted),
        "blocked_target_count": len(blocked),
        "coordinate_derived_source_count": coord,
        "internal_runtime_source_count": internal,
        "answer_key_used_for_prediction": answer_key_used,
        "native_metrics_used_for_selection": False,
        "coordinate_truth_used_before_selection": False,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
        "control_count": len(controls),
        "passed_control_count": sum(1 for control in controls if control["passed"]),
        "failed_checks": failed,
        "next_action": next_action,
        "controls": controls,
        "target_validation_results": results,
        "locked_interpretation": (
            "V40 passing means the non-coordinate mechanism grammar identifies causal perturbation pressure for KcsA, XCL1, and SNCA: "
            "which operators should break, weaken, shift, or block the mechanism, while keeping coordinates, MD, answer-key leakage, and folding-solved claims disabled."
        ),
    }


def _answer_text(result: dict[str, Any]) -> str:
    target = result["target"]
    supported = ", ".join(result["supported_perturbation_buckets"]) or "none"
    if target == "KcsA":
        answer = "Filter and ion-selectivity perturbations are causal pressure points; generic channel annotation alone is insufficient."
    elif target == "XCL1_lymphotactin":
        answer = "State separation is causal pressure; deleting either state or pooling states weakens or invalidates the metamorphic mechanism."
    else:
        answer = "Disorder evidence and context boundaries are causal pressure; compact single-fold and native-fold overpromotion are blocked."
    return f"{answer} Supported perturbation buckets: {supported}."


def _write_scientist_answer(path: Path, result: dict[str, Any]) -> None:
    lines = [
        f"# V40 Scientist Question Answer: {result['target']}",
        "",
        f"Question: {result['scientist_question']}",
        "",
        f"Validation level: `{result['validation_level']}`",
        "",
        _answer_text(result),
        "",
        "This remains mechanism perturbation pressure, not atomistic structure prediction and not a solved-folding claim.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V40 Mechanism Perturbation Pressure Tests",
        "",
        f"Status: `{cert['control_status']}`",
        f"Supported targets: `{cert['supported_target_count']}` / `{cert['target_count']}`",
        f"Controls: `{cert['passed_control_count']}` / `{cert['control_count']}`",
        f"Coordinate-derived source count: `{cert['coordinate_derived_source_count']}`",
        f"Internal runtime source count: `{cert['internal_runtime_source_count']}`",
        f"Answer key used for prediction: `{cert['answer_key_used_for_prediction']}`",
        "",
        "## Target Summary",
    ]
    for result in cert["target_validation_results"]:
        lines.extend([
            f"### {result['target']}",
            f"Validation level: `{result['validation_level']}`",
            f"Perturbations: `{result['perturbation_count']}`",
            f"Supported perturbation buckets: `{result['supported_perturbation_buckets']}`",
            f"Scientist answer: {_answer_text(result)}",
            "",
        ])
    lines.extend(["## Plain English Interpretation", cert["locked_interpretation"]])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    cert = build_v40()
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v40_mechanism_perturbation_pressure_tests_certificate.json"
    report_path = out_dir / "V40_MECHANISM_PERTURBATION_PRESSURE_TESTS_REPORT.md"
    decision_path = out_dir / "v40_mechanism_perturbation_pressure_tests_next_decision.json"
    decision = {
        "kind": "V40_MECHANISM_PERTURBATION_PRESSURE_TESTS_NEXT_DECISION_v0",
        "decision_status": cert["control_status"],
        "next_action": cert["next_action"],
        "claim_allowed": False,
        "new_MD_allowed": False,
        "folding_problem_solved": False,
    }
    cert = {
        **cert,
        "next_decision": decision,
        "artifacts": {"certificate": str(cert_path), "report": str(report_path), "decision": str(decision_path)},
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, decision)
    _write_report(report_path, cert)
    for result in cert["target_validation_results"]:
        target_dir = VALIDATION_ROOT / result["target"]
        _write_json(target_dir / "perturbation_validation_result.json", result)
        _write_scientist_answer(target_dir / "scientist_question_answer.md", result)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V40 mechanism perturbation pressure tests.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = write_outputs(args.out_dir)
    cert = _read_json(paths["certificate"], "V40 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "control_status": cert["control_status"],
        "validation_level_by_target": cert["validation_level_by_target"],
        "perturbation_count_by_target": cert["perturbation_count_by_target"],
        "validated_perturbation_count_by_target": cert["validated_perturbation_count_by_target"],
        "control_count": cert["control_count"],
        "passed_control_count": cert["passed_control_count"],
        "coordinate_derived_source_count": cert["coordinate_derived_source_count"],
        "internal_runtime_source_count": cert["internal_runtime_source_count"],
        "answer_key_used_for_prediction": cert["answer_key_used_for_prediction"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
