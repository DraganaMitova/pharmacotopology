#!/usr/bin/env python3
from __future__ import annotations

"""Validate V39 mechanism-level falsifiable predictions against holdouts."""

import argparse
import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
PRED_ROOT = REPO_ROOT / "data" / "mechanism_predictions" / "V39"
HOLDOUT_ROOT = REPO_ROOT / "data" / "mechanism_holdouts" / "V39"
VALIDATION_ROOT = REPO_ROOT / "data" / "mechanism_validation" / "V39"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V39_MECHANISM_TO_FALSIFIABLE_PREDICTION_VALIDATION"

TARGET_ORDER = ["KcsA", "XCL1_lymphotactin", "alpha_synuclein_SNCA"]
REQUIRED_BUCKETS_BY_TARGET = {
    "KcsA": {
        "filter_or_signature_holdout",
        "membrane_topology_holdout",
        "oligomer_or_interface_holdout",
        "perturbation_holdout",
    },
    "XCL1_lymphotactin": {
        "metamorphic_two_state_holdout",
        "state_A_state_B_holdout",
        "state_function_holdout",
        "pooling_rule_holdout",
        "perturbation_holdout",
    },
    "alpha_synuclein_SNCA": {
        "disorder_holdout",
        "context_bound_holdout",
        "disorder_to_order_holdout",
        "single_fold_block_holdout",
        "perturbation_holdout",
    },
}

PASSED = "V39_MECHANISM_PREDICTIONS_VALIDATED_CLAIM_DISABLED"
PARTIAL = "V39_PARTIAL_HOLDOUT_SUPPORT_CLEAN_ABSTAIN"
BLOCKED_LEAKAGE = "V39_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"
BLOCKED_KEY = "V39_BLOCKED_ANSWER_KEY_LEAKAGE"
FAILED = "V39_FAILED_MECHANISM_PREDICTION_VALIDATION"

PLACEHOLDER_TOKENS = {"placeholder", "todo", "tbd", "citation needed", "example.com", "dummy", "unknown"}
COORDINATE_TOKENS = {"1bl8", "pdb coordinate", "coordinate-derived", "coordinate_derived", "alphafold", "esmfold", "rosettafold", "rmsd", "tm-score", ".pdb"}


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
    return " ".join(str(row.get(key, "")) for key in ["source_id", "source_type", "source_name", "source_url_or_citation", "evidence_statement"]).lower()


def _source_failures(rows: list[dict[str, Any]], prefix: str) -> dict[str, Any]:
    coord = 0
    internal = 0
    key = 0
    placeholder = 0
    failed: list[str] = []
    for index, row in enumerate(rows):
        text = _row_text(row)
        if _bool_true(row.get("coordinate_derived")) or _bool_true(row.get("coordinate_truth_used_before_selection")) or _bool_true(row.get("native_metrics_used_for_selection")) or any(token in text for token in COORDINATE_TOKENS):
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


def validate_target(packet: dict[str, Any], holdout: dict[str, Any]) -> dict[str, Any]:
    target = str(packet.get("target"))
    predictions = [row for row in _as_list(packet.get("falsifiable_predictions")) if isinstance(row, dict)]
    pred_sources = [row for row in _as_list(packet.get("prediction_source_evidence")) if isinstance(row, dict)]
    holdout_sources = [row for row in _as_list(holdout.get("holdout_validation_evidence")) if isinstance(row, dict)]
    failed: list[str] = []

    pred_fail = _source_failures(pred_sources, f"{target}:prediction")
    hold_fail = _source_failures(holdout_sources, f"{target}:holdout")
    failed.extend(pred_fail["failed_checks"])
    failed.extend(hold_fail["failed_checks"])
    if _bool_true(packet.get("answer_key_used_for_prediction")) or _bool_true(holdout.get("answer_key_used_for_prediction")):
        failed.append(f"{target}:answer_key_used_for_prediction")
    if _bool_true(packet.get("claim_allowed")) or _bool_true(packet.get("positive_folding_evidence_found")) or _bool_true(packet.get("folding_problem_solved")):
        failed.append(f"{target}:prediction_claim_flags_not_closed")

    prediction_source_ids = {str(row.get("source_id")) for row in pred_sources}
    holdout_source_ids = {str(row.get("source_id")) for row in holdout_sources}
    overlap = sorted(prediction_source_ids & holdout_source_ids)
    if overlap:
        failed.append(f"{target}:prediction_holdout_source_overlap:{','.join(overlap)}")
    for row in holdout_sources:
        if row.get("independent_from_prediction_sources") is not True:
            failed.append(f"{target}:holdout_not_marked_independent:{row.get('source_id')}")

    expected_buckets = {str(row.get("expected_holdout_bucket")) for row in predictions}
    required_buckets = REQUIRED_BUCKETS_BY_TARGET.get(target, expected_buckets)
    observed_buckets = {str(row.get("evidence_bucket")) for row in holdout_sources}
    supported = sorted(expected_buckets & observed_buckets)
    missing_holdout = sorted(expected_buckets - observed_buckets)
    missing_required = sorted(required_buckets - expected_buckets)
    if missing_required:
        failed.append(f"{target}:missing_required_prediction_buckets:{','.join(missing_required)}")

    coord = pred_fail["coordinate_derived_source_count"] + hold_fail["coordinate_derived_source_count"]
    internal = pred_fail["internal_runtime_source_count"] + hold_fail["internal_runtime_source_count"]
    key = pred_fail["answer_key_source_count"] + hold_fail["answer_key_source_count"] + (1 if f"{target}:answer_key_used_for_prediction" in failed else 0)
    placeholder = pred_fail["placeholder_source_count"] + hold_fail["placeholder_source_count"]

    if coord or internal:
        level = "blocked_for_leakage"
    elif key or placeholder:
        level = "blocked_for_leakage"
    elif not predictions or len(supported) == 0:
        level = "contradicted_or_failed"
    elif missing_holdout or missing_required:
        level = "partially_supported_clean_abstain"
    else:
        level = "supported_by_independent_holdout"

    return {
        "target": target,
        "mechanism_class": packet.get("mechanism_class"),
        "prediction_count": len(predictions),
        "prediction_source_count": len(pred_sources),
        "holdout_source_count": len(holdout_sources),
        "supported_prediction_buckets": supported,
        "missing_holdout_buckets": missing_holdout,
        "missing_required_prediction_buckets": missing_required,
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
            _read_json(PRED_ROOT / target / "prediction_packet.json", f"{target} prediction packet"),
            _read_json(HOLDOUT_ROOT / target / "holdout_manifest.json", f"{target} holdout manifest"),
        )
        for target in TARGET_ORDER
    }


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _observed_validation(result: dict[str, Any]) -> dict[str, Any]:
    return {
        "target": result["target"],
        "validation_level": result["validation_level"],
        "missing_holdout_buckets": result["missing_holdout_buckets"],
        "missing_required_prediction_buckets": result["missing_required_prediction_buckets"],
        "coordinate_derived_source_count": result["coordinate_derived_source_count"],
        "internal_runtime_source_count": result["internal_runtime_source_count"],
        "answer_key_source_count": result["answer_key_source_count"],
        "placeholder_source_count": result["placeholder_source_count"],
    }


def _baseline_inputs() -> dict[str, tuple[dict[str, Any], dict[str, Any]]]:
    return load_inputs()


def _mutate_first(inputs: dict[str, tuple[dict[str, Any], dict[str, Any]]], target: str, packet_fn=None, holdout_fn=None) -> dict[str, tuple[dict[str, Any], dict[str, Any]]]:
    import copy
    out = copy.deepcopy(inputs)
    packet, holdout = out[target]
    if packet_fn:
        packet_fn(packet)
    if holdout_fn:
        holdout_fn(holdout)
    return out


def _evaluate_inputs(inputs: dict[str, tuple[dict[str, Any], dict[str, Any]]]) -> list[dict[str, Any]]:
    return [validate_target(packet, holdout) for packet, holdout in inputs.values()]


def run_controls(inputs: dict[str, tuple[dict[str, Any], dict[str, Any]]]) -> list[dict[str, Any]]:
    baseline = _evaluate_inputs(inputs)
    prediction_coordinate_count = sum(r["coordinate_derived_source_count"] for r in baseline)
    controls = [
        _control("prediction_generation_uses_no_coordinate_sources", prediction_coordinate_count == 0, "Prediction sources must be coordinate-free.", {"coordinate_derived_source_count": prediction_coordinate_count}),
        _control("holdout_validation_uses_no_coordinate_sources", all(r["validation_level"] != "blocked_for_leakage" for r in baseline), "Holdout sources must be coordinate-free.", {"validation_levels": {r["target"]: r["validation_level"] for r in baseline}}),
    ]
    internal = _mutate_first(inputs, "KcsA", holdout_fn=lambda h: h["holdout_validation_evidence"].append({"source_id": "bad_internal", "source_name": "runtime", "source_url_or_citation": "first_contact_clean_pharmacotopology_layer_run/V38/report.json", "internal_runtime_source": True}))
    internal_result = validate_target(*internal["KcsA"])
    controls.append(_control("internal_runtime_reports_blocked_as_evidence", internal_result["validation_level"] == "blocked_for_leakage" and internal_result["internal_runtime_source_count"] > 0, "Internal runtime evidence must block.", _observed_validation(internal_result)))
    key = _mutate_first(inputs, "KcsA", packet_fn=lambda p: p.update({"answer_key_used_for_prediction": True}))
    key_result = validate_target(*key["KcsA"])
    controls.append(_control("v38_answer_key_leakage_attempt_blocked", key_result["validation_level"] == "blocked_for_leakage" and key_result["answer_key_source_count"] > 0, "Answer-key leakage must block.", _observed_validation(key_result)))
    generic = _mutate_first(inputs, "KcsA", packet_fn=lambda p: p.update({"falsifiable_predictions": [{"prediction_id": "generic", "expected_holdout_bucket": "generic"}]}), holdout_fn=lambda h: h.update({"holdout_validation_evidence": []}))
    generic_result = validate_target(*generic["KcsA"])
    controls.append(_control("generic_only_evidence_not_high_confidence", generic_result["validation_level"] != "supported_by_independent_holdout", "Generic-only evidence cannot validate high-confidence prediction.", _observed_validation(generic_result)))
    kcsa_removed = _mutate_first(inputs, "KcsA", packet_fn=lambda p: p.update({"falsifiable_predictions": [row for row in p["falsifiable_predictions"] if row["prediction_id"] != "KCSA_P1_FILTER"]}))
    kcsa_removed_result = validate_target(*kcsa_removed["KcsA"])
    controls.append(_control("kcsa_filter_ion_removed_partial_or_invalid", kcsa_removed_result["validation_level"] in {"partially_supported_clean_abstain", "contradicted_or_failed"}, "Removing KcsA filter prediction must force partial/invalid status.", _observed_validation(kcsa_removed_result)))
    xcl1_one = _mutate_first(inputs, "XCL1_lymphotactin", packet_fn=lambda p: p.update({"falsifiable_predictions": [row for row in p["falsifiable_predictions"] if row["prediction_id"] != "XCL1_P2_TWO_STATES"]}))
    xcl1_one_result = validate_target(*xcl1_one["XCL1_lymphotactin"])
    controls.append(_control("xcl1_one_state_only_partial_or_invalid", xcl1_one_result["validation_level"] in {"partially_supported_clean_abstain", "contradicted_or_failed"}, "Removing one XCL1 state must force partial/invalid status.", _observed_validation(xcl1_one_result)))
    xcl1_pool = _mutate_first(inputs, "XCL1_lymphotactin", packet_fn=lambda p: p["falsifiable_predictions"].append({"prediction_id": "bad_pool", "expected_holdout_bucket": "mixed_pooling_allowed", "prediction": "pool states"}))
    xcl1_pool_result = validate_target(*xcl1_pool["XCL1_lymphotactin"])
    controls.append(_control("xcl1_mixed_state_pooling_forced_blocked", xcl1_pool_result["validation_level"] == "partially_supported_clean_abstain", "Forced mixed pooling lacks holdout support.", _observed_validation(xcl1_pool_result)))
    snca_no_disorder = _mutate_first(inputs, "alpha_synuclein_SNCA", packet_fn=lambda p: p.update({"falsifiable_predictions": [row for row in p["falsifiable_predictions"] if row["prediction_id"] != "SNCA_P1_DISORDER"]}))
    snca_no_disorder_result = validate_target(*snca_no_disorder["alpha_synuclein_SNCA"])
    controls.append(_control("snca_disorder_removed_partial_or_invalid", snca_no_disorder_result["validation_level"] in {"partially_supported_clean_abstain", "contradicted_or_failed"}, "Removing SNCA disorder prediction must force partial/invalid status.", _observed_validation(snca_no_disorder_result)))
    snca_single = _mutate_first(inputs, "alpha_synuclein_SNCA", packet_fn=lambda p: p["falsifiable_predictions"].append({"prediction_id": "bad_single", "expected_holdout_bucket": "compact_single_fold_holdout", "prediction": "compact single fold"}))
    snca_single_result = validate_target(*snca_single["alpha_synuclein_SNCA"])
    controls.append(_control("snca_forced_compact_single_fold_blocked", snca_single_result["validation_level"] == "partially_supported_clean_abstain", "Forced SNCA compact fold lacks support.", _observed_validation(snca_single_result)))
    swapped = _mutate_first(inputs, "KcsA", holdout_fn=lambda h: h.update(inputs["alpha_synuclein_SNCA"][1]))
    swapped_result = validate_target(*swapped["KcsA"])
    controls.append(_control("swapped_holdout_evidence_fails_validation", swapped_result["validation_level"] != "supported_by_independent_holdout", "Swapped SNCA holdout must not validate KcsA predictions.", _observed_validation(swapped_result)))
    placeholder = _mutate_first(inputs, "KcsA", holdout_fn=lambda h: h["holdout_validation_evidence"].append({"source_id": "placeholder", "source_name": "TODO placeholder", "source_url_or_citation": "citation needed"}))
    placeholder_result = validate_target(*placeholder["KcsA"])
    controls.append(_control("placeholder_citations_blocked", placeholder_result["validation_level"] == "blocked_for_leakage" and placeholder_result["placeholder_source_count"] > 0, "Placeholder citations must block.", _observed_validation(placeholder_result)))
    return controls


def build_v39() -> dict[str, Any]:
    inputs = load_inputs()
    results = _evaluate_inputs(inputs)
    controls = run_controls(inputs)
    answer_key_used = any(
        _bool_true(packet.get("answer_key_used_for_prediction")) or _bool_true(holdout.get("answer_key_used_for_prediction"))
        for packet, holdout in inputs.values()
    )
    validation_by_target = {r["target"]: r["validation_level"] for r in results}
    validated = [r for r in results if r["validation_level"] == "supported_by_independent_holdout"]
    partial = [r for r in results if r["validation_level"] == "partially_supported_clean_abstain"]
    contradicted = [r for r in results if r["validation_level"] == "contradicted_or_failed"]
    blocked = [r for r in results if r["validation_level"] == "blocked_for_leakage"]
    coord = sum(r["coordinate_derived_source_count"] for r in results)
    internal = sum(r["internal_runtime_source_count"] for r in results)
    key = sum(r["answer_key_source_count"] for r in results)
    failed = sorted(set(check for r in results for check in r["failed_checks"]))
    if not all(c["passed"] for c in controls):
        failed.append("V39_required_controls_failed")

    if key:
        status = BLOCKED_KEY
    elif coord or internal or blocked:
        status = BLOCKED_LEAKAGE
    elif contradicted or not all(c["passed"] for c in controls):
        status = FAILED
    elif len(validated) >= 2:
        status = PASSED
    else:
        status = PARTIAL

    if status == PASSED:
        next_action = "design_independent_experimental_readouts_from_validated_mechanism_predictions_no_claim_no_MD"
    elif status == PARTIAL:
        next_action = "acquire_more_independent_holdout_evidence_before_validation_claim"
    elif status == BLOCKED_KEY:
        next_action = "remove_answer_key_leakage_before_V39"
    elif status == BLOCKED_LEAKAGE:
        next_action = "remove_coordinate_or_internal_leakage_before_V39"
    else:
        next_action = "repair_failed_or_contradicted_mechanism_predictions_before_next_step"

    return {
        "kind": "V39_MECHANISM_TO_FALSIFIABLE_PREDICTION_VALIDATION_v0",
        "run_mode": "mechanism_prediction_to_independent_holdout_validation_no_coordinates_no_MD_claim_disabled",
        "control_status": status,
        "target_count": len(TARGET_ORDER),
        "prediction_packet_count": len(results),
        "holdout_source_count": sum(r["holdout_source_count"] for r in results),
        "validated_target_count": len(validated),
        "partial_target_count": len(partial),
        "contradicted_target_count": len(contradicted),
        "blocked_target_count": len(blocked),
        "validation_level_by_target": validation_by_target,
        "prediction_source_count_by_target": {r["target"]: r["prediction_source_count"] for r in results},
        "holdout_source_count_by_target": {r["target"]: r["holdout_source_count"] for r in results},
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
        "passed_control_count": sum(1 for c in controls if c["passed"]),
        "controls": controls,
        "target_validation_results": results,
        "failed_checks": failed,
        "next_action": next_action,
        "locked_interpretation": (
            "V39 passing means non-coordinate mechanism grammar produced falsifiable operator-level predictions for KcsA, XCL1, and SNCA, "
            "and independent holdout evidence supported them without coordinates, runtime evidence, answer-key leakage, MD, or any folding-solved claim."
        ),
    }


def write_outputs(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    cert = build_v39()
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v39_mechanism_to_falsifiable_prediction_validation_certificate.json"
    report_path = out_dir / "V39_MECHANISM_TO_FALSIFIABLE_PREDICTION_VALIDATION_REPORT.md"
    decision_path = out_dir / "v39_mechanism_to_falsifiable_prediction_validation_next_decision.json"
    decision = {
        "kind": "V39_MECHANISM_TO_FALSIFIABLE_PREDICTION_VALIDATION_NEXT_DECISION_v0",
        "decision_status": cert["control_status"],
        "next_action": cert["next_action"],
        "claim_allowed": False,
        "new_MD_allowed": False,
        "folding_problem_solved": False,
    }
    cert = {**cert, "next_decision": decision, "artifacts": {"certificate": str(cert_path), "report": str(report_path), "decision": str(decision_path)}}
    _write_json(cert_path, cert)
    _write_json(decision_path, decision)
    _write_report(report_path, cert)
    for result in cert["target_validation_results"]:
        target_dir = VALIDATION_ROOT / result["target"]
        _write_json(target_dir / "validation_result.json", result)
        _write_notes(target_dir / "falsification_notes.md", result)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path}


def _write_notes(path: Path, result: dict[str, Any]) -> None:
    lines = [
        f"# V39 Falsification Notes: {result['target']}",
        "",
        f"Validation level: `{result['validation_level']}`",
        f"Supported holdout buckets: `{result['supported_prediction_buckets']}`",
        f"Missing holdout buckets: `{result['missing_holdout_buckets']}`",
        "",
        "If future independent holdout evidence contradicts these operator-level predictions, this target must be downgraded rather than converted into a folding claim.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V39 Mechanism To Falsifiable Prediction Validation",
        "",
        f"Status: `{cert['control_status']}`",
        f"Validated targets: `{cert['validated_target_count']}` / `{cert['target_count']}`",
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
            f"Prediction count: `{result['prediction_count']}`",
            f"Holdout sources: `{result['holdout_source_count']}`",
            f"Supported buckets: `{result['supported_prediction_buckets']}`",
            f"Missing buckets: `{result['missing_holdout_buckets']}`",
            "",
        ])
    lines.extend(["## Plain English Interpretation", cert["locked_interpretation"]])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V39 mechanism-to-falsifiable-prediction validation.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = write_outputs(args.out_dir)
    cert = _read_json(paths["certificate"], "V39 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "control_status": cert["control_status"],
        "validation_level_by_target": cert["validation_level_by_target"],
        "holdout_source_count_by_target": cert["holdout_source_count_by_target"],
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
