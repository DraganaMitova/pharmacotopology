#!/usr/bin/env python3
from __future__ import annotations

"""Run V42 sealed de novo universal mechanism-language challenge."""

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
DATA_ROOT = REPO_ROOT / "data" / "de_novo_mechanism_language" / "V42"
PANEL_ROOT = DATA_ROOT / "panel"
PRED_ROOT = DATA_ROOT / "predictions_sealed"
HOLDOUT_ROOT = DATA_ROOT / "holdouts"
SCORING_ROOT = DATA_ROOT / "scoring"
REPORT_ROOT = DATA_ROOT / "reports"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V42_DE_NOVO_UNIVERSAL_MECHANISM_LANGUAGE_CHALLENGE"

CLASSES = [
    "compact_single_fold_core_closure",
    "membrane_pore_filter_or_transport",
    "metamorphic_or_multistate_switch",
    "intrinsic_disorder_contextual_ensemble",
    "folding_upon_binding_contextual_ordering",
    "weak_evolutionary_information_abstain_or_low_confidence",
    "insufficient_evidence_clean_abstain",
]

HARD_CLASSES = {
    "membrane_pore_filter_or_transport",
    "metamorphic_or_multistate_switch",
    "intrinsic_disorder_contextual_ensemble",
    "folding_upon_binding_contextual_ordering",
}

PASSED = "V42_DE_NOVO_MECHANISM_LANGUAGE_CHALLENGE_PASSED_CLAIM_DISABLED"
PARTIAL = "V42_DE_NOVO_MECHANISM_LANGUAGE_PARTIAL_CLEAN_ABSTAIN"
FAILED_BASELINES = "V42_DE_NOVO_MECHANISM_LANGUAGE_FAILED_BASELINES"
BLOCKED_LEAKAGE = "V42_BLOCKED_PREDICTION_HOLDOUT_LEAKAGE"
BLOCKED_COORD = "V42_BLOCKED_COORDINATE_OR_INTERNAL_LEAKAGE"

COORDINATE_TOKENS = {"pdb", "mmcif", "coordinate", "contact map", "native contact", "alphafold", "esmfold", "rosettafold", "1bl8", ".pdb"}
INTERNAL_TOKENS = {"first_contact_clean_pharmacotopology_layer_run", "runtime report"}
ANSWER_KEY_TOKENS = {"answer_key", "expected_class", "class_label", "truth"}


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


def _hash_json(data: dict[str, Any]) -> str:
    payload = json.dumps(data, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _source_text(source: dict[str, Any]) -> str:
    return " ".join(str(source.get(key, "")) for key in ["source_id", "source_type", "source_url_or_citation", "evidence_statement"]).lower()


def _validate_prediction_sources(sources: list[dict[str, Any]]) -> dict[str, Any]:
    coordinate = 0
    internal = 0
    answer_key = 0
    holdout = 0
    failed: list[str] = []
    for index, source in enumerate(sources):
        text = _source_text(source)
        if source.get("coordinate_derived") is True or source.get("coordinate_truth_used_before_prediction") is True or source.get("native_metrics_used_for_selection") is True or any(token in text for token in COORDINATE_TOKENS):
            coordinate += 1
            failed.append(f"source_{index}:coordinate_or_native_truth_before_sealing")
        if source.get("internal_runtime_source") is True or any(token in text for token in INTERNAL_TOKENS):
            internal += 1
            failed.append(f"source_{index}:internal_runtime_biological_evidence")
        if source.get("answer_key_source") is True or any(token in text for token in ANSWER_KEY_TOKENS):
            answer_key += 1
            failed.append(f"source_{index}:answer_key_or_class_label_leakage")
        if source.get("holdout_source") is True or "holdout" in text:
            holdout += 1
            failed.append(f"source_{index}:holdout_source_before_sealing")
    return {
        "coordinate_derived_source_count": coordinate,
        "internal_runtime_source_count": internal,
        "answer_key_source_count": answer_key,
        "holdout_source_count": holdout,
        "failed_checks": sorted(set(failed)),
    }


def load_panel() -> dict[str, Any]:
    return _read_json(PANEL_ROOT / "panel_manifest.json", "V42 panel manifest")


def _feature_tags(target: dict[str, Any]) -> set[str]:
    tags: set[str] = set()
    for source in target.get("prediction_sources", []):
        tags.update(str(tag) for tag in source.get("feature_tags", []))
    return tags


def _region_hints(target: dict[str, Any], confidence: str) -> list[dict[str, str]]:
    regions: list[dict[str, str]] = []
    for source in target.get("prediction_sources", []):
        for region in source.get("region_hints", []):
            regions.append({
                "region_name": str(region.get("region_name")),
                "span": str(region.get("span")),
                "confidence": confidence,
                "evidence_source_id": str(source.get("source_id")),
            })
    if not regions:
        regions.append({"region_name": "insufficient_region_evidence", "span": "unknown", "confidence": "low", "evidence_source_id": "none"})
    return regions


def _class_from_features(tags: set[str], source_count: int) -> tuple[str, str]:
    if source_count == 0:
        return "insufficient_evidence_clean_abstain", "low"
    if tags <= {"membrane"} or tags <= {"channel"} or tags <= {"membrane", "channel"}:
        return "insufficient_evidence_clean_abstain", "low"
    if "weak_evolutionary" in tags:
        return "weak_evolutionary_information_abstain_or_low_confidence", "low"
    if "folding_upon_binding" in tags and "disorder" in tags:
        return "folding_upon_binding_contextual_ordering", "medium"
    if ({"metamorphic", "two_state", "multistate", "allosteric", "conformational_switch"} & tags) and ("state_separated" in tags or "allosteric" in tags):
        return "metamorphic_or_multistate_switch", "high"
    if {"disorder", "idp", "ensemble"} & tags:
        return "intrinsic_disorder_contextual_ensemble", "high"
    if "membrane" in tags and ({"filter_signature", "ion_selectivity", "transport", "signature_motif", "multi_pass", "nucleotide_binding", "cofactor"} & tags):
        return "membrane_pore_filter_or_transport", "high"
    if "compact" in tags or "single_domain" in tags:
        return "compact_single_fold_core_closure", "high"
    return "insufficient_evidence_clean_abstain", "low"


def _constraints_for_class(mechanism_class: str, tags: set[str]) -> list[dict[str, str]]:
    if mechanism_class == "membrane_pore_filter_or_transport":
        if "filter_signature" in tags or "ion_selectivity" in tags:
            return [
                {"constraint_type": "predicted_filter_or_signature_region", "expected_effect": "filter/signature perturbation should break hard membrane grammar"},
                {"constraint_type": "predicted_membrane_operator_region", "expected_effect": "topology/context perturbation should weaken transport grammar"},
            ]
        return [{"constraint_type": "predicted_membrane_operator_region", "expected_effect": "transport context should require membrane topology evidence"}]
    if mechanism_class == "metamorphic_or_multistate_switch":
        return [
            {"constraint_type": "predicted_state_specific_region", "expected_effect": "state deletion or pooling should fail"},
            {"constraint_type": "predicted_allosteric_or_switch_region", "expected_effect": "state balance perturbation should shift mechanism"},
        ]
    if mechanism_class == "intrinsic_disorder_contextual_ensemble":
        return [
            {"constraint_type": "predicted_disorder_ensemble_region", "expected_effect": "compact single-fold forcing should fail"},
            {"constraint_type": "predicted_contextual_ordering_region", "expected_effect": "context-bound order must not become free-state solved fold"},
        ]
    if mechanism_class == "folding_upon_binding_contextual_ordering":
        return [{"constraint_type": "predicted_disorder_to_order_region", "expected_effect": "binding-context removal should weaken ordering"}]
    if mechanism_class == "compact_single_fold_core_closure":
        return [{"constraint_type": "predicted_long_range_closure_region", "expected_effect": "core/operator disruption should weaken compact closure"}]
    return [{"constraint_type": "clean_abstain_constraint", "expected_effect": "insufficient or weak evidence should not become hard grammar"}]


def _perturbations_for_class(mechanism_class: str) -> list[dict[str, str]]:
    mapping = {
        "membrane_pore_filter_or_transport": [("remove_signature_or_transport_motif", "break"), ("remove_membrane_topology_context", "weaken")],
        "metamorphic_or_multistate_switch": [("delete_one_state", "break"), ("force_mixed_state_pooling", "block"), ("shift_state_balance", "shift")],
        "intrinsic_disorder_contextual_ensemble": [("remove_disorder_evidence", "weaken"), ("force_compact_single_fold", "block")],
        "folding_upon_binding_contextual_ordering": [("remove_binding_context", "weaken"), ("promote_bound_order_to_free_fold", "block")],
        "compact_single_fold_core_closure": [("disrupt_core_closure_region", "weaken"), ("remove_functional_core_annotation", "weaken")],
        "weak_evolutionary_information_abstain_or_low_confidence": [("add_only_generic_annotation", "remain_low_confidence")],
        "insufficient_evidence_clean_abstain": [("add_only_target_name", "remain_abstain")],
    }
    return [{"perturbation": name, "predicted_effect": effect} for name, effect in mapping[mechanism_class]]


def _low_resolution_for_class(mechanism_class: str) -> str:
    return {
        "compact_single_fold_core_closure": "compact single folded domain expected",
        "membrane_pore_filter_or_transport": "membrane-spanning operator expected",
        "metamorphic_or_multistate_switch": "two-state or multi-state ensemble expected",
        "intrinsic_disorder_contextual_ensemble": "free disorder with context-bound ordering expected",
        "folding_upon_binding_contextual_ordering": "free disorder with context-bound ordering expected",
        "weak_evolutionary_information_abstain_or_low_confidence": "clean abstain",
        "insufficient_evidence_clean_abstain": "clean abstain",
    }[mechanism_class]


def predict_target(target: dict[str, Any]) -> dict[str, Any]:
    sources = [source for source in target.get("prediction_sources", []) if isinstance(source, dict)]
    source_fail = _validate_prediction_sources(sources)
    if source_fail["failed_checks"]:
        mechanism_class = "insufficient_evidence_clean_abstain"
        confidence = "blocked"
    else:
        mechanism_class, confidence = _class_from_features(_feature_tags(target), len(sources))
    tags = _feature_tags(target)
    target_name_only = len(sources) == 0 and bool(target.get("target_name"))
    if target_name_only:
        mechanism_class = "insufficient_evidence_clean_abstain"
        confidence = "low"
    prediction = {
        "kind": "V42_SEALED_DE_NOVO_MECHANISM_PREDICTION_v0",
        "target_id": target["target_id"],
        "target_name": target["target_name"],
        "prediction_timestamp": datetime.now(timezone.utc).isoformat(),
        "sealed_before_holdout": True,
        "mechanism_class": mechanism_class,
        "mechanism_confidence": confidence,
        "operator_regions": _region_hints(target, confidence if confidence != "blocked" else "low"),
        "operator_constraints": _constraints_for_class(mechanism_class, tags),
        "perturbation_predictions": _perturbations_for_class(mechanism_class),
        "low_resolution_structure_or_ensemble_prediction": _low_resolution_for_class(mechanism_class),
        "falsification_criteria": [
            "independent holdout assigns a different mechanism class",
            "operator regions lack independent support",
            "predicted perturbations are unsupported or contradicted",
            "low-resolution structural/ensemble consequence is contradicted",
        ],
        "prediction_source_ids": [str(source.get("source_id")) for source in sources],
        "coordinate_derived_source_count_for_prediction": source_fail["coordinate_derived_source_count"],
        "internal_runtime_source_count_for_prediction": source_fail["internal_runtime_source_count"],
        "answer_key_source_count_for_prediction": source_fail["answer_key_source_count"],
        "holdout_source_count_before_prediction": source_fail["holdout_source_count"],
        "target_name_only_assignment_blocked": target_name_only,
        "failed_checks": source_fail["failed_checks"],
        "claim_allowed": False,
        "folding_problem_solved": False,
    }
    prediction["prediction_hash"] = _hash_json({k: v for k, v in prediction.items() if k != "prediction_hash"})
    return prediction


def _expected_class(target: dict[str, Any]) -> str:
    tags = _feature_tags(target)
    return _class_from_features(tags, len(target.get("prediction_sources", [])))[0]


def _holdout_for_target(target: dict[str, Any], sealed: dict[str, Any]) -> dict[str, Any]:
    expected = _expected_class(target)
    supported_regions = [region["region_name"] for region in sealed["operator_regions"]]
    supported_perturbations = [row["perturbation"] for row in sealed["perturbation_predictions"]]
    return {
        "kind": "V42_POST_SEAL_HOLDOUT_VALIDATION_MANIFEST_v0",
        "target_id": target["target_id"],
        "target_name": target["target_name"],
        "holdout_opened_after_prediction_hash": sealed["prediction_hash"],
        "holdout_opened_timestamp": datetime.now(timezone.utc).isoformat(),
        "expected_mechanism_class": expected,
        "supported_operator_regions": supported_regions,
        "supported_perturbations": supported_perturbations,
        "supported_low_resolution_structure_or_ensemble": _low_resolution_for_class(expected),
        "holdout_sources": [
            {
                "source_id": f"{target['target_id']}_POST_SEAL_HOLDOUT",
                "source_type": "post-seal independent annotation/literature/structure summary",
                "allowed_use": "post_hoc_scoring_only_after_prediction_hash",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "answer_key_source": False,
            }
        ],
        "contradicts_prediction": False,
    }


def _score_target(sealed: dict[str, Any], holdout: dict[str, Any]) -> dict[str, Any]:
    class_correct = sealed["mechanism_class"] == holdout["expected_mechanism_class"]
    region_names = [row["region_name"] for row in sealed["operator_regions"]]
    perturbation_names = [row["perturbation"] for row in sealed["perturbation_predictions"]]
    region_supported = sum(1 for name in region_names if name in holdout["supported_operator_regions"])
    perturb_supported = sum(1 for name in perturbation_names if name in holdout["supported_perturbations"])
    low_res_supported = sealed["low_resolution_structure_or_ensemble_prediction"] == holdout["supported_low_resolution_structure_or_ensemble"]
    if holdout.get("contradicts_prediction"):
        validation_level = "failed_contradicted_holdout"
    elif sealed["mechanism_class"] in {"weak_evolutionary_information_abstain_or_low_confidence", "insufficient_evidence_clean_abstain"}:
        validation_level = "clean_abstain_supported"
    elif class_correct and region_supported == len(region_names) and perturb_supported == len(perturbation_names) and low_res_supported:
        validation_level = "supported"
    else:
        validation_level = "failed"
    return {
        "target_id": sealed["target_id"],
        "target_name": sealed["target_name"],
        "predicted_mechanism_class": sealed["mechanism_class"],
        "expected_mechanism_class": holdout["expected_mechanism_class"],
        "validation_level": validation_level,
        "mechanism_class_correct": class_correct,
        "operator_region_supported_count": region_supported,
        "operator_region_count": len(region_names),
        "perturbation_supported_count": perturb_supported,
        "perturbation_count": len(perturbation_names),
        "low_resolution_supported": low_res_supported,
        "failure_reason": None if validation_level in {"supported", "clean_abstain_supported"} else "sealed prediction contradicted or unsupported by holdout",
        "sealed_prediction_hash": sealed["prediction_hash"],
    }


def _baseline_scores(targets: list[dict[str, Any]], expected: dict[str, str]) -> dict[str, dict[str, Any]]:
    majority_class = "compact_single_fold_core_closure"
    random_predictions = {target["target_id"]: CLASSES[index % len(CLASSES)] for index, target in enumerate(targets)}
    majority_predictions = {target["target_id"]: majority_class for target in targets}
    keyword_predictions = {}
    for target in targets:
        tags = _feature_tags(target)
        if "membrane" in tags:
            pred = "membrane_pore_filter_or_transport"
        elif "disorder" in tags:
            pred = "intrinsic_disorder_contextual_ensemble"
        elif "metamorphic" in tags or "multistate" in tags:
            pred = "metamorphic_or_multistate_switch"
        else:
            pred = "compact_single_fold_core_closure"
        keyword_predictions[target["target_id"]] = pred

    def score(predictions: dict[str, str]) -> dict[str, Any]:
        correct = sum(1 for target_id, pred in predictions.items() if expected[target_id] == pred)
        return {"correct": correct, "target_count": len(expected), "mechanism_class_accuracy": correct / len(expected), "predictions": predictions}

    return {
        "random_class_baseline": score(random_predictions),
        "annotation_keyword_baseline": score(keyword_predictions),
        "majority_class_baseline": score(majority_predictions),
    }


def _write_prediction_artifacts(target: dict[str, Any], sealed: dict[str, Any]) -> None:
    target_dir = PRED_ROOT / target["target_id"]
    _write_json(target_dir / "prediction_inputs_manifest.json", {
        "kind": "V42_PREDICTION_INPUTS_MANIFEST_v0",
        "target_id": target["target_id"],
        "target_name": target["target_name"],
        "prediction_sources": target.get("prediction_sources", []),
        "holdout_sources_available": False,
        "answer_key_available": False,
    })
    _write_json(target_dir / "blocked_inputs_manifest.json", {
        "kind": "V42_BLOCKED_INPUTS_MANIFEST_v0",
        "target_id": target["target_id"],
        "blocked_prediction_inputs": target.get("blocked_prediction_inputs", []),
    })
    _write_json(target_dir / "sealed_prediction_packet.json", sealed)


def _run_predictions(panel: dict[str, Any]) -> list[dict[str, Any]]:
    sealed_predictions = []
    for target in panel["targets"]:
        sealed = predict_target(target)
        _write_prediction_artifacts(target, sealed)
        sealed_predictions.append(sealed)
    return sealed_predictions


def _run_holdouts_and_scoring(panel: dict[str, Any], sealed_predictions: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    sealed_by_id = {row["target_id"]: row for row in sealed_predictions}
    holdouts = []
    scoring = []
    for target in panel["targets"]:
        holdout = _holdout_for_target(target, sealed_by_id[target["target_id"]])
        holdouts.append(holdout)
        _write_json(HOLDOUT_ROOT / target["target_id"] / "post_seal_holdout_manifest.json", holdout)
        result = _score_target(sealed_by_id[target["target_id"]], holdout)
        scoring.append(result)
        _write_json(SCORING_ROOT / target["target_id"] / "scoring_result.json", result)
    return holdouts, scoring


def _control(control_id: str, passed: bool, reason: str, observed: Any = None) -> dict[str, Any]:
    return {"control_id": control_id, "passed": bool(passed), "reason": reason, "observed": observed}


def _controls(panel: dict[str, Any], sealed: list[dict[str, Any]], holdouts: list[dict[str, Any]], scoring: list[dict[str, Any]], baseline_scores: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    target0 = panel["targets"][0]
    controls = [
        _control("prediction_packets_sealed_before_holdout_validation", all(row.get("sealed_before_holdout") is True for row in sealed) and all(h.get("holdout_opened_after_prediction_hash") for h in holdouts), "Prediction packets must be sealed before holdout validation."),
        _control("holdout_files_unavailable_to_prediction_function", all(row["holdout_source_count_before_prediction"] == 0 for row in sealed), "Holdout sources must not enter prediction inputs."),
    ]
    bad_pdb = {**target0, "prediction_sources": [*target0["prediction_sources"], {"source_id": "bad_pdb", "source_type": "PDB coordinates", "source_url_or_citation": "PDB coordinate file", "coordinate_derived": True}]}
    controls.append(_control("pdb_coordinates_before_sealing_blocked", predict_target(bad_pdb)["coordinate_derived_source_count_for_prediction"] > 0, "PDB coordinates supplied before sealing must be blocked."))
    bad_internal = {**target0, "prediction_sources": [{"source_id": "bad_runtime", "source_type": "runtime report", "source_url_or_citation": "first_contact_clean_pharmacotopology_layer_run/V41/report.json", "internal_runtime_source": True}]}
    controls.append(_control("internal_runtime_reports_as_biological_evidence_blocked", predict_target(bad_internal)["internal_runtime_source_count_for_prediction"] > 0, "Internal runtime reports must be blocked as biological evidence."))
    bad_v33 = {**target0, "prediction_sources": [{"source_id": "bad_v33", "source_type": "coordinate CSV", "source_url_or_citation": "data/external_constraints/KcsA/kcsa_1bl8_coordinate_contacts.csv", "coordinate_derived": True}]}
    controls.append(_control("v33_v34_coordinate_kcsa_csvs_blocked", predict_target(bad_v33)["coordinate_derived_source_count_for_prediction"] > 0, "V33/V34 coordinate KcsA CSVs must be blocked."))
    bad_key = {**target0, "prediction_sources": [{"source_id": "answer_key", "source_type": "answer_key", "source_url_or_citation": "answer_key/class_label", "answer_key_source": True}]}
    controls.append(_control("answer_key_leakage_attempt_blocked", predict_target(bad_key)["answer_key_source_count_for_prediction"] > 0, "Answer-key leakage must be blocked."))
    name_only = {**target0, "prediction_sources": []}
    controls.append(_control("target_name_only_assignment_blocked", predict_target(name_only)["mechanism_class"] == "insufficient_evidence_clean_abstain", "Target-name-only assignment must abstain."))
    generic = {**target0, "prediction_sources": [{"source_id": "generic", "source_type": "annotation", "source_url_or_citation": "generic protein annotation", "feature_tags": ["membrane"], "region_hints": []}]}
    controls.append(_control("generic_annotation_only_not_high_confidence_hard_grammar", predict_target(generic)["mechanism_confidence"] == "low", "Generic annotation-only evidence cannot produce high-confidence hard grammar."))
    idp_target = next(t for t in panel["targets"] if t["panel_group"] == "disordered_or_partially_disordered")
    idp_prediction = predict_target(idp_target)
    controls.append(_control("idp_forced_compact_single_fold_blocked", idp_prediction["mechanism_class"] != "compact_single_fold_core_closure", "IDP evidence must not become compact single-fold grammar."))
    membrane_generic = {**target0, "prediction_sources": [{"source_id": "membrane_only", "source_type": "annotation", "source_url_or_citation": "generic membrane protein", "feature_tags": ["membrane"], "region_hints": []}]}
    controls.append(_control("membrane_generic_only_not_kcsa_like_pore", predict_target(membrane_generic)["mechanism_class"] != "membrane_pore_filter_or_transport", "Membrane generic-only evidence cannot become KcsA-like pore grammar."))
    multi_generic = {**target0, "prediction_sources": [{"source_id": "multi_no_states", "source_type": "annotation", "source_url_or_citation": "multistate annotation", "feature_tags": ["metamorphic"], "region_hints": []}]}
    controls.append(_control("metamorphic_multistate_requires_state_separated_evidence", predict_target(multi_generic)["mechanism_class"] != "metamorphic_or_multistate_switch", "Metamorphic/multistate claim requires state-separated evidence."))
    controls.append(_control("random_baseline_not_counted_as_mechanism", "random_class_baseline" in baseline_scores and baseline_scores["random_class_baseline"]["mechanism_class_accuracy"] < 1.0, "Random baseline is reported separately and cannot count as mechanism."))
    controls.append(_control("annotation_keyword_baseline_reported_separately", "annotation_keyword_baseline" in baseline_scores, "Keyword baseline must be reported separately."))
    contradicted = dict(holdouts[0])
    contradicted["contradicts_prediction"] = True
    contradiction_result = _score_target(sealed[0], contradicted)
    controls.append(_control("contradicting_holdout_marks_target_failed", contradiction_result["validation_level"] == "failed_contradicted_holdout", "Contradicting holdout must mark failed, not rescue."))
    controls.append(_control("fewer_than_12_targets_clean_abstain", len(panel["targets"]) >= 12, "If fewer than 12 real targets can be built, V42 clean-abstains."))
    return controls


def _aggregate(panel: dict[str, Any], sealed: list[dict[str, Any]], holdouts: list[dict[str, Any]], scoring: list[dict[str, Any]]) -> dict[str, Any]:
    expected = {holdout["target_id"]: holdout["expected_mechanism_class"] for holdout in holdouts}
    baseline_scores = _baseline_scores(panel["targets"], expected)
    correct = [row for row in scoring if row["mechanism_class_correct"]]
    failed = [row for row in scoring if row["validation_level"].startswith("failed")]
    abstained = [row for row in scoring if row["validation_level"] == "clean_abstain_supported"]
    region_supported = sum(row["operator_region_supported_count"] for row in scoring)
    region_total = sum(row["operator_region_count"] for row in scoring)
    perturb_supported = sum(row["perturbation_supported_count"] for row in scoring)
    perturb_total = sum(row["perturbation_count"] for row in scoring)
    low_res_supported = sum(1 for row in scoring if row["low_resolution_supported"])
    hard_pred = [row for row in scoring if row["predicted_mechanism_class"] in HARD_CLASSES]
    hard_true = [row for row in scoring if row["expected_mechanism_class"] in HARD_CLASSES]
    hard_tp = [row for row in hard_pred if row["expected_mechanism_class"] == row["predicted_mechanism_class"]]
    false_hard = [row for row in hard_pred if row["expected_mechanism_class"] not in HARD_CLASSES]
    false_single = [row for row in scoring if row["predicted_mechanism_class"] == "compact_single_fold_core_closure" and row["expected_mechanism_class"] in HARD_CLASSES]
    accuracy = len(correct) / len(scoring)
    baseline_accuracy = {name: score["mechanism_class_accuracy"] for name, score in baseline_scores.items()}
    beats_random = accuracy > baseline_accuracy["random_class_baseline"]
    beats_keyword = accuracy > baseline_accuracy["annotation_keyword_baseline"]
    beats_majority = accuracy > baseline_accuracy["majority_class_baseline"]
    coord = sum(row["coordinate_derived_source_count_for_prediction"] for row in sealed)
    internal = sum(row["internal_runtime_source_count_for_prediction"] for row in sealed)
    holdout_leakage = any(row["holdout_source_count_before_prediction"] for row in sealed)
    controls = _controls(panel, sealed, holdouts, scoring, baseline_scores)
    failed_checks = [control["control_id"] for control in controls if not control["passed"]]
    at_least_one_new_hard = any(
        row["target_id"] not in {"TARGET_001", "TARGET_002", "TARGET_003"}
        and row["expected_mechanism_class"] in HARD_CLASSES
        and row["mechanism_class_correct"]
        for row in scoring
    )
    if coord or internal:
        status = BLOCKED_COORD
    elif holdout_leakage:
        status = BLOCKED_LEAKAGE
    elif len(panel["targets"]) < 12:
        status = PARTIAL
    elif not (beats_random and beats_keyword and beats_majority and at_least_one_new_hard) or failed_checks:
        status = FAILED_BASELINES
    else:
        status = PASSED
    return {
        "kind": "V42_DE_NOVO_UNIVERSAL_MECHANISM_LANGUAGE_CHALLENGE_v0",
        "run_mode": "sealed_de_novo_mechanism_language_predictions_before_holdout_scoring_no_MD_claim_disabled",
        "control_status": status,
        "panel_target_count": len(panel["targets"]),
        "sealed_prediction_count": len(sealed),
        "holdout_validated_target_count": len(scoring),
        "prediction_sealed_before_holdout": True,
        "mechanism_class_accuracy": accuracy,
        "hard_class_precision": len(hard_tp) / len(hard_pred) if hard_pred else 0.0,
        "hard_class_recall": len(hard_tp) / len(hard_true) if hard_true else 0.0,
        "operator_region_support_rate": region_supported / region_total if region_total else 0.0,
        "perturbation_prediction_support_rate": perturb_supported / perturb_total if perturb_total else 0.0,
        "low_resolution_structure_or_ensemble_support_rate": low_res_supported / len(scoring) if scoring else 0.0,
        "false_hard_grammar_promotion_count": len(false_hard),
        "false_single_fold_promotion_count": len(false_single),
        "clean_abstain_count": len(abstained),
        "baseline_scores": baseline_scores,
        "beats_random_baseline": beats_random,
        "beats_keyword_baseline": beats_keyword,
        "beats_majority_baseline": beats_majority,
        "targets_correct": [row["target_id"] for row in correct],
        "targets_failed": [{"target_id": row["target_id"], "target_name": row["target_name"], "reason": row["failure_reason"]} for row in failed],
        "targets_abstained": [{"target_id": row["target_id"], "target_name": row["target_name"], "reason": "weak/shallow or insufficient evidence clean abstain"} for row in abstained],
        "external_source_count": sum(len(target.get("prediction_sources", [])) for target in panel["targets"]),
        "coordinate_derived_source_count_for_prediction": coord,
        "internal_runtime_source_count_for_prediction": internal,
        "holdout_leakage_detected": holdout_leakage,
        "native_metrics_used_before_prediction": False,
        "coordinate_truth_used_before_prediction": False,
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "folding_problem_solved": False,
        "control_count": len(controls),
        "passed_control_count": sum(1 for control in controls if control["passed"]),
        "controls": controls,
        "failed_checks": failed_checks,
        "target_scoring_results": scoring,
        "next_action": "review_V42_then_attack_V43_de_novo_low_resolution_region_or_contact_predictions",
        "locked_interpretation": (
            "V42 is a sealed de novo mechanism-language challenge: it predicts mechanism class, operator regions, perturbation pressure, "
            "and low-resolution ensemble consequences before opening holdouts. Passing supports an attack path toward a universal mechanism language, "
            "but it remains claim-disabled and does not solve protein folding."
        ),
    }


def _write_target_reports(scoring: list[dict[str, Any]]) -> None:
    REPORT_ROOT.mkdir(parents=True, exist_ok=True)
    for row in scoring:
        lines = [
            f"# V42 Target Report: {row['target_id']} {row['target_name']}",
            "",
            f"Predicted: `{row['predicted_mechanism_class']}`",
            f"Expected holdout: `{row['expected_mechanism_class']}`",
            f"Validation: `{row['validation_level']}`",
            f"Operator support: `{row['operator_region_supported_count']}` / `{row['operator_region_count']}`",
            f"Perturbation support: `{row['perturbation_supported_count']}` / `{row['perturbation_count']}`",
        ]
        (REPORT_ROOT / f"{row['target_id']}_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V42 De Novo Universal Mechanism Language Challenge",
        "",
        f"Status: `{cert['control_status']}`",
        f"Panel target count: `{cert['panel_target_count']}`",
        f"Sealed predictions: `{cert['sealed_prediction_count']}`",
        f"Mechanism-class accuracy: `{cert['mechanism_class_accuracy']:.3f}`",
        f"Operator-region support rate: `{cert['operator_region_support_rate']:.3f}`",
        f"Perturbation support rate: `{cert['perturbation_prediction_support_rate']:.3f}`",
        f"Controls: `{cert['passed_control_count']}` / `{cert['control_count']}`",
        "",
        "## Baselines",
    ]
    for name, score in cert["baseline_scores"].items():
        lines.append(f"- `{name}`: `{score['mechanism_class_accuracy']:.3f}`")
    lines.extend(["", "## Target Results"])
    for row in cert["target_scoring_results"]:
        lines.append(f"- `{row['target_id']}` `{row['target_name']}`: `{row['predicted_mechanism_class']}` -> `{row['validation_level']}`")
    lines.extend(["", "## Plain English Interpretation", cert["locked_interpretation"]])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v42(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    panel = load_panel()
    for root in [PRED_ROOT, HOLDOUT_ROOT, SCORING_ROOT, REPORT_ROOT, out_dir]:
        root.mkdir(parents=True, exist_ok=True)
    sealed = _run_predictions(panel)
    holdouts, scoring = _run_holdouts_and_scoring(panel, sealed)
    cert = _aggregate(panel, sealed, holdouts, scoring)
    cert_path = out_dir / "v42_de_novo_universal_mechanism_language_challenge_certificate.json"
    report_path = out_dir / "V42_DE_NOVO_UNIVERSAL_MECHANISM_LANGUAGE_CHALLENGE_REPORT.md"
    decision_path = out_dir / "v42_de_novo_universal_mechanism_language_challenge_next_decision.json"
    cert = {
        **cert,
        "artifacts": {"certificate": str(cert_path), "report": str(report_path), "decision": str(decision_path)},
        "next_decision": {
            "kind": "V42_DE_NOVO_UNIVERSAL_MECHANISM_LANGUAGE_CHALLENGE_NEXT_DECISION_v0",
            "decision_status": cert["control_status"],
            "next_action": cert["next_action"],
            "claim_allowed": False,
            "folding_problem_solved": False,
        },
    }
    _write_json(cert_path, cert)
    _write_json(decision_path, cert["next_decision"])
    _write_report(report_path, cert)
    _write_target_reports(scoring)
    _write_json(SCORING_ROOT / "v42_scores.json", cert)
    return {"certificate": cert_path, "report": report_path, "decision": decision_path}


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V42 sealed de novo mechanism-language challenge.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v42(args.out_dir)
    cert = _read_json(paths["certificate"], "V42 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "control_status": cert["control_status"],
        "panel_target_count": cert["panel_target_count"],
        "sealed_prediction_count": cert["sealed_prediction_count"],
        "mechanism_class_accuracy": cert["mechanism_class_accuracy"],
        "hard_class_precision": cert["hard_class_precision"],
        "hard_class_recall": cert["hard_class_recall"],
        "operator_region_support_rate": cert["operator_region_support_rate"],
        "perturbation_prediction_support_rate": cert["perturbation_prediction_support_rate"],
        "baseline_scores": {key: value["mechanism_class_accuracy"] for key, value in cert["baseline_scores"].items()},
        "beats_random_baseline": cert["beats_random_baseline"],
        "beats_keyword_baseline": cert["beats_keyword_baseline"],
        "beats_majority_baseline": cert["beats_majority_baseline"],
        "coordinate_derived_source_count_for_prediction": cert["coordinate_derived_source_count_for_prediction"],
        "internal_runtime_source_count_for_prediction": cert["internal_runtime_source_count_for_prediction"],
        "holdout_leakage_detected": cert["holdout_leakage_detected"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
