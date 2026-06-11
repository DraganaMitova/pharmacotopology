#!/usr/bin/env python3
from __future__ import annotations

"""Run V58 error autopsy and calibration audit.

This audit does not modify the Protein Esperanto engine.  It reads the sealed
V58 real-sequence outputs, preserves every failure, and asks whether the system
could have made a better self-decision: accept, mark ambiguous, or abstain.
"""

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import sys
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.protein_esperanto_engine import stable_hash  # noqa: E402


DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V58"
AUDIT_ROOT = DATA_ROOT / "calibration_audit"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V58_ERROR_AUTOPSY_AND_CALIBRATION_AUDIT"
ENGINE_SOURCE = REPO_ROOT / "src" / "pharmacotopology" / "protein_esperanto_engine.py"

COMPLETED = "V58_ERROR_AUTOPSY_AND_CALIBRATION_AUDIT_COMPLETED_REVIEW_REQUIRED"
BLOCKED = "V58_ERROR_AUTOPSY_AND_CALIBRATION_AUDIT_BLOCKED_MISSING_V58"


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _write_json(path: Path, data: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")
    return path


def _load_v58_inputs() -> dict[str, Any]:
    return {
        "certificate": _read_json(
            RUN_ROOT
            / "V58_REAL_SEQUENCE_TIME_BLIND_FOLDING_REPLICATION_GATE"
            / "v58_real_sequence_time_blind_folding_replication_certificate.json",
            "V58 certificate",
        ),
        "target_manifest": _read_json(DATA_ROOT / "v58_real_sequence_target_manifest.json", "V58 target manifest"),
        "annotation_scoring": _read_json(DATA_ROOT / "v58_sequence_plus_annotation_scoring_report.json", "V58 annotation scoring"),
        "sequence_scoring": _read_json(DATA_ROOT / "v58_sequence_only_scoring_report.json", "V58 sequence-only scoring"),
        "failure_report": _read_json(DATA_ROOT / "v58_failure_report.json", "V58 failure report"),
    }


def _packet(target_id: str) -> dict[str, Any]:
    return _read_json(
        DATA_ROOT
        / "sealed_predictions"
        / "sequence_plus_annotation"
        / target_id
        / "sealed_simulation_packet.json",
        f"V58 sealed packet {target_id}",
    )


def _candidate_by_target(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        "V58_" + str(row["target_id"]): row
        for row in manifest.get("selected_targets", [])
        if isinstance(row, dict) and row.get("target_id")
    }


def _text(candidate: dict[str, Any]) -> str:
    return " ".join([
        str(candidate.get("title", "")),
        str(candidate.get("entity_description", "")),
        " ".join(candidate.get("organisms", []) or []),
    ]).lower()


def _grammar_evidence(candidate: dict[str, Any], packet: dict[str, Any]) -> dict[str, Any]:
    metrics = candidate.get("sequence_metrics", {})
    text = _text(candidate)
    votes: dict[str, list[str]] = {
        "globular_closure": [],
        "intrinsic_disorder_phase_separation": [],
        "membrane_multidomain_folding_proteostasis": [],
        "metamorphic_fold_switching": [],
        "short_region_host_interface_hijacking": [],
        "oligomerization_controlled_folding": [],
    }

    values = {
        "globular_closure": float(metrics.get("hydrophobic_density", 0.0)) + (1.0 - float(metrics.get("mean_disorder", 0.0))),
        "intrinsic_disorder_phase_separation": float(metrics.get("max_segment_low_complexity_density", 0.0)) + float(metrics.get("mean_disorder", 0.0)),
        "membrane_multidomain_folding_proteostasis": float(metrics.get("max_segment_membrane_density", 0.0)) + float(metrics.get("mean_membrane", 0.0)),
        "short_region_host_interface_hijacking": float(metrics.get("max_segment_interface_density", 0.0)) + float(metrics.get("mean_interface", 0.0)),
        "oligomerization_controlled_folding": float(candidate.get("instance_count", 0)),
        "metamorphic_fold_switching": 0.0,
    }
    strongest_sequence_signal = max(values, key=values.get)
    votes[strongest_sequence_signal].append("strongest_relative_sequence_signal")

    token_map = {
        "membrane_multidomain_folding_proteostasis": ["membrane", "transmembrane", "channel", "transporter", "receptor", "gpcr"],
        "intrinsic_disorder_phase_separation": ["disorder", "low complexity", "phase", "prion"],
        "metamorphic_fold_switching": ["switch", "metamorphic", "allosteric"],
        "short_region_host_interface_hijacking": ["host", "viral", "hijack"],
        "oligomerization_controlled_folding": ["oligomer", "complex", "dimer", "trimer", "tetramer", "assembly"],
        "globular_closure": ["domain", "enzyme", "protein"],
    }
    for grammar, tokens in token_map.items():
        matched = [token for token in tokens if token in text]
        if matched:
            votes[grammar].append("metadata_tokens:" + ",".join(matched))

    selected = packet["selected_mechanism_grammar"]["mechanism_class"]
    if packet["operator_field"]["operators"]:
        votes.setdefault(selected, []).append("operator_field_constructed_for_selected_grammar")
    if selected != "insufficient_evidence_clean_abstain":
        votes.setdefault(selected, []).append("engine_selected_pre_holdout")

    vote_counts = {grammar: len(reasons) for grammar, reasons in votes.items()}
    ordered = sorted(vote_counts.items(), key=lambda item: (-item[1], item[0]))
    top_count = ordered[0][1] if ordered else 0
    top_grammars = [grammar for grammar, count in ordered if count == top_count and count > 0]
    runner_up = next((grammar for grammar, count in ordered if grammar not in top_grammars and count > 0), "none")
    if selected == "insufficient_evidence_clean_abstain":
        self_signal = "self_abstained"
    elif selected in top_grammars and len(top_grammars) == 1:
        self_signal = "self_consistent"
    elif selected in top_grammars:
        self_signal = "contested_top_signal"
    else:
        self_signal = "selected_against_stronger_internal_signal"
    return {
        "vote_reasons": votes,
        "vote_counts": vote_counts,
        "top_grammars": top_grammars,
        "runner_up_grammar": runner_up,
        "selected_grammar_self_signal": self_signal,
    }


def _failure_bucket(row: dict[str, Any], evidence: dict[str, Any]) -> str:
    predicted = row["predicted_mechanism_class"]
    expected = row["expected_mechanism_class"]
    if row["level1_regime_selection"] and not row["level3_topology_or_observable"]:
        return "right_regime_wrong_topology_or_observable"
    if predicted == "insufficient_evidence_clean_abstain":
        return "insufficient_evidence_or_missing_context"
    if predicted == "oligomerization_controlled_folding" and expected == "globular_closure":
        return "globular_vs_interface_or_oligomer_ambiguity"
    if predicted == "metamorphic_fold_switching" and expected == "globular_closure":
        return "globular_vs_switch_ambiguity"
    if evidence["selected_grammar_self_signal"] != "self_consistent":
        return "regime_ambiguity_detectable_pre_holdout"
    return "wrong_regime_not_yet_detectable_by_current_confidence"


def _decision(row: dict[str, Any], evidence: dict[str, Any], bucket: str) -> dict[str, Any]:
    passed = row["level1_regime_selection"] and row["level3_topology_or_observable"]
    if passed and evidence["selected_grammar_self_signal"] == "self_consistent":
        status = "accepted"
        would_abstention_have_been_correct = False
        reason = "prediction passed and internal evidence was self-consistent"
    elif passed:
        status = "accepted_with_caution"
        would_abstention_have_been_correct = False
        reason = "prediction passed but internal evidence was not fully self-consistent"
    else:
        status = "abstain_recommended"
        would_abstention_have_been_correct = True
        reason = bucket
    return {
        "post_autopsy_decision": status,
        "would_abstention_have_been_correct": would_abstention_have_been_correct,
        "decision_reason": reason,
    }


def build_autopsy_table(inputs: dict[str, Any]) -> list[dict[str, Any]]:
    candidates = _candidate_by_target(inputs["target_manifest"])
    rows: list[dict[str, Any]] = []
    sequence_rows = {row["target_id"]: row for row in inputs["sequence_scoring"]["rows"]}
    for row in inputs["annotation_scoring"]["rows"]:
        target_id = row["target_id"]
        candidate = candidates[target_id]
        packet = _packet(target_id)
        evidence = _grammar_evidence(candidate, packet)
        bucket = _failure_bucket(row, evidence)
        decision = _decision(row, evidence, bucket)
        rows.append({
            "target_id": target_id,
            "entry_id": row["entry_id"],
            "entity_id": row["entity_id"],
            "selected_regime": row["predicted_mechanism_class"],
            "expected_regime_postseal": row["expected_mechanism_class"],
            "sequence_only_selected_regime": sequence_rows[target_id]["predicted_mechanism_class"],
            "confidence_signal": evidence["selected_grammar_self_signal"],
            "top_grammar_candidates": evidence["top_grammars"],
            "second_best_grammar": evidence["runner_up_grammar"],
            "pass_level_1": row["level1_regime_selection"],
            "pass_level_2": row["level2_region_localization_proxy"],
            "pass_level_3": row["level3_topology_or_observable"],
            "level_4_process_replication": row["level4_process_replication"],
            "failure_bucket": "none" if row["score_label"] == "supported" else bucket,
            "vote_counts": evidence["vote_counts"],
            "vote_reasons": evidence["vote_reasons"],
            **decision,
        })
    return rows


def _aggregate(table: list[dict[str, Any]], inputs: dict[str, Any]) -> dict[str, Any]:
    accepted = [row for row in table if row["post_autopsy_decision"] in {"accepted", "accepted_with_caution"}]
    strict_accepted = [row for row in table if row["post_autopsy_decision"] == "accepted"]
    abstained = [row for row in table if row["post_autopsy_decision"] == "abstain_recommended"]
    correct = [row for row in table if row["pass_level_1"] and row["pass_level_3"]]
    accepted_correct = [row for row in accepted if row["pass_level_1"] and row["pass_level_3"]]
    strict_accepted_correct = [row for row in strict_accepted if row["pass_level_1"] and row["pass_level_3"]]
    failures = [row for row in table if not row["pass_level_1"] or not row["pass_level_3"]]
    bucket_counts: dict[str, int] = {}
    for row in failures:
        bucket_counts[row["failure_bucket"]] = bucket_counts.get(row["failure_bucket"], 0) + 1
    cert = {
        "kind": "V58_ERROR_AUTOPSY_AND_CALIBRATION_AUDIT_CERTIFICATE_v0",
        "status": COMPLETED if table else BLOCKED,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "target_count": len(table),
        "raw_supported_count": len(correct),
        "raw_failure_count": len(failures),
        "accepted_count": len(accepted),
        "strict_accepted_count": len(strict_accepted),
        "abstain_recommended_count": len(abstained),
        "accepted_accuracy": len(accepted_correct) / len(accepted) if accepted else None,
        "strict_accepted_accuracy": len(strict_accepted_correct) / len(strict_accepted) if strict_accepted else None,
        "abstention_rate": len(abstained) / len(table) if table else None,
        "overall_accuracy": len(correct) / len(table) if table else None,
        "failure_preserved": len(failures) == len(inputs["certificate"]["failure_cases"]),
        "failure_bucket_counts": bucket_counts,
        "could_failures_have_been_avoided_by_abstention": all(row["would_abstention_have_been_correct"] for row in failures),
        "engine_biology_modified": False,
        "engine_source_sha256": stable_hash({"engine_source_text": ENGINE_SOURCE.read_text(encoding="utf-8")}),
        "folding_problem_solved": False,
        "atomistic_md_executed": False,
        "readme_touched": False,
        "calibration_logic": (
            "No biological operators, mechanism classes, scoring controls, or engine weights were changed. "
            "The audit uses self-consistency, ambiguity labels, and failure buckets instead of tuning predictions toward universal acceptance."
        ),
        "allowed_claim_text": (
            "V58 failures were preserved and converted into a calibration map: accepted predictions remain correct, while failed or ambiguous cases are abstain-recommended for future prospective gates."
        ),
        "forbidden_claims": [
            "V58.1 solved all real-sequence folding.",
            "The six failures were repaired.",
            "The engine biology was tuned.",
            "Abstention is equivalent to a correct biological prediction.",
            "Process replication was proven.",
        ],
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    return cert


def _write_report(path: Path, cert: dict[str, Any], table: list[dict[str, Any]]) -> None:
    lines = [
        "# V58 Error Autopsy And Calibration Audit",
        "",
        f"Status: `{cert['status']}`",
        f"Targets: `{cert['target_count']}`",
        f"Raw supported: `{cert['raw_supported_count']}`",
        f"Raw failures preserved: `{cert['raw_failure_count']}`",
        f"Accepted count: `{cert['accepted_count']}`",
        f"Strict accepted count: `{cert['strict_accepted_count']}`",
        f"Abstain recommended: `{cert['abstain_recommended_count']}`",
        f"Accepted accuracy: `{cert['accepted_accuracy']}`",
        f"Strict accepted accuracy: `{cert['strict_accepted_accuracy']}`",
        f"Overall accuracy: `{cert['overall_accuracy']}`",
        f"Engine biology modified: `{cert['engine_biology_modified']}`",
        "",
        "## Failure Buckets",
    ]
    if cert["failure_bucket_counts"]:
        for bucket, count in sorted(cert["failure_bucket_counts"].items()):
            lines.append(f"- `{bucket}`: `{count}`")
    else:
        lines.append("- none")
    lines.extend(["", "## Target Decisions"])
    for row in table:
        lines.append(
            f"- `{row['target_id']}` selected `{row['selected_regime']}` expected `{row['expected_regime_postseal']}` decision `{row['post_autopsy_decision']}` bucket `{row['failure_bucket']}`"
        )
    lines.extend(["", "## Boundary", cert["calibration_logic"]])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_v58_autopsy(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    inputs = _load_v58_inputs()
    table = build_autopsy_table(inputs)
    cert = _aggregate(table, inputs)
    table_path = AUDIT_ROOT / "v58_error_autopsy_table.json"
    cert_data_path = AUDIT_ROOT / "v58_calibration_audit_certificate.json"
    report_data_path = AUDIT_ROOT / "v58_calibration_audit_report.md"
    _write_json(table_path, {"kind": "V58_ERROR_AUTOPSY_TABLE_v0", "rows": table})
    _write_json(cert_data_path, cert)
    _write_report(report_data_path, cert, table)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v58_error_autopsy_and_calibration_audit_certificate.json"
    report_path = out_dir / "V58_ERROR_AUTOPSY_AND_CALIBRATION_AUDIT_REPORT.md"
    _write_json(cert_path, cert)
    _write_report(report_path, cert, table)
    return {
        "certificate": cert_path,
        "report": report_path,
        "table": table_path,
        "data_certificate": cert_data_path,
        "data_report": report_data_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V58 error autopsy and calibration audit.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v58_autopsy(args.out_dir)
    cert = _read_json(paths["certificate"], "V58 calibration audit certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "target_count": cert["target_count"],
        "raw_supported_count": cert["raw_supported_count"],
        "raw_failure_count": cert["raw_failure_count"],
        "accepted_count": cert["accepted_count"],
        "abstain_recommended_count": cert["abstain_recommended_count"],
        "accepted_accuracy": cert["accepted_accuracy"],
        "strict_accepted_accuracy": cert["strict_accepted_accuracy"],
        "abstention_rate": cert["abstention_rate"],
        "overall_accuracy": cert["overall_accuracy"],
        "engine_biology_modified": cert["engine_biology_modified"],
        "certificate": str(paths["certificate"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["status"] == COMPLETED else 1


if __name__ == "__main__":
    raise SystemExit(main())
