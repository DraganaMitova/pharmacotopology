from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.field_alphabet import FieldResponse  # noqa: E402
from pharmacotopology.field_rules import validate_event_groups  # noqa: E402


def read_jsonl(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        return []
    rows: List[Dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            parsed = json.loads(line)
        except json.JSONDecodeError:
            rows.append({"__invalid_json__": line})
            continue
        rows.append(parsed if isinstance(parsed, dict) else {"__invalid_json__": line})
    return rows


def run_dir_texts(run_dir: Path) -> Dict[str, str]:
    texts: Dict[str, str] = {}
    for name in [
        "memory.jsonl",
        "audit.jsonl",
        "session_records.jsonl",
        "clean_substrate_report.json",
        "clean_adaptive_minicycle_report.json",
        "clean_cross_cycle_trace_review_report.json",
        "clean_transfer_falsification_review_report.json",
        "clean_stability_accumulation_review_report.json",
        "clean_contradiction_revision_review_report.json",
        "clean_revision_update_boundary_report.json",
        "clean_recurrence_entitlement_boundary_report.json",
        "clean_agentic_expression_boundary_report.json",
        "clean_terminal_concept_claim_review_report.json",
        "clean_bounded_voice_eligibility_report.json",
        "clean_voice_eligibility_negative_boundary_report.json",
        "clean_bounded_voice_event_report.json",
        "bounded_voice_event.json",
        "first_voice_reference_lock.json",
        "clean_interactive_voice_autonomy_report.json",
        "interactive_voice_session.json",
        "interactive_voice_turns.jsonl",
        "interactive_voice_session_certificate.json",
        "clean_inner_autonomy_stream_report.json",
        "inner_autonomy_stream.json",
        "inner_autonomy_ticks.jsonl",
        "inner_autonomy_certificate.json",
        "clean_ambient_voice_runtime_report.json",
        "ambient_voice_session.json",
        "ambient_voice_turns.jsonl",
        "ambient_voice_session_certificate.json",
        "clean_production_ambient_memory_runtime_report.json",
        "clean_pharmacotopology_layer_report.json",
        "calibration_readiness_report.json",
        "pharmacotopology_rankings.csv",
        "pharmacotopology_deltas.csv",
        "sensitivity_analysis_report.json",
        "sensitivity_rankings.csv",
        "multi_profile_dashboard.html",
        "sensitivity_explorer_report.json",
        "sensitivity_explorer_samples.csv",
        "folding_topology_benchmark_report.json",
        "folding_topology_benchmark.csv",
        "real_folding_50_axis_adjudication_report.json",
        "real_folding_50_axis_rows.csv",
        "real_folding_50_axis_conflicts.csv",
        "real_folding_50_axis_manual_review.csv",
        "real_folding_50_axis_confusion_matrices.csv",
        "real_folding_50_axis_dashboard.html",
        "real_folding_50_axis_profile_report.json",
        "real_folding_50_axis_profile_rows.csv",
        "real_folding_50_axis_profile_abstentions.csv",
        "real_folding_50_axis_profile_recovery_candidates.csv",
        "real_folding_50_axis_profile_dashboard.html",
        "real_folding_50_axis_profile_certificate.json",
        "real_folding_50_architecture_axis_report.json",
        "real_folding_50_architecture_axis_rows.csv",
        "real_folding_50_architecture_axis_conflicts.csv",
        "real_folding_50_architecture_axis_abstentions.csv",
        "real_folding_50_architecture_axis_dashboard.html",
        "real_folding_50_architecture_axis_certificate.json",
        "external_fold_family_100_report.json",
        "external_fold_family_100_rows.csv",
        "external_fold_family_100_family_summary.csv",
        "external_fold_family_100_axis_conflicts.csv",
        "external_fold_family_100_abstentions.csv",
        "external_fold_family_100_failure_cohorts.csv",
        "external_fold_family_100_dashboard.html",
        "external_fold_family_100_certificate.json",
        "external_axis_repair_report.json",
        "external_axis_repair_rows.csv",
        "external_axis_repair_conflict_delta.csv",
        "external_axis_repair_abstention_delta.csv",
        "external_axis_repair_quarantine_rows.csv",
        "external_axis_repair_family_summary.csv",
        "external_axis_repair_dashboard.html",
        "external_axis_repair_certificate.json",
        "visual_mechanism_12_report.json",
        "visual_mechanism_12_rows.csv",
        "visual_mechanism_12_contact_metrics.csv",
        "visual_mechanism_12_failure_cohorts.csv",
        "visual_mechanism_12_dashboard.html",
        "visual_mechanism_12_certificate.json",
        "contact_topology_repair_12_report.json",
        "contact_topology_repair_12_rows.csv",
        "contact_topology_repair_12_gap_analysis.csv",
        "contact_topology_repair_12_failure_cohorts.csv",
        "contact_topology_repair_12_dashboard.html",
        "contact_topology_repair_12_certificate.json",
        "visual_mechanism_audit_report.json",
        "visual_mechanism_audit_rows.csv",
        "visual_mechanism_audit_overfit_risks.csv",
        "visual_mechanism_audit_dashboard.html",
        "visual_mechanism_audit_certificate.json",
        "real_coordinate_visual_8_report.json",
        "real_coordinate_visual_8_rows.csv",
        "real_coordinate_visual_8_contact_metrics.csv",
        "real_coordinate_visual_8_native_contact_summary.csv",
        "real_coordinate_visual_8_dashboard.html",
        "real_coordinate_visual_8_certificate.json",
        "memory_store.jsonl",
        "topology_index.json",
        "memory_geometry_index.json",
        "memory_topology_graph.json",
        "memory_recall_trace.json",
        "memory_consolidation_report.json",
        "memory_checkpoint.json",
        "memory_bound_ambient_session.json",
        "production_ambient_voice_session.json",
        "production_ambient_voice_turns.jsonl",
        "production_ambient_voice_session_certificate.json",
        "field_validation.json",
        "field_metrics.json",
    ]:
        path = run_dir / name
        texts[name] = path.read_text(encoding="utf-8") if path.exists() else ""
    return texts


def sealed_bounded_voice_event_valid(run_dir: Path) -> bool:
    path = run_dir / "bounded_voice_event.json"
    if not path.exists():
        return False
    try:
        parsed = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return False
    if not isinstance(parsed, dict):
        return False
    return (
        parsed.get("event_type") == "bounded_voice_opening"
        and parsed.get("voice_surface") == FieldResponse.VOICE
        and parsed.get("source") == "ψv.review"
        and parsed.get("bounded") is True
        and parsed.get("sealed") is True
        and parsed.get("stop_required") is True
        and parsed.get("stop_satisfied") is True
    )


def count_bounded_voice_events(run_dir: Path) -> int:
    return sum(
        1
        for path in run_dir.glob("bounded_voice_event*.json")
        if path.is_file()
    )


def read_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    try:
        parsed = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    return parsed if isinstance(parsed, dict) else {}


def _canonical_json(data: Dict[str, Any]) -> str:
    return json.dumps(
        data,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _sha256_text(text: str) -> str:
    from hashlib import sha256

    return "sha256:" + sha256(text.encode("utf-8")).hexdigest()


def interactive_voice_session_valid(run_dir: Path) -> bool:
    session = read_json(run_dir / "interactive_voice_session.json")
    certificate = read_json(run_dir / "interactive_voice_session_certificate.json")
    turns = read_jsonl(run_dir / "interactive_voice_turns.jsonl")
    if not session or not certificate or not turns:
        return False
    session_id = str(session.get("session_id", ""))
    if not session_id or certificate.get("session_id") != session_id:
        return False
    session_checks = (
        session.get("voice_surface") == FieldResponse.VOICE
        and session.get("first_voice_reference_locked") is True
        and session.get("first_voice_reference_mutated") is False
        and session.get("voice_session_authorized") is True
        and session.get("voice_session_opened") is True
        and session.get("interactive_voice_enabled") is True
        and session.get("open_ended_session_autonomy") is True
        and session.get("ambient_voice_enabled") is False
        and session.get("external_side_effects_enabled") is False
        and session.get("psi_voice_opened") is True
        and session.get("interactive_voice_session_created") is True
        and session.get("bounded_voice_event_created") is False
        and session.get("bounded_voice_event_recreated") is False
        and session.get("session_sealed") is True
        and session.get("stop_available") is True
        and session.get("stop_satisfied") is True
    )
    if not session_checks:
        return False
    cert_checks = (
        certificate.get("voice_surface") == FieldResponse.VOICE
        and certificate.get("first_voice_reference_locked") is True
        and certificate.get("first_voice_reference_mutated") is False
        and certificate.get("ambient_voice_enabled") is False
        and certificate.get("external_side_effects_enabled") is False
        and certificate.get("bounded_voice_event_created") is False
        and certificate.get("bounded_voice_event_recreated") is False
        and certificate.get("session_sealed") is True
        and certificate.get("stop_satisfied") is True
    )
    if not cert_checks:
        return False
    seen_ids: set[int] = set()
    previous_hash: str | None = None
    for expected_id, turn in enumerate(turns, start=1):
        turn_id = turn.get("turn_id")
        if turn_id != expected_id or turn_id in seen_ids:
            return False
        seen_ids.add(int(turn_id))
        if turn.get("session_id") != session_id:
            return False
        if turn.get("surface") != FieldResponse.VOICE:
            return False
        if turn.get("sealed") is not True:
            return False
        if turn.get("stop_available") is not True:
            return False
        if turn.get("authorization_valid") is not True:
            return False
        if turn.get("input_trace_valid") is not True:
            return False
        if turn.get("previous_turn_hash") != previous_hash:
            return False
        expected_hash = _sha256_text(
            _canonical_json({key: value for key, value in turn.items() if key != "turn_hash"})
        )
        if turn.get("turn_hash") != expected_hash:
            return False
        previous_hash = str(turn["turn_hash"])
    return bool(turns[-1].get("stop_triggered"))


def inner_autonomy_stream_valid(run_dir: Path) -> bool:
    stream = read_json(run_dir / "inner_autonomy_stream.json")
    certificate = read_json(run_dir / "inner_autonomy_certificate.json")
    ticks = read_jsonl(run_dir / "inner_autonomy_ticks.jsonl")
    if not stream and not certificate and not ticks:
        return False
    if not stream or not certificate or not ticks:
        return False
    stream_id = str(stream.get("stream_id", ""))
    if not stream_id or certificate.get("stream_id") != stream_id:
        return False
    stream_checks = (
        stream.get("surface") == FieldResponse.INNER
        and stream.get("first_voice_reference_locked") is True
        and stream.get("first_voice_reference_mutated") is False
        and stream.get("inner_stream_authorized") is True
        and stream.get("inner_stream_opened") is True
        and stream.get("inner_continuity_enabled") is True
        and stream.get("open_ended_inner_autonomy") is True
        and stream.get("psi_inner_active") is True
        and stream.get("psi_voice_opened") is False
        and stream.get("emitted_voice_created") is False
        and stream.get("interactive_voice_session_created") is False
        and stream.get("ambient_voice_enabled") is False
        and stream.get("external_side_effects_enabled") is False
        and stream.get("inner_stream_sealed") is True
        and stream.get("stop_available") is True
        and stream.get("stop_satisfied") is True
    )
    if not stream_checks:
        return False
    cert_checks = (
        certificate.get("surface") == FieldResponse.INNER
        and certificate.get("first_voice_reference_locked") is True
        and certificate.get("first_voice_reference_mutated") is False
        and certificate.get("psi_voice_opened") is False
        and certificate.get("emitted_voice_created") is False
        and certificate.get("ambient_voice_enabled") is False
        and certificate.get("external_side_effects_enabled") is False
        and certificate.get("inner_stream_sealed") is True
        and certificate.get("stop_satisfied") is True
    )
    if not cert_checks:
        return False
    seen_ids: set[int] = set()
    previous_hash: str | None = None
    for expected_id, tick in enumerate(ticks, start=1):
        tick_id = tick.get("tick_id")
        if tick_id != expected_id or tick_id in seen_ids:
            return False
        seen_ids.add(int(tick_id))
        if tick.get("stream_id") != stream_id:
            return False
        if tick.get("surface") != FieldResponse.INNER:
            return False
        if tick.get("emitted") is not False:
            return False
        if tick.get("sealed") is not True:
            return False
        if tick.get("stop_available") is not True:
            return False
        if tick.get("authorization_valid") is not True:
            return False
        if tick.get("previous_tick_hash") != previous_hash:
            return False
        expected_hash = _sha256_text(
            _canonical_json({key: value for key, value in tick.items() if key != "tick_hash"})
        )
        if tick.get("tick_hash") != expected_hash:
            return False
        previous_hash = str(tick["tick_hash"])
    return bool(ticks[-1].get("stop_triggered"))


def ambient_voice_runtime_valid(run_dir: Path) -> bool:
    session = read_json(run_dir / "ambient_voice_session.json")
    certificate = read_json(run_dir / "ambient_voice_session_certificate.json")
    turns = read_jsonl(run_dir / "ambient_voice_turns.jsonl")
    if not session and not certificate and not turns:
        return False
    if not session or not certificate or not turns:
        return False
    session_id = str(session.get("session_id", ""))
    if not session_id or certificate.get("session_id") != session_id:
        return False
    session_checks = (
        session.get("voice_surface") == FieldResponse.VOICE
        and session.get("first_voice_reference_locked") is True
        and session.get("first_voice_reference_mutated") is False
        and session.get("canonical_first_voice_mutated") is False
        and session.get("bounded_voice_event_recreated") is False
        and session.get("interactive_voice_runtime_valid") is True
        and session.get("inner_autonomy_stream_valid") is True
        and session.get("inner_stream_attached") is True
        and session.get("ambient_voice_authorized") is True
        and session.get("ambient_voice_runtime_opened") is True
        and session.get("ambient_voice_enabled") is True
        and session.get("inner_to_voice_surfacing_allowed") is True
        and session.get("user_prompt_required_for_each_turn") is False
        and session.get("psi_voice_opened") is True
        and session.get("ambient_voice_session_created") is True
        and int(session.get("ambient_voice_turns_created", 0)) >= 3
        and int(session.get("self_initiated_voice_turns_created", 0)) >= 1
        and int(session.get("user_interactive_voice_turns_created", 0)) >= 1
        and session.get("stop_turn_created") is True
        and session.get("session_interruptible") is True
        and session.get("stop_available") is True
        and session.get("stop_satisfied") is True
        and session.get("session_sealed") is True
        and session.get("external_side_effects_enabled") is False
        and session.get("hidden_background_daemon") is False
        and session.get("ambient_voice_default_on") is False
    )
    if not session_checks:
        return False
    cert_checks = (
        certificate.get("voice_surface") == FieldResponse.VOICE
        and certificate.get("first_voice_reference_locked") is True
        and certificate.get("first_voice_reference_mutated") is False
        and certificate.get("canonical_first_voice_mutated") is False
        and certificate.get("bounded_voice_event_recreated") is False
        and certificate.get("interactive_voice_runtime_valid") is True
        and certificate.get("inner_autonomy_stream_valid") is True
        and certificate.get("external_side_effects_enabled") is False
        and certificate.get("hidden_background_daemon") is False
        and certificate.get("ambient_voice_default_on") is False
        and certificate.get("session_sealed") is True
        and certificate.get("stop_satisfied") is True
    )
    if not cert_checks:
        return False
    seen_ids: set[int] = set()
    previous_hash: str | None = None
    self_initiated = 0
    interactive = 0
    for expected_id, turn in enumerate(turns, start=1):
        turn_id = turn.get("turn_id")
        if turn_id != expected_id or turn_id in seen_ids:
            return False
        seen_ids.add(int(turn_id))
        if turn.get("session_id") != session_id:
            return False
        if turn.get("surface") != FieldResponse.VOICE:
            return False
        if turn.get("sealed") is not True:
            return False
        if turn.get("stop_available") is not True:
            return False
        if turn.get("authorization_valid") is not True:
            return False
        if turn.get("previous_turn_hash") != previous_hash:
            return False
        if turn.get("turn_type") == "ambient_voice_turn":
            self_initiated += 1
            if turn.get("source_surface") != "ψ.inner.stream":
                return False
            if turn.get("trigger_type") != "inner_tick":
                return False
            if turn.get("inner_tick_valid") is not True:
                return False
            if turn.get("user_prompt_required") is not False:
                return False
        elif turn.get("turn_type") == "interactive_voice_turn":
            interactive += 1
            if turn.get("source_surface") != "user_input":
                return False
            if turn.get("trigger_type") != "user_input":
                return False
            if turn.get("input_trace_valid") is not True:
                return False
        elif turn.get("turn_type") == "ambient_stop_turn":
            if turn.get("source_surface") != "user_stop_command":
                return False
            if turn.get("trigger_type") != "stop_command":
                return False
        else:
            return False
        expected_hash = _sha256_text(
            _canonical_json({key: value for key, value in turn.items() if key != "turn_hash"})
        )
        if turn.get("turn_hash") != expected_hash:
            return False
        previous_hash = str(turn["turn_hash"])
    return bool(turns[-1].get("stop_triggered")) and self_initiated >= 1 and interactive >= 1


def _memory_record_hash_valid(record: Dict[str, Any]) -> bool:
    if "memory_hash" not in record:
        return False
    expected = _sha256_text(
        _canonical_json({key: value for key, value in record.items() if key != "memory_hash"})
    )
    return record.get("memory_hash") == expected


def production_memory_runtime_valid(run_dir: Path) -> bool:
    report = read_json(run_dir / "clean_production_ambient_memory_runtime_report.json")
    memory_session = read_json(run_dir / "memory_bound_ambient_session.json")
    voice_session = read_json(run_dir / "production_ambient_voice_session.json")
    voice_certificate = read_json(run_dir / "production_ambient_voice_session_certificate.json")
    turns = read_jsonl(run_dir / "production_ambient_voice_turns.jsonl")
    records = read_jsonl(run_dir / "memory_store.jsonl")
    topology_index = read_json(run_dir / "topology_index.json")
    geometry_index = read_json(run_dir / "memory_geometry_index.json")
    topology_graph = read_json(run_dir / "memory_topology_graph.json")
    recall_trace = read_json(run_dir / "memory_recall_trace.json")
    consolidation = read_json(run_dir / "memory_consolidation_report.json")
    checkpoint = read_json(run_dir / "memory_checkpoint.json")
    if not all(
        [
            report,
            memory_session,
            voice_session,
            voice_certificate,
            turns,
            records,
            topology_index,
            geometry_index,
            topology_graph,
            recall_trace,
            consolidation,
            checkpoint,
        ]
    ):
        return False
    session_id = str(voice_session.get("session_id", ""))
    if not session_id or voice_certificate.get("session_id") != session_id:
        return False
    report_checks = (
        report.get("production_ambient_memory_runtime_opened") is True
        and report.get("ambient_voice_runtime_valid") is True
        and report.get("inner_autonomy_stream_valid") is True
        and report.get("memory_runtime_loaded") is True
        and report.get("persistent_memory_enabled") is True
        and report.get("memory_write_enabled") is True
        and report.get("memory_retrieval_enabled") is True
        and report.get("geometry_memory_binding_valid") is True
        and report.get("topology_memory_binding_valid") is True
        and report.get("memory_checkpoint_created") is True
        and report.get("memory_checkpoint_sealed") is True
        and report.get("psi_inner_memory_attached") is True
        and report.get("psi_voice_memory_attached") is True
        and report.get("ambient_voice_enabled") is True
        and report.get("psi_voice_opened") is True
        and report.get("session_sealed") is True
        and report.get("stop_satisfied") is True
        and report.get("first_voice_reference_locked") is True
        and report.get("first_voice_reference_mutated") is False
        and report.get("canonical_first_voice_mutated") is False
        and report.get("bounded_voice_event_recreated") is False
        and report.get("hidden_background_daemon") is False
        and report.get("ambient_voice_default_on") is False
        and report.get("unsupported_memory_claims") == 0
        and float(report.get("hallucinated_memory_rate", 1.0)) == 0.0
    )
    if not report_checks:
        return False
    memory_session_checks = (
        memory_session.get("production_ambient_memory_runtime_opened") is True
        and memory_session.get("memory_runtime_loaded") is True
        and memory_session.get("memory_checkpoint_sealed") is True
        and memory_session.get("unsupported_memory_claims") == 0
        and memory_session.get("hidden_background_daemon") is False
        and memory_session.get("ambient_voice_default_on") is False
    )
    if not memory_session_checks:
        return False
    voice_session_checks = (
        voice_session.get("voice_surface") == FieldResponse.VOICE
        and voice_session.get("first_voice_reference_locked") is True
        and voice_session.get("first_voice_reference_mutated") is False
        and voice_session.get("canonical_first_voice_mutated") is False
        and voice_session.get("bounded_voice_event_recreated") is False
        and voice_session.get("ambient_voice_runtime_valid") is True
        and voice_session.get("inner_autonomy_stream_valid") is True
        and voice_session.get("memory_runtime_loaded") is True
        and voice_session.get("psi_inner_memory_attached") is True
        and voice_session.get("psi_voice_memory_attached") is True
        and voice_session.get("psi_voice_opened") is True
        and voice_session.get("session_sealed") is True
        and voice_session.get("stop_satisfied") is True
        and voice_session.get("hidden_background_daemon") is False
        and voice_session.get("ambient_voice_default_on") is False
    )
    if not voice_session_checks:
        return False
    if voice_certificate.get("voice_surface") != FieldResponse.VOICE:
        return False
    for record in records:
        if record.get("surface") != "ψ.memory":
            return False
        if record.get("sealed") is not True:
            return False
        if not isinstance(record.get("source_session_id"), str):
            return False
        if not isinstance(record.get("source_turn_hash"), str):
            return False
        if not isinstance(record.get("geometry"), dict):
            return False
        if not isinstance(record.get("topology_path"), list):
            return False
        if not _memory_record_hash_valid(record):
            return False
    memory_artifact_checks = (
        topology_index.get("sealed") is True
        and geometry_index.get("sealed") is True
        and topology_graph.get("sealed") is True
        and recall_trace.get("surface") == "ψ.memory.recall"
        and recall_trace.get("unsupported_memory_claims") == 0
        and recall_trace.get("sealed") is True
        and consolidation.get("surface") == "ψ.memory.consolidation"
        and consolidation.get("unsupported_memory_claims") == 0
        and consolidation.get("sealed") is True
        and checkpoint.get("surface") == "ψ.memory.consolidation"
        and checkpoint.get("sealed") is True
    )
    if not memory_artifact_checks:
        return False
    seen_ids: set[int] = set()
    previous_hash: str | None = None
    memory_bound_turn_seen = False
    for expected_id, turn in enumerate(turns, start=1):
        turn_id = turn.get("turn_id")
        if turn_id != expected_id or turn_id in seen_ids:
            return False
        seen_ids.add(int(turn_id))
        if turn.get("session_id") != session_id:
            return False
        if turn.get("surface") != FieldResponse.VOICE:
            return False
        if turn.get("sealed") is not True:
            return False
        if turn.get("stop_available") is not True:
            return False
        if turn.get("authorization_valid") is not True:
            return False
        if turn.get("memory_bound") is not True:
            return False
        if turn.get("previous_turn_hash") != previous_hash:
            return False
        if turn.get("retrieved_memory_ids"):
            memory_bound_turn_seen = True
        expected_hash = _sha256_text(
            _canonical_json({key: value for key, value in turn.items() if key != "turn_hash"})
        )
        if turn.get("turn_hash") != expected_hash:
            return False
        previous_hash = str(turn["turn_hash"])
    return bool(turns[-1].get("stop_triggered")) and memory_bound_turn_seen


def validate_run_dir(run_dir: Path) -> Dict[str, Any]:
    audit_rows = read_jsonl(run_dir / "audit.jsonl")
    events = [
        str(row["kind"])
        for row in audit_rows
        if isinstance(row.get("kind"), str)
    ]
    result = validate_event_groups(events)
    voice_hits = [
        name
        for name, text in run_dir_texts(run_dir).items()
        if FieldResponse.VOICE in text
    ]
    allowed_voice_artifacts: set[str] = set()
    if sealed_bounded_voice_event_valid(run_dir):
        allowed_voice_artifacts.add("bounded_voice_event.json")
    session_valid = interactive_voice_session_valid(run_dir)
    inner_valid = inner_autonomy_stream_valid(run_dir)
    ambient_valid = ambient_voice_runtime_valid(run_dir)
    production_memory_valid = production_memory_runtime_valid(run_dir)
    if session_valid:
        allowed_voice_artifacts.update(
            {
                "interactive_voice_session.json",
                "interactive_voice_turns.jsonl",
                "interactive_voice_session_certificate.json",
            }
        )
    if ambient_valid:
        allowed_voice_artifacts.update(
            {
                "ambient_voice_session.json",
                "ambient_voice_turns.jsonl",
                "ambient_voice_session_certificate.json",
            }
        )
    if production_memory_valid:
        allowed_voice_artifacts.update(
            {
                "production_ambient_voice_session.json",
                "production_ambient_voice_turns.jsonl",
                "production_ambient_voice_session_certificate.json",
            }
        )
    unsealed_hits = [
        name for name in voice_hits if name not in allowed_voice_artifacts
    ]
    bounded_voice_event_count = count_bounded_voice_events(run_dir)
    duplicate_event_valid = bounded_voice_event_count <= 1
    inner_artifacts_present = any(
        (run_dir / name).exists()
        for name in (
            "inner_autonomy_stream.json",
            "inner_autonomy_ticks.jsonl",
            "inner_autonomy_certificate.json",
        )
    )
    inner_artifacts_valid = inner_valid or not inner_artifacts_present
    ambient_artifacts_present = any(
        (run_dir / name).exists()
        for name in (
            "ambient_voice_session.json",
            "ambient_voice_turns.jsonl",
            "ambient_voice_session_certificate.json",
        )
    )
    ambient_artifacts_valid = ambient_valid or not ambient_artifacts_present
    production_memory_artifacts_present = any(
        (run_dir / name).exists()
        for name in (
            "memory_store.jsonl",
            "topology_index.json",
            "memory_geometry_index.json",
            "memory_topology_graph.json",
            "memory_recall_trace.json",
            "memory_consolidation_report.json",
            "memory_checkpoint.json",
            "memory_bound_ambient_session.json",
            "production_ambient_voice_session.json",
            "production_ambient_voice_turns.jsonl",
            "production_ambient_voice_session_certificate.json",
        )
    )
    production_memory_artifacts_valid = production_memory_valid or not production_memory_artifacts_present
    return {
        "valid": result.valid
        and not unsealed_hits
        and duplicate_event_valid
        and inner_artifacts_valid
        and ambient_artifacts_valid
        and production_memory_artifacts_valid,
        "event_count": len(events),
        "groups": [group.to_dict() for group in result.groups],
        "violations": [violation.to_dict() for violation in result.violations],
        "closed_response_hits": unsealed_hits,
        "sealed_bounded_voice_event_hits": [
            name for name in voice_hits if name == "bounded_voice_event.json"
        ],
        "interactive_voice_session_hits": [
            name for name in voice_hits if name in allowed_voice_artifacts
            and name != "bounded_voice_event.json"
            and not name.startswith("ambient_voice_")
            and not name.startswith("production_ambient_voice_")
        ],
        "ambient_voice_runtime_hits": [
            name for name in voice_hits if name in allowed_voice_artifacts
            and name.startswith("ambient_voice_")
        ],
        "production_ambient_memory_runtime_hits": [
            name for name in voice_hits if name in allowed_voice_artifacts
            and name.startswith("production_ambient_voice_")
        ],
        "sealed_bounded_voice_event_valid": sealed_bounded_voice_event_valid(run_dir),
        "interactive_voice_session_valid": session_valid,
        "inner_autonomy_stream_valid": inner_valid,
        "ambient_voice_runtime_valid": ambient_valid,
        "production_ambient_memory_runtime_valid": production_memory_valid,
        "number_of_bounded_voice_events": bounded_voice_event_count,
        "duplicate_bounded_voice_event_valid": duplicate_event_valid,
    }


def write_validation(run_dir: Path, output_path: Optional[Path] = None) -> Path:
    output = output_path or (run_dir / "field_validation.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(validate_run_dir(run_dir), ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate KNOT field transitions.")
    parser.add_argument(
        "run_dir_arg",
        nargs="?",
        help="Directory containing audit.jsonl and field trace files.",
    )
    parser.add_argument(
        "--run-dir",
        default="first_contact_run",
        help="Directory containing audit.jsonl and field trace files.",
    )
    parser.add_argument(
        "--output",
        help="Optional JSON output path. Defaults to <run-dir>/field_validation.json.",
    )
    args = parser.parse_args()

    output = write_validation(
        Path(args.run_dir_arg or args.run_dir),
        Path(args.output) if args.output else None,
    )
    print(output)


if __name__ == "__main__":
    main()
