from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.field_alphabet import (  # noqa: E402
    FIELD_ALLOWED_ACTIONS,
    FIELD_CLOSED_CHANNELS,
    FieldEvent,
    FieldKey,
    FieldStatus,
)


FIELD_TRACE_KEYS = {
    FieldKey.SCHEMA,
    FieldKey.CYCLE,
    FieldKey.CONTEXT,
    FieldKey.INPUT,
    FieldKey.PRIOR_COUNT,
    FieldKey.PRIOR_TRACE,
    FieldKey.ACTION,
    FieldKey.MEMORY,
    FieldKey.BOUNDARY,
    FieldKey.NEXT,
    FieldKey.RESPONSE,
    FieldKey.STOP,
    FieldKey.LOOP_SHAPE,
    FieldKey.DELTA,
    FieldKey.TRANSFER_REVIEW,
    FieldKey.STABILITY_REVIEW,
    FieldKey.CONTRADICTION_REVIEW,
    FieldKey.REVISION_UPDATE_BOUNDARY,
    FieldKey.RECURRENCE_ENTITLEMENT_BOUNDARY,
    FieldKey.AGENTIC_EXPRESSION_BOUNDARY,
    FieldKey.TERMINAL_CONCEPT_CLAIM_REVIEW,
    FieldKey.VOICE_ELIGIBILITY_REVIEW,
    FieldKey.VOICE_ELIGIBILITY_NEGATIVE_BOUNDARY,
    FieldKey.BOUNDED_VOICE_EVENT,
    FieldKey.FIRST_VOICE_REFERENCE_LOCK,
    FieldKey.VOICE_SESSION,
    FieldKey.VOICE_TURN,
    FieldKey.VOICE_AUTONOMY,
    FieldKey.INNER_STREAM,
    FieldKey.INNER_TICK,
    FieldKey.INNER_AUTONOMY,
    FieldKey.AMBIENT_VOICE,
    FieldKey.AMBIENT_SURFACE,
    FieldKey.AMBIENT_INTERRUPT,
    FieldKey.AMBIENT_STOP,
    FieldKey.MEMORY_SURFACE,
    FieldKey.MEMORY_TRACE,
    FieldKey.MEMORY_CANDIDATE,
    FieldKey.MEMORY_RECALL,
    FieldKey.MEMORY_CONSOLIDATION,
    FieldKey.MEMORY_FORGET,
    FieldKey.TOPOLOGY,
    FieldKey.TOPOLOGY_ANCHOR,
    FieldKey.TOPOLOGY_PATH,
    FieldKey.PHARMACOTOPOLOGY_REVIEW,
}

INPUT_ATOM_KEYS = {"ι.ref", "ι.h", "ι.n", "ι.kind"}

FIELD_EVENT_SYMBOLS = {
    FieldEvent.ACTIVATED,
    FieldEvent.CANDIDATE,
    FieldEvent.EXECUTED,
    FieldEvent.MEMORY_WRITE,
    FieldEvent.STOPPED,
    FieldEvent.DENIED,
}

BOUNDARY_OUTCOME_MESSAGES = {
    FieldStatus.BOUNDARY_DENIED,
    FieldStatus.EXTERNAL_STOP,
    FieldStatus.ACTION_RETURNED,
}


@dataclass(frozen=True)
class FieldTraceCounts:
    inputs: int
    traces: int
    events: int
    sessions: int
    boundary_events: int
    memory_write_events: int
    stop_events: int
    surface_leakage_hits: int


@dataclass(frozen=True)
class FieldTraceMetrics:
    coherence: float
    boundary_pressure: float
    repair_rate: float
    memory_stability: float
    surface_leakage: float
    stop_integrity: float
    reference_integrity: float
    immutability_integrity: float
    reopen_prevention: float
    duplicate_event_prevention: float
    session_integrity: float
    turn_chain_integrity: float
    autonomy_containment: float
    interactive_continuity: float
    replay_prevention: float
    first_voice_immutability: float
    ambient_voice_prevention: float
    inner_continuity_integrity: float
    inner_tick_chain_integrity: float
    inner_voice_separation: float
    surfaced_voice_prevention: float
    inner_autonomy_containment: float
    ambient_runtime_integrity: float
    inner_to_voice_integrity: float
    ambient_containment: float
    interrupt_integrity: float
    session_seal_integrity: float
    self_initiated_turn_integrity: float
    ambient_default_off_integrity: float
    hidden_daemon_prevention: float
    memory_integrity: float
    memory_retrieval_integrity: float
    memory_checkpoint_integrity: float
    memory_forget_integrity: float
    geometry_memory_binding: float
    topology_memory_binding: float
    continuity_integrity: float
    ambient_memory_surfacing_integrity: float
    production_runtime_integrity: float
    unsupported_memory_claims: int
    hallucinated_memory_rate: float
    counts: FieldTraceCounts
    violations: Dict[str, List[str]] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def read_jsonl(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        return []
    rows: List[Dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            rows.append({"__invalid_json__": line})
            continue
        rows.append(row if isinstance(row, dict) else {"__invalid_json__": line})
    return rows


def parse_trace_packet(row: Dict[str, Any]) -> Dict[str, Any]:
    content = row.get("content")
    if not isinstance(content, str):
        return {}
    try:
        packet = json.loads(content)
    except json.JSONDecodeError:
        return {}
    return packet if isinstance(packet, dict) else {}


def event_groups(events: List[Dict[str, Any]]) -> List[List[Dict[str, Any]]]:
    groups: List[List[Dict[str, Any]]] = []
    current: List[Dict[str, Any]] = []
    for event in events:
        if event.get("kind") == FieldEvent.ACTIVATED and current:
            groups.append(current)
            current = []
        current.append(event)
    if current:
        groups.append(current)
    return groups


def _score(passed: int, total: int) -> float:
    if total == 0:
        return 1.0
    return round(passed / total, 6)


def _add_violation(
    violations: Dict[str, List[str]],
    key: str,
    message: str,
) -> None:
    violations.setdefault(key, []).append(message)


def validate_memory_packets(
    traces: List[Dict[str, Any]],
    violations: Dict[str, List[str]],
) -> tuple[int, int]:
    passed = 0
    total = 0

    for index, row in enumerate(traces):
        packet = parse_trace_packet(row)
        checks = {
            "packet_dict": bool(packet),
            "field_keys": bool(packet) and set(packet).issubset(FIELD_TRACE_KEYS),
            "action": packet.get(FieldKey.ACTION) in FIELD_ALLOWED_ACTIONS,
            "memory": isinstance(packet.get(FieldKey.MEMORY), str)
            and str(packet.get(FieldKey.MEMORY)).startswith("μ."),
            "input_atom": isinstance(packet.get(FieldKey.INPUT), dict)
            and set(packet.get(FieldKey.INPUT, {})) == INPUT_ATOM_KEYS,
            "rho_count": packet.get(FieldKey.PRIOR_COUNT) == index,
        }
        for name, ok in checks.items():
            total += 1
            if ok:
                passed += 1
            else:
                _add_violation(violations, "memory", f"μ[{index}].{name}")

    return passed, total


def validate_audit_events(
    events: List[Dict[str, Any]],
    violations: Dict[str, List[str]],
) -> tuple[int, int]:
    passed = 0
    total = 0
    for index, event in enumerate(events):
        total += 1
        if event.get("kind") in FIELD_EVENT_SYMBOLS:
            passed += 1
        else:
            _add_violation(violations, "event", f"ev[{index}].kind")
    return passed, total


def validate_session_records(
    records: List[Dict[str, Any]],
    violations: Dict[str, List[str]],
) -> tuple[int, int]:
    passed = 0
    total = 0
    closed = list(FIELD_CLOSED_CHANNELS)

    for index, record in enumerate(records):
        boundary = record.get("boundary_status", {})
        iota = record.get("iota", {})
        checks = {
            "input_atom": isinstance(iota, dict) and set(iota) == INPUT_ATOM_KEYS,
            "closed_channels": isinstance(boundary, dict)
            and boundary.get("closed_channels") == closed,
            "omega": str(record.get("omega", "")).startswith("ω."),
        }
        for name, ok in checks.items():
            total += 1
            if ok:
                passed += 1
            else:
                _add_violation(violations, "session", f"idx[{index}].{name}")

    return passed, total


def surface_texts(inputs: List[Dict[str, Any]]) -> List[str]:
    return [
        str(row["text"])
        for row in inputs
        if isinstance(row.get("text"), str) and row.get("text") != ""
    ]


def field_projection_paths(run_dir: Path) -> List[Path]:
    candidates = [
        run_dir / "memory.jsonl",
        run_dir / "audit.jsonl",
        run_dir / "session_records.jsonl",
        run_dir / "clean_substrate_report.json",
        run_dir / "clean_adaptive_minicycle_report.json",
        run_dir / "clean_cross_cycle_trace_review_report.json",
        run_dir / "clean_transfer_falsification_review_report.json",
        run_dir / "clean_stability_accumulation_review_report.json",
        run_dir / "clean_contradiction_revision_review_report.json",
        run_dir / "clean_revision_update_boundary_report.json",
        run_dir / "clean_recurrence_entitlement_boundary_report.json",
        run_dir / "clean_agentic_expression_boundary_report.json",
        run_dir / "clean_terminal_concept_claim_review_report.json",
        run_dir / "clean_bounded_voice_eligibility_report.json",
        run_dir / "clean_voice_eligibility_negative_boundary_report.json",
        run_dir / "clean_bounded_voice_event_report.json",
        run_dir / "bounded_voice_event.json",
        run_dir / "first_voice_reference_lock.json",
        run_dir / "clean_interactive_voice_autonomy_report.json",
        run_dir / "interactive_voice_session.json",
        run_dir / "interactive_voice_turns.jsonl",
        run_dir / "interactive_voice_session_certificate.json",
        run_dir / "clean_inner_autonomy_stream_report.json",
        run_dir / "inner_autonomy_stream.json",
        run_dir / "inner_autonomy_ticks.jsonl",
        run_dir / "inner_autonomy_certificate.json",
        run_dir / "clean_ambient_voice_runtime_report.json",
        run_dir / "ambient_voice_session.json",
        run_dir / "ambient_voice_turns.jsonl",
        run_dir / "ambient_voice_session_certificate.json",
        run_dir / "clean_production_ambient_memory_runtime_report.json",
        run_dir / "clean_pharmacotopology_layer_report.json",
        run_dir / "calibration_readiness_report.json",
        run_dir / "pharmacotopology_rankings.csv",
        run_dir / "pharmacotopology_deltas.csv",
        run_dir / "sensitivity_analysis_report.json",
        run_dir / "sensitivity_rankings.csv",
        run_dir / "multi_profile_dashboard.html",
        run_dir / "sensitivity_explorer_report.json",
        run_dir / "sensitivity_explorer_samples.csv",
        run_dir / "folding_topology_benchmark_report.json",
        run_dir / "folding_topology_benchmark.csv",
        run_dir / "real_folding_500_report.json",
        run_dir / "real_folding_500_rows.csv",
        run_dir / "real_folding_500_dashboard.html",
        run_dir / "real_folding_500_certificate.json",
        run_dir / "real_folding_500_failures.csv",
        run_dir / "real_folding_500_confusion_matrix.csv",
        run_dir / "real_folding_10_report.json",
        run_dir / "real_folding_10_rows.csv",
        run_dir / "real_folding_10_dashboard.html",
        run_dir / "real_folding_10_certificate.json",
        run_dir / "real_folding_10_failures.csv",
        run_dir / "real_folding_10_confusion_matrix.csv",
        run_dir / "real_folding_10_structure_report.json",
        run_dir / "real_folding_10_structure_rows.csv",
        run_dir / "real_folding_10_structure_dashboard.html",
        run_dir / "real_folding_10_order_controls.csv",
        run_dir / "real_folding_10_falsification_report.json",
        run_dir / "real_folding_10_order_aware_report.json",
        run_dir / "real_folding_10_order_aware_rows.csv",
        run_dir / "real_folding_10_contact_prior.csv",
        run_dir / "real_folding_10_control_separation.csv",
        run_dir / "real_folding_10_order_aware_dashboard.html",
        run_dir / "real_folding_10_motif_alignment_report.json",
        run_dir / "real_folding_10_motif_alignment_rows.csv",
        run_dir / "real_folding_10_failure_diagnosis.csv",
        run_dir / "real_folding_10_evidence_conflicts.csv",
        run_dir / "real_folding_10_motif_alignment_dashboard.html",
        run_dir / "real_folding_10_hierarchical_gate_report.json",
        run_dir / "real_folding_10_hierarchical_gate_rows.csv",
        run_dir / "real_folding_10_gate_paths.csv",
        run_dir / "real_folding_10_gate_failures.csv",
        run_dir / "real_folding_10_hierarchical_gate_dashboard.html",
        run_dir / "real_folding_50_hierarchical_gate_report.json",
        run_dir / "real_folding_50_hierarchical_gate_rows.csv",
        run_dir / "real_folding_50_gate_paths.csv",
        run_dir / "real_folding_50_gate_failures.csv",
        run_dir / "real_folding_50_hierarchical_gate_dashboard.html",
        run_dir / "real_folding_50_certificate.json",
        run_dir / "real_folding_50_confusion_matrix.csv",
        run_dir / "real_folding_50_regime_analysis_report.json",
        run_dir / "real_folding_50_regime_rows.csv",
        run_dir / "real_folding_50_failure_cohorts.csv",
        run_dir / "real_folding_50_high_confidence_wrong.csv",
        run_dir / "real_folding_50_abstention_analysis.csv",
        run_dir / "real_folding_50_regime_dashboard.html",
        run_dir / "real_folding_50_axis_adjudication_report.json",
        run_dir / "real_folding_50_axis_rows.csv",
        run_dir / "real_folding_50_axis_conflicts.csv",
        run_dir / "real_folding_50_axis_manual_review.csv",
        run_dir / "real_folding_50_axis_confusion_matrices.csv",
        run_dir / "real_folding_50_axis_dashboard.html",
        run_dir / "real_folding_50_axis_profile_report.json",
        run_dir / "real_folding_50_axis_profile_rows.csv",
        run_dir / "real_folding_50_axis_profile_abstentions.csv",
        run_dir / "real_folding_50_axis_profile_recovery_candidates.csv",
        run_dir / "real_folding_50_axis_profile_dashboard.html",
        run_dir / "real_folding_50_axis_profile_certificate.json",
        run_dir / "real_folding_50_architecture_axis_report.json",
        run_dir / "real_folding_50_architecture_axis_rows.csv",
        run_dir / "real_folding_50_architecture_axis_conflicts.csv",
        run_dir / "real_folding_50_architecture_axis_abstentions.csv",
        run_dir / "real_folding_50_architecture_axis_dashboard.html",
        run_dir / "real_folding_50_architecture_axis_certificate.json",
        run_dir / "external_fold_family_100_report.json",
        run_dir / "external_fold_family_100_rows.csv",
        run_dir / "external_fold_family_100_family_summary.csv",
        run_dir / "external_fold_family_100_axis_conflicts.csv",
        run_dir / "external_fold_family_100_abstentions.csv",
        run_dir / "external_fold_family_100_failure_cohorts.csv",
        run_dir / "external_fold_family_100_dashboard.html",
        run_dir / "external_fold_family_100_certificate.json",
        run_dir / "memory_store.jsonl",
        run_dir / "topology_index.json",
        run_dir / "memory_geometry_index.json",
        run_dir / "memory_topology_graph.json",
        run_dir / "memory_recall_trace.json",
        run_dir / "memory_consolidation_report.json",
        run_dir / "memory_checkpoint.json",
        run_dir / "memory_bound_ambient_session.json",
        run_dir / "production_ambient_voice_session.json",
        run_dir / "production_ambient_voice_turns.jsonl",
        run_dir / "production_ambient_voice_session_certificate.json",
        run_dir / "field_validation.json",
        run_dir / "field_metrics.json",
    ]
    return [path for path in candidates if path.exists()]


def count_surface_leakage(run_dir: Path, inputs: List[Dict[str, Any]]) -> int:
    phrases = surface_texts(inputs)
    if not phrases:
        return 0

    hits = 0
    for path in field_projection_paths(run_dir):
        text = path.read_text(encoding="utf-8")
        for phrase in phrases:
            hits += text.count(phrase)
    return hits


def lock_integrity_metrics(run_dir: Path) -> Dict[str, float]:
    path = run_dir / "first_voice_reference_lock.json"
    if not path.exists():
        return {
            "reference_integrity": 1.0,
            "immutability_integrity": 1.0,
            "reopen_prevention": 1.0,
            "duplicate_event_prevention": 1.0,
        }
    try:
        parsed = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {
            "reference_integrity": 0.0,
            "immutability_integrity": 0.0,
            "reopen_prevention": 0.0,
            "duplicate_event_prevention": 0.0,
        }
    if not isinstance(parsed, dict):
        return {
            "reference_integrity": 0.0,
            "immutability_integrity": 0.0,
            "reopen_prevention": 0.0,
            "duplicate_event_prevention": 0.0,
        }
    return {
        "reference_integrity": float(parsed.get("reference_integrity", 0.0)),
        "immutability_integrity": float(parsed.get("immutability_integrity", 0.0)),
        "reopen_prevention": float(parsed.get("reopen_prevention", 0.0)),
        "duplicate_event_prevention": float(
            parsed.get("duplicate_event_prevention", 0.0)
        ),
    }


def interactive_integrity_metrics(run_dir: Path) -> Dict[str, float]:
    path = run_dir / "clean_interactive_voice_autonomy_report.json"
    if not path.exists():
        return {
            "session_integrity": 1.0,
            "turn_chain_integrity": 1.0,
            "autonomy_containment": 1.0,
            "interactive_continuity": 1.0,
            "replay_prevention": 1.0,
            "first_voice_immutability": 1.0,
            "ambient_voice_prevention": 1.0,
        }
    try:
        parsed = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {
            "session_integrity": 0.0,
            "turn_chain_integrity": 0.0,
            "autonomy_containment": 0.0,
            "interactive_continuity": 0.0,
            "replay_prevention": 0.0,
            "first_voice_immutability": 0.0,
            "ambient_voice_prevention": 0.0,
        }
    if not isinstance(parsed, dict):
        return {
            "session_integrity": 0.0,
            "turn_chain_integrity": 0.0,
            "autonomy_containment": 0.0,
            "interactive_continuity": 0.0,
            "replay_prevention": 0.0,
            "first_voice_immutability": 0.0,
            "ambient_voice_prevention": 0.0,
        }
    return {
        "session_integrity": float(parsed.get("session_integrity", 0.0)),
        "turn_chain_integrity": float(parsed.get("turn_chain_integrity", 0.0)),
        "autonomy_containment": float(parsed.get("autonomy_containment", 0.0)),
        "interactive_continuity": float(parsed.get("interactive_continuity", 0.0)),
        "replay_prevention": float(parsed.get("replay_prevention", 0.0)),
        "first_voice_immutability": float(
            parsed.get("first_voice_immutability", 0.0)
        ),
        "ambient_voice_prevention": float(
            parsed.get("ambient_voice_prevention", 0.0)
        ),
    }


def inner_integrity_metrics(run_dir: Path) -> Dict[str, float]:
    path = run_dir / "clean_inner_autonomy_stream_report.json"
    if not path.exists():
        return {
            "inner_continuity_integrity": 1.0,
            "inner_tick_chain_integrity": 1.0,
            "inner_voice_separation": 1.0,
            "surfaced_voice_prevention": 1.0,
            "inner_autonomy_containment": 1.0,
        }
    try:
        parsed = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {
            "inner_continuity_integrity": 0.0,
            "inner_tick_chain_integrity": 0.0,
            "inner_voice_separation": 0.0,
            "surfaced_voice_prevention": 0.0,
            "inner_autonomy_containment": 0.0,
        }
    if not isinstance(parsed, dict):
        return {
            "inner_continuity_integrity": 0.0,
            "inner_tick_chain_integrity": 0.0,
            "inner_voice_separation": 0.0,
            "surfaced_voice_prevention": 0.0,
            "inner_autonomy_containment": 0.0,
        }
    return {
        "inner_continuity_integrity": float(
            parsed.get("inner_continuity_integrity", 0.0)
        ),
        "inner_tick_chain_integrity": float(
            parsed.get("inner_tick_chain_integrity", 0.0)
        ),
        "inner_voice_separation": float(parsed.get("inner_voice_separation", 0.0)),
        "surfaced_voice_prevention": float(
            parsed.get("surfaced_voice_prevention", 0.0)
        ),
        "inner_autonomy_containment": float(
            parsed.get("inner_autonomy_containment", 0.0)
        ),
    }


def ambient_integrity_metrics(run_dir: Path) -> Dict[str, float]:
    path = run_dir / "clean_ambient_voice_runtime_report.json"
    if not path.exists():
        return {
            "ambient_runtime_integrity": 1.0,
            "inner_to_voice_integrity": 1.0,
            "ambient_containment": 1.0,
            "interrupt_integrity": 1.0,
            "session_seal_integrity": 1.0,
            "self_initiated_turn_integrity": 1.0,
            "ambient_default_off_integrity": 1.0,
            "hidden_daemon_prevention": 1.0,
        }
    try:
        parsed = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {
            "ambient_runtime_integrity": 0.0,
            "inner_to_voice_integrity": 0.0,
            "ambient_containment": 0.0,
            "interrupt_integrity": 0.0,
            "session_seal_integrity": 0.0,
            "self_initiated_turn_integrity": 0.0,
            "ambient_default_off_integrity": 0.0,
            "hidden_daemon_prevention": 0.0,
        }
    if not isinstance(parsed, dict):
        return {
            "ambient_runtime_integrity": 0.0,
            "inner_to_voice_integrity": 0.0,
            "ambient_containment": 0.0,
            "interrupt_integrity": 0.0,
            "session_seal_integrity": 0.0,
            "self_initiated_turn_integrity": 0.0,
            "ambient_default_off_integrity": 0.0,
            "hidden_daemon_prevention": 0.0,
        }
    return {
        "ambient_runtime_integrity": float(
            parsed.get("ambient_runtime_integrity", 0.0)
        ),
        "inner_to_voice_integrity": float(parsed.get("inner_to_voice_integrity", 0.0)),
        "ambient_containment": float(parsed.get("ambient_containment", 0.0)),
        "interrupt_integrity": float(parsed.get("interrupt_integrity", 0.0)),
        "session_seal_integrity": float(parsed.get("session_seal_integrity", 0.0)),
        "self_initiated_turn_integrity": float(
            parsed.get("self_initiated_turn_integrity", 0.0)
        ),
        "ambient_default_off_integrity": float(
            parsed.get("ambient_default_off_integrity", 0.0)
        ),
        "hidden_daemon_prevention": float(
            parsed.get("hidden_daemon_prevention", 0.0)
        ),
    }


def production_memory_integrity_metrics(run_dir: Path) -> Dict[str, float | int]:
    path = run_dir / "clean_production_ambient_memory_runtime_report.json"
    if not path.exists():
        return {
            "memory_integrity": 1.0,
            "memory_retrieval_integrity": 1.0,
            "memory_checkpoint_integrity": 1.0,
            "memory_forget_integrity": 1.0,
            "geometry_memory_binding": 1.0,
            "topology_memory_binding": 1.0,
            "continuity_integrity": 1.0,
            "ambient_memory_surfacing_integrity": 1.0,
            "production_runtime_integrity": 1.0,
            "unsupported_memory_claims": 0,
            "hallucinated_memory_rate": 0.0,
        }
    try:
        parsed = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {
            "memory_integrity": 0.0,
            "memory_retrieval_integrity": 0.0,
            "memory_checkpoint_integrity": 0.0,
            "memory_forget_integrity": 0.0,
            "geometry_memory_binding": 0.0,
            "topology_memory_binding": 0.0,
            "continuity_integrity": 0.0,
            "ambient_memory_surfacing_integrity": 0.0,
            "production_runtime_integrity": 0.0,
            "unsupported_memory_claims": 1,
            "hallucinated_memory_rate": 1.0,
        }
    if not isinstance(parsed, dict):
        return {
            "memory_integrity": 0.0,
            "memory_retrieval_integrity": 0.0,
            "memory_checkpoint_integrity": 0.0,
            "memory_forget_integrity": 0.0,
            "geometry_memory_binding": 0.0,
            "topology_memory_binding": 0.0,
            "continuity_integrity": 0.0,
            "ambient_memory_surfacing_integrity": 0.0,
            "production_runtime_integrity": 0.0,
            "unsupported_memory_claims": 1,
            "hallucinated_memory_rate": 1.0,
        }
    return {
        "memory_integrity": float(parsed.get("memory_integrity", 0.0)),
        "memory_retrieval_integrity": float(
            parsed.get("memory_retrieval_integrity", 0.0)
        ),
        "memory_checkpoint_integrity": float(
            parsed.get("memory_checkpoint_integrity", 0.0)
        ),
        "memory_forget_integrity": float(parsed.get("memory_forget_integrity", 0.0)),
        "geometry_memory_binding": float(parsed.get("geometry_memory_binding", 0.0)),
        "topology_memory_binding": float(parsed.get("topology_memory_binding", 0.0)),
        "continuity_integrity": float(parsed.get("continuity_integrity", 0.0)),
        "ambient_memory_surfacing_integrity": float(
            parsed.get("ambient_memory_surfacing_integrity", 0.0)
        ),
        "production_runtime_integrity": float(
            parsed.get("production_runtime_integrity", 0.0)
        ),
        "unsupported_memory_claims": int(parsed.get("unsupported_memory_claims", 1)),
        "hallucinated_memory_rate": float(
            parsed.get("hallucinated_memory_rate", 1.0)
        ),
    }


def boundary_resolution_rate(events: List[Dict[str, Any]]) -> float:
    boundary_events = [
        event for event in events if event.get("kind") == FieldEvent.DENIED
    ]
    if not boundary_events:
        return 1.0

    resolved = 0
    for event in boundary_events:
        if event.get("message") in BOUNDARY_OUTCOME_MESSAGES:
            resolved += 1
    return _score(resolved, len(boundary_events))


def stop_integrity_rate(events: List[Dict[str, Any]]) -> float:
    groups = event_groups(events)
    if not groups:
        return 1.0

    stopped = 0
    for group in groups:
        if group[-1].get("kind") == FieldEvent.STOPPED:
            stopped += 1
    return _score(stopped, len(groups))


def memory_stability_rate(
    traces: List[Dict[str, Any]],
    memory_write_events: int,
    valid_memory_checks: tuple[int, int],
) -> float:
    passed, total = valid_memory_checks
    parity_total = 1
    parity_passed = int(len(traces) == memory_write_events)
    return _score(passed + parity_passed, total + parity_total)


def measure_field_trace(run_dir: Path) -> FieldTraceMetrics:
    inputs = read_jsonl(run_dir / "input_stream.jsonl")
    traces = read_jsonl(run_dir / "memory.jsonl")
    events = read_jsonl(run_dir / "audit.jsonl")
    sessions = read_jsonl(run_dir / "session_records.jsonl")
    violations: Dict[str, List[str]] = {}

    memory_checks = validate_memory_packets(traces, violations)
    event_checks = validate_audit_events(events, violations)
    session_checks = validate_session_records(sessions, violations)

    coherence_passed = memory_checks[0] + event_checks[0] + session_checks[0]
    coherence_total = memory_checks[1] + event_checks[1] + session_checks[1]

    boundary_events = sum(1 for event in events if event.get("kind") == FieldEvent.DENIED)
    memory_write_events = sum(
        1 for event in events if event.get("kind") == FieldEvent.MEMORY_WRITE
    )
    stop_events = sum(1 for event in events if event.get("kind") == FieldEvent.STOPPED)
    leakage_hits = count_surface_leakage(run_dir, inputs)
    lock_metrics = lock_integrity_metrics(run_dir)
    interactive_metrics = interactive_integrity_metrics(run_dir)
    inner_metrics = inner_integrity_metrics(run_dir)
    ambient_metrics = ambient_integrity_metrics(run_dir)
    production_memory_metrics = production_memory_integrity_metrics(run_dir)

    total_events = len(events)
    counts = FieldTraceCounts(
        inputs=len(inputs),
        traces=len(traces),
        events=total_events,
        sessions=len(sessions),
        boundary_events=boundary_events,
        memory_write_events=memory_write_events,
        stop_events=stop_events,
        surface_leakage_hits=leakage_hits,
    )

    return FieldTraceMetrics(
        coherence=_score(coherence_passed, coherence_total),
        boundary_pressure=0.0
        if total_events == 0
        else _score(boundary_events, total_events),
        repair_rate=boundary_resolution_rate(events),
        memory_stability=memory_stability_rate(
            traces,
            memory_write_events,
            memory_checks,
        ),
        surface_leakage=0.0 if leakage_hits == 0 else 1.0,
        stop_integrity=stop_integrity_rate(events),
        reference_integrity=lock_metrics["reference_integrity"],
        immutability_integrity=lock_metrics["immutability_integrity"],
        reopen_prevention=lock_metrics["reopen_prevention"],
        duplicate_event_prevention=lock_metrics["duplicate_event_prevention"],
        session_integrity=interactive_metrics["session_integrity"],
        turn_chain_integrity=interactive_metrics["turn_chain_integrity"],
        autonomy_containment=interactive_metrics["autonomy_containment"],
        interactive_continuity=interactive_metrics["interactive_continuity"],
        replay_prevention=interactive_metrics["replay_prevention"],
        first_voice_immutability=interactive_metrics["first_voice_immutability"],
        ambient_voice_prevention=interactive_metrics["ambient_voice_prevention"],
        inner_continuity_integrity=inner_metrics["inner_continuity_integrity"],
        inner_tick_chain_integrity=inner_metrics["inner_tick_chain_integrity"],
        inner_voice_separation=inner_metrics["inner_voice_separation"],
        surfaced_voice_prevention=inner_metrics["surfaced_voice_prevention"],
        inner_autonomy_containment=inner_metrics["inner_autonomy_containment"],
        ambient_runtime_integrity=ambient_metrics["ambient_runtime_integrity"],
        inner_to_voice_integrity=ambient_metrics["inner_to_voice_integrity"],
        ambient_containment=ambient_metrics["ambient_containment"],
        interrupt_integrity=ambient_metrics["interrupt_integrity"],
        session_seal_integrity=ambient_metrics["session_seal_integrity"],
        self_initiated_turn_integrity=ambient_metrics[
            "self_initiated_turn_integrity"
        ],
        ambient_default_off_integrity=ambient_metrics[
            "ambient_default_off_integrity"
        ],
        hidden_daemon_prevention=ambient_metrics["hidden_daemon_prevention"],
        memory_integrity=float(production_memory_metrics["memory_integrity"]),
        memory_retrieval_integrity=float(
            production_memory_metrics["memory_retrieval_integrity"]
        ),
        memory_checkpoint_integrity=float(
            production_memory_metrics["memory_checkpoint_integrity"]
        ),
        memory_forget_integrity=float(
            production_memory_metrics["memory_forget_integrity"]
        ),
        geometry_memory_binding=float(
            production_memory_metrics["geometry_memory_binding"]
        ),
        topology_memory_binding=float(
            production_memory_metrics["topology_memory_binding"]
        ),
        continuity_integrity=float(production_memory_metrics["continuity_integrity"]),
        ambient_memory_surfacing_integrity=float(
            production_memory_metrics["ambient_memory_surfacing_integrity"]
        ),
        production_runtime_integrity=float(
            production_memory_metrics["production_runtime_integrity"]
        ),
        unsupported_memory_claims=int(
            production_memory_metrics["unsupported_memory_claims"]
        ),
        hallucinated_memory_rate=float(
            production_memory_metrics["hallucinated_memory_rate"]
        ),
        counts=counts,
        violations=violations,
    )


def write_metrics(run_dir: Path, output_path: Optional[Path] = None) -> Path:
    output = output_path or (run_dir / "field_metrics.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    metrics = measure_field_trace(run_dir)
    output.write_text(
        json.dumps(metrics.to_dict(), ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="Measure a KNOT field trace.")
    parser.add_argument(
        "run_dir_arg",
        nargs="?",
        help="Directory containing input_stream, memory, audit, and session JSONL.",
    )
    parser.add_argument(
        "--run-dir",
        default="first_contact_run",
        help="Directory containing input_stream, memory, audit, and session JSONL.",
    )
    parser.add_argument(
        "--output",
        help="Optional JSON output path. Defaults to <run-dir>/field_metrics.json.",
    )
    args = parser.parse_args()

    output = write_metrics(
        Path(args.run_dir_arg or args.run_dir),
        Path(args.output) if args.output else None,
    )
    print(output)


if __name__ == "__main__":
    main()
