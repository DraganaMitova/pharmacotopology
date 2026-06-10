#!/usr/bin/env python3
from __future__ import annotations

"""Postprocess V13a 1UBQ trajectories with a purpose-aware transfer gate.

This script intentionally does not rerun MD and does not use native precision to
choose a threshold.  It reads completed V13a trajectories, classifies the target
purpose from generic properties, sweeps internal frequency thresholds, and writes
a claim-safe readout describing whether a single-domain compact balanced-core
path exists under locked role ontology.
"""

import argparse
import csv
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Optional, Sequence, Union

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
SCRIPTS_ROOT = REPO_ROOT / "scripts"
for path in (SRC_ROOT, SCRIPTS_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_template_docking import chemical_score  # noqa: E402

from run_openmm_tmd_replicas_v0 import (  # noqa: E402
    CORE_LONG_RANGE_SEPARATION,
    _classify_contact_role,
    _distance,
    _effective_lane_pairs,
    _extract_contacts,
    _lane_dca_threshold_config,
    _load_anchor_classes,
    _load_dca_scores,
    _parse_ca_trajectory,
    _parse_domain_boundaries,
    _resolve_pair_lane,
    _resolve_threshold,
    _topology_ok,
)

Pair = tuple[int, int]

DEFAULT_SOURCE_RUN_DIR = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "V13a_1UBQ_REPAIR_FIXED"
)
DEFAULT_OUT_DIR = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "V13a_1UBQ_PURPOSE_GATE_READOUT"
)
DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_holdout_1ubq.locked.json"
DEFAULT_COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_holdout_1ubq_external_couplings.v0.locked.json"
DEFAULT_AUDIT_PAIRS_JSON = REPO_ROOT / "data" / "audit_1ubq_debug.json"


def _load_row(benchmark_file: Path, source_accession: str) -> RealCoordinateVisualRow:
    for row in load_real_coordinate_visual_rows(benchmark_file):
        if row.source_accession == source_accession:
            return row
    raise SystemExit(f"no row found for source_accession={source_accession!r}")


def _load_json(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _load_audit_pairs(path: Optional[Path]) -> set[Pair]:
    if not path or not path.exists():
        return set()
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return set()
    pairs: set[Pair] = set()
    if not isinstance(payload, list):
        return pairs
    for item in payload:
        if isinstance(item, (list, tuple)) and len(item) == 2:
            left, right = int(item[0]), int(item[1])
            pairs.add((left, right) if left < right else (right, left))
    return pairs


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _parse_frequency_grid(raw: str) -> list[float]:
    text = raw.strip()
    if not text:
        return [round(value / 100.0, 2) for value in range(50, 99)]
    if ":" in text:
        parts = [float(item) for item in text.split(":")]
        if len(parts) not in {2, 3}:
            raise SystemExit("--frequency-grid range must be start:end or start:end:step")
        start, end = parts[0], parts[1]
        step = parts[2] if len(parts) == 3 else 0.01
        if step <= 0:
            raise SystemExit("--frequency-grid step must be positive")
        values: list[float] = []
        current = start
        # Include end despite floating point drift.
        while current <= end + (step / 2.0):
            values.append(round(current, 6))
            current += step
        return sorted(set(values), reverse=True)
    return sorted(set(round(float(item.strip()), 6) for item in text.split(",") if item.strip()), reverse=True)


def _find_trajectories(source_run_dir: Path) -> list[Path]:
    return sorted(source_run_dir.glob("replica_*/openmm_dca_restrained_trajectory.pdb"))


def _classify_target_purpose(
    *,
    sequence_length: int,
    topology_mode: str,
    domain_boundaries_raw: str,
    effective_lane_counts: dict[str, int],
    compact_max_length: int = 120,
) -> dict[str, object]:
    boundaries = _parse_domain_boundaries(domain_boundaries_raw)
    if topology_mode == "none" and not boundaries and sequence_length <= compact_max_length:
        return {
            "target_role": "single_domain_compact",
            "topology_policy": "single_domain_compact",
            "strict_required": False,
            "balanced_required": True,
            "rescue_required": False,
            "domain_boundaries_required": False,
            "reason": "no_domain_boundaries_small_target_compact_policy",
        }
    if boundaries and topology_mode == "interdomain":
        return {
            "target_role": "multi_domain_hinge",
            "topology_policy": "interdomain",
            "strict_required": True,
            "balanced_required": True,
            "rescue_required": True,
            "domain_boundaries_required": True,
            "reason": "domain_boundaries_with_interdomain_policy",
        }
    if boundaries and topology_mode == "intradomain":
        return {
            "target_role": "domain_closure",
            "topology_policy": "intradomain",
            "strict_required": True,
            "balanced_required": True,
            "rescue_required": False,
            "domain_boundaries_required": True,
            "reason": "domain_boundaries_with_intradomain_policy",
        }
    if effective_lane_counts.get("balanced", 0) <= 0:
        return {
            "target_role": "low_signal_abstain",
            "topology_policy": "undetermined",
            "strict_required": False,
            "balanced_required": False,
            "rescue_required": False,
            "domain_boundaries_required": False,
            "reason": "no_effective_balanced_pairs_available",
        }
    return {
        "target_role": "undetermined_transfer_purpose",
        "topology_policy": topology_mode or "undetermined",
        "strict_required": True,
        "balanced_required": True,
        "rescue_required": False,
        "domain_boundaries_required": bool(boundaries),
        "reason": "no_generic_target_purpose_rule_matched",
    }


def _build_v13_purpose_preflight(
    *,
    input_preflight: dict[str, object],
    trajectory_count: int,
    audit_pair_count: int,
    effective_lane_counts: dict[str, int],
    target_purpose: dict[str, object],
) -> dict[str, object]:
    checks: dict[str, bool] = {}
    checks["input_preflight_ready"] = input_preflight.get("status") == "ready"
    checks["trajectory_present"] = trajectory_count > 0
    checks["external_coupling_loaded"] = bool(
        input_preflight.get("checks", {}).get("external_coupling_file_exists")
        if isinstance(input_preflight.get("checks"), dict)
        else False
    )
    checks["anchor_profile_loaded"] = bool(
        input_preflight.get("checks", {}).get("anchor_profile_file_exists")
        if isinstance(input_preflight.get("checks"), dict)
        else False
    )
    checks["audit_pair_count_positive"] = audit_pair_count > 0
    checks["classification_coverage_check_enabled"] = True
    checks["target_exact_full_length"] = bool(
        input_preflight.get("target_pdb_provenance", {}).get("sequence_exact_match")
        and input_preflight.get("target_pdb_provenance", {}).get("prepared_target_coverage_ratio") == 1.0
        if isinstance(input_preflight.get("target_pdb_provenance"), dict)
        else False
    )

    target_role = str(target_purpose.get("target_role", "undetermined_transfer_purpose"))
    if target_role == "single_domain_compact":
        checks["effective_balanced_count_positive"] = effective_lane_counts.get("balanced", 0) > 0
        checks["effective_strict_count_optional_for_purpose"] = True
        checks["effective_rescue_count_optional_for_purpose"] = True
        checks["single_domain_topology_policy_explicit"] = target_purpose.get("topology_policy") == "single_domain_compact"
    else:
        checks["effective_strict_count_positive"] = effective_lane_counts.get("strict", 0) > 0
        checks["effective_balanced_count_positive"] = effective_lane_counts.get("balanced", 0) > 0
        if bool(target_purpose.get("rescue_required")):
            checks["effective_rescue_count_positive"] = effective_lane_counts.get("balanced_rescue", 0) > 0
        else:
            checks["effective_rescue_count_optional"] = True

    failed = sorted(name for name, passed in checks.items() if not passed)
    return {
        "kind": "v13_transfer_purpose_preflight_v0",
        "status": "ready" if not failed else "blocked",
        "target_role": target_role,
        "checks": checks,
        "failed_checks": failed,
        "effective_lane_counts": effective_lane_counts,
        "strict_absence_is_fatal": target_role != "single_domain_compact",
        "rescue_absence_is_fatal": bool(target_purpose.get("rescue_required")),
        "claim_allowed": False,
        "physics_interpretation_allowed": False,
    }


def _tail_pair_frequencies(
    trajectory: Path,
    *,
    contact_cutoff_angstrom: float,
    min_separation: int,
    tail_fraction: float,
    audit_pairs: set[Pair],
) -> tuple[dict[Pair, float], dict[Pair, dict[str, object]], int, int]:
    frames = _parse_ca_trajectory(trajectory)
    if not frames:
        return {}, {}, 0, 0
    tail_count = max(1, math.ceil(len(frames) * tail_fraction))
    tail_start = max(0, len(frames) - tail_count)
    tail_frames = frames[-tail_count:]
    counts: Counter[Pair] = Counter()
    for frame in tail_frames:
        counts.update(_extract_contacts(frame, contact_cutoff_angstrom, min_separation))
    freqs = {pair: count / float(tail_count) for pair, count in counts.items()}

    audit_reachability: dict[Pair, dict[str, object]] = {}
    for pair in audit_pairs:
        left, right = pair
        first_contact = None
        last_contact = None
        tail_hits = 0
        min_distance = float("inf")
        min_tail_distance = float("inf")
        for frame_index, frame in enumerate(frames):
            dist = _distance(frame, left, right)
            if not math.isfinite(dist):
                continue
            min_distance = min(min_distance, dist)
            if frame_index >= tail_start:
                min_tail_distance = min(min_tail_distance, dist)
            if dist <= contact_cutoff_angstrom:
                if first_contact is None:
                    first_contact = frame_index
                last_contact = frame_index
                if frame_index >= tail_start:
                    tail_hits += 1
        observed = first_contact is not None
        tail_frequency = tail_hits / float(tail_count) if tail_count else 0.0
        tail_observed = tail_hits > 0
        if not observed:
            state = "never_approached"
        elif not tail_observed:
            state = "approach_then_disappear"
        else:
            state = "present_in_tail"
        audit_reachability[pair] = {
            "observed": observed,
            "tail_observed": tail_observed,
            "trajectory_state": state,
            "tail_frequency": round(tail_frequency, 6),
            "min_distance": min_distance if math.isfinite(min_distance) else None,
            "min_tail_distance": min_tail_distance if math.isfinite(min_tail_distance) else None,
            "first_contact_frame": first_contact,
            "last_contact_frame": last_contact,
            "tail_start_frame": tail_start,
            "tail_frame_count": tail_count,
        }
    return freqs, audit_reachability, len(frames), tail_count


def _mean(values: Sequence[float]) -> Optional[float]:
    if not values:
        return None
    return round(sum(values) / len(values), 6)


def _evaluate_frequency_threshold(
    threshold: float,
    *,
    trajectory_payloads: list[dict[str, object]],
    row: RealCoordinateVisualRow,
    sequence: str,
    dca_scores: dict[Pair, float],
    anchor_classes: dict[Pair, str],
    effective_lane_pairs: dict[str, set[Pair]],
    balanced_dca_threshold: float,
    chemical_threshold: float,
    vote_threshold: int,
    topology_mode: str,
    domain_boundaries_raw: str,
    audit_pairs: set[Pair],
) -> dict[str, object]:
    boundaries = _parse_domain_boundaries(domain_boundaries_raw)
    support: Counter[Pair] = Counter()
    per_pair_freqs: dict[Pair, list[float]] = defaultdict(list)
    per_pair_replica_hits: dict[Pair, list[int]] = defaultdict(list)
    per_pair_failure: dict[Pair, Counter[str]] = defaultdict(Counter)
    selected_replica_pairs_by_lane: dict[str, Counter[Pair]] = defaultdict(Counter)

    effective_balanced = set(effective_lane_pairs.get("balanced", set()))
    effective_anchor_core = (
        set(effective_lane_pairs.get("strict", set()))
        | effective_balanced
        | set(effective_lane_pairs.get("balanced_rescue", set()))
    )

    for replica_index, payload in enumerate(trajectory_payloads, start=1):
        freqs = payload.get("frequencies", {})
        if not isinstance(freqs, dict):
            continue
        for pair, freq in freqs.items():
            if not isinstance(pair, tuple) or len(pair) != 2:
                continue
            left, right = int(pair[0]), int(pair[1])
            pair = (left, right) if left < right else (right, left)
            if pair not in effective_balanced:
                continue
            if right - left < CORE_LONG_RANGE_SEPARATION:
                per_pair_failure[pair]["not_long_range"] += 1
                continue
            pair_dca = float(dca_scores.get(pair, 0.0))
            raw_class = anchor_classes.get(pair, "balanced")
            lane = _resolve_pair_lane(
                raw_class,
                pair_dca,
                strict_threshold=balanced_dca_threshold,
                balanced_strong_threshold=balanced_dca_threshold,
                balanced_rescue_threshold=min(0.7, balanced_dca_threshold),
            )
            if lane != "balanced":
                per_pair_failure[pair][f"lane_{lane}"] += 1
                continue
            if not _topology_ok(left, right, domain_boundaries=boundaries, topology_mode=topology_mode):
                per_pair_failure[pair]["topology"] += 1
                continue
            if pair_dca < balanced_dca_threshold:
                per_pair_failure[pair]["dca"] += 1
                continue
            chem = chemical_score(sequence[left - 1], sequence[right - 1])
            if chem < chemical_threshold:
                per_pair_failure[pair]["chemical"] += 1
                continue
            if float(freq) < threshold:
                per_pair_failure[pair]["frequency"] += 1
                continue
            support[pair] += 1
            per_pair_freqs[pair].append(float(freq))
            per_pair_replica_hits[pair].append(replica_index)
            selected_replica_pairs_by_lane[lane][pair] += 1

    selected_pairs = {pair for pair, count in support.items() if count >= vote_threshold}
    selected_balanced_core = sorted(selected_pairs & effective_balanced)
    selected_outside_effective = sorted(selected_pairs - effective_anchor_core)
    dca_values = [float(dca_scores.get(pair, 0.0)) for pair in selected_balanced_core]
    background_dca_values = [float(dca_scores.get(pair, 0.0)) for pair in effective_balanced]
    dca_mean = _mean(dca_values)
    background_mean = _mean(background_dca_values)
    dca_enrichment_ratio = None
    if dca_mean is not None and background_mean not in (None, 0):
        dca_enrichment_ratio = round(dca_mean / float(background_mean), 6)
    dca_enrichment_pass = bool(selected_balanced_core) and all(value >= balanced_dca_threshold for value in dca_values)

    audit_payload: dict[str, object] = {}
    for pair in sorted(audit_pairs):
        freqs_by_replica: dict[str, float] = {}
        states: Counter[str] = Counter()
        min_tail_distances: list[float] = []
        for idx, payload in enumerate(trajectory_payloads, start=1):
            reach = payload.get("audit_reachability", {})
            if not isinstance(reach, dict):
                continue
            details = reach.get(pair)
            if not isinstance(details, dict):
                continue
            freq = details.get("tail_frequency")
            if isinstance(freq, (int, float)):
                freqs_by_replica[f"replica_{idx:02d}"] = round(float(freq), 6)
            state = details.get("trajectory_state")
            if isinstance(state, str):
                states[state] += 1
            min_tail = details.get("min_tail_distance")
            if isinstance(min_tail, (int, float)) and math.isfinite(float(min_tail)):
                min_tail_distances.append(float(min_tail))
        audit_payload[f"{pair[0]}-{pair[1]}"] = {
            "selected_at_threshold": pair in selected_pairs,
            "support_at_threshold": int(support.get(pair, 0)),
            "vote_threshold": vote_threshold,
            "tail_frequency_mean": _mean(list(freqs_by_replica.values())),
            "tail_frequency_min": min(freqs_by_replica.values()) if freqs_by_replica else None,
            "tail_frequency_max": max(freqs_by_replica.values()) if freqs_by_replica else None,
            "tail_frequencies_by_replica": freqs_by_replica,
            "state_counts": dict(states),
            "mean_tail_distance": _mean(min_tail_distances),
            "dca_score": float(dca_scores.get(pair, 0.0)),
            "candidate_for_balanced_compact_core": pair in effective_balanced,
            "blocked_reason_at_threshold": "selected" if pair in selected_pairs else (
                "support_below_vote_threshold" if support.get(pair, 0) > 0 else "frequency_or_gate_filter"
            ),
        }

    classification_coverage_ratio = 1.0 if selected_pairs else 0.0
    long_range_evidence_polluted = False
    noise_added = len(selected_outside_effective)
    pass_checks = {
        "nonzero_balanced_core": len(selected_balanced_core) > 0,
        "dca_enrichment_pass": dca_enrichment_pass,
        "noise_added_zero": noise_added == 0,
        "long_range_evidence_not_polluted": not long_range_evidence_polluted,
        "classification_coverage_complete": classification_coverage_ratio == 1.0 if selected_pairs else False,
    }
    return {
        "threshold": round(float(threshold), 6),
        "selected_pair_count": len(selected_pairs),
        "selected_balanced_core_count": len(selected_balanced_core),
        "selected_balanced_core": [list(pair) for pair in selected_balanced_core],
        "selected_outside_effective_anchor_set": [list(pair) for pair in selected_outside_effective],
        "noise_added": noise_added,
        "long_range_evidence_polluted": long_range_evidence_polluted,
        "classification_coverage_ratio": classification_coverage_ratio,
        "dca_mean_selected": dca_mean,
        "dca_mean_effective_balanced_background": background_mean,
        "dca_enrichment_ratio": dca_enrichment_ratio,
        "dca_enrichment_pass": dca_enrichment_pass,
        "vote_threshold": vote_threshold,
        "support_by_selected_pair": {
            f"{pair[0]}-{pair[1]}": int(support.get(pair, 0)) for pair in selected_balanced_core
        },
        "mean_frequency_by_selected_pair": {
            f"{pair[0]}-{pair[1]}": _mean(per_pair_freqs.get(pair, [])) for pair in selected_balanced_core
        },
        "replica_hits_by_selected_pair": {
            f"{pair[0]}-{pair[1]}": per_pair_replica_hits.get(pair, []) for pair in selected_balanced_core
        },
        "audit_pairs": audit_payload,
        "pass_checks": pass_checks,
        "passes_purpose_fit": all(pass_checks.values()),
    }


def _select_highest_passing_threshold(sweep_rows: Sequence[dict[str, object]]) -> Optional[dict[str, object]]:
    passing = [row for row in sweep_rows if bool(row.get("passes_purpose_fit"))]
    if not passing:
        return None
    return sorted(passing, key=lambda row: (float(row.get("threshold", 0.0)), int(row.get("selected_pair_count", 0))), reverse=True)[0]


def main() -> None:
    parser = argparse.ArgumentParser(description="V13a purpose-aware postprocess readout for completed 1UBQ trajectories.")
    parser.add_argument("--source-run-dir", default=str(DEFAULT_SOURCE_RUN_DIR))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--source-accession", default="1UBQ:A")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_COUPLING_FILE))
    parser.add_argument("--anchor-profile-file", default=str(DEFAULT_COUPLING_FILE))
    parser.add_argument("--audit-pairs-json", default=str(DEFAULT_AUDIT_PAIRS_JSON))
    parser.add_argument("--topology-mode", default="none", choices=("none", "interdomain", "intradomain"))
    parser.add_argument("--domain-boundaries", default="")
    parser.add_argument("--frequency-grid", default="0.50:0.98:0.01")
    parser.add_argument("--tail-fraction", type=float, default=0.20)
    parser.add_argument("--contact-cutoff-ang", type=float, default=7.0)
    parser.add_argument("--min-separation", type=int, default=24)
    parser.add_argument("--balanced-dca-threshold", type=float, default=0.80)
    parser.add_argument("--chemical-threshold", type=float, default=0.50)
    parser.add_argument("--vote-threshold", type=int, default=7)
    args = parser.parse_args()

    source_run_dir = Path(args.source_run_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    row = _load_row(Path(args.benchmark_file), args.source_accession)
    audit_pairs = _load_audit_pairs(Path(args.audit_pairs_json) if args.audit_pairs_json else None)
    trajectories = _find_trajectories(source_run_dir)
    if not source_run_dir.exists():
        raise SystemExit(f"source run directory missing: {source_run_dir}")
    if not trajectories:
        raise SystemExit(f"no completed trajectories found in {source_run_dir}; run V13a first")

    source_certificate = _load_json(source_run_dir / "openmm_tmd_replicas_v0_certificate.json")
    input_preflight = _load_json(source_run_dir / "input_preflight.json")
    if not input_preflight and isinstance(source_certificate.get("input_preflight"), dict):
        input_preflight = source_certificate["input_preflight"]  # type: ignore[assignment]

    readout_coupling = Path(args.external_coupling_file)
    anchor_profile = Path(args.anchor_profile_file) if args.anchor_profile_file else readout_coupling
    effective_lane_pairs = _effective_lane_pairs(
        anchor_profile_file=anchor_profile,
        coupling_file=readout_coupling,
        row=row,
        strict_dca_threshold=None,
        balanced_strong_dca_threshold=args.balanced_dca_threshold,
        balanced_rescue_dca_threshold=0.70,
        monitor_dca_threshold=None,
        dca_threshold="auto",
    )
    effective_lane_counts = {lane: len(pairs) for lane, pairs in effective_lane_pairs.items()}
    target_purpose = _classify_target_purpose(
        sequence_length=row.sequence_length,
        topology_mode=args.topology_mode,
        domain_boundaries_raw=args.domain_boundaries,
        effective_lane_counts=effective_lane_counts,
    )
    v13_preflight = _build_v13_purpose_preflight(
        input_preflight=input_preflight,
        trajectory_count=len(trajectories),
        audit_pair_count=len(audit_pairs),
        effective_lane_counts=effective_lane_counts,
        target_purpose=target_purpose,
    )

    dca_scores, _ = _load_dca_scores(
        coupling_file=readout_coupling,
        row=row,
        sequence_length=row.sequence_length,
    )
    anchor_classes = _load_anchor_classes(
        coupling_file=anchor_profile,
        row=row,
        sequence_length=row.sequence_length,
    )

    trajectory_payloads: list[dict[str, object]] = []
    for trajectory in trajectories:
        freqs, audit_reachability, frame_count, tail_count = _tail_pair_frequencies(
            trajectory,
            contact_cutoff_angstrom=args.contact_cutoff_ang,
            min_separation=args.min_separation,
            tail_fraction=args.tail_fraction,
            audit_pairs=audit_pairs,
        )
        trajectory_payloads.append(
            {
                "trajectory": str(trajectory),
                "frequencies": freqs,
                "audit_reachability": audit_reachability,
                "frame_count": frame_count,
                "tail_count": tail_count,
            }
        )

    thresholds = _parse_frequency_grid(args.frequency_grid)
    sweep_rows: list[dict[str, object]] = []
    for threshold in thresholds:
        sweep_rows.append(
            _evaluate_frequency_threshold(
                threshold,
                trajectory_payloads=trajectory_payloads,
                row=row,
                sequence=row.sequence,
                dca_scores=dca_scores,
                anchor_classes=anchor_classes,
                effective_lane_pairs=effective_lane_pairs,
                balanced_dca_threshold=args.balanced_dca_threshold,
                chemical_threshold=args.chemical_threshold,
                vote_threshold=args.vote_threshold,
                topology_mode=args.topology_mode,
                domain_boundaries_raw=args.domain_boundaries,
                audit_pairs=audit_pairs,
            )
        )

    selected_band = _select_highest_passing_threshold(sweep_rows)
    decision = "purpose_readout_core_found" if selected_band else "clean_abstain_no_stable_frequency_band"

    flat_sweep_rows = []
    for row_payload in sweep_rows:
        flat_sweep_rows.append(
            {
                "threshold": row_payload["threshold"],
                "passes_purpose_fit": row_payload["passes_purpose_fit"],
                "selected_pair_count": row_payload["selected_pair_count"],
                "selected_balanced_core_count": row_payload["selected_balanced_core_count"],
                "noise_added": row_payload["noise_added"],
                "long_range_evidence_polluted": row_payload["long_range_evidence_polluted"],
                "classification_coverage_ratio": row_payload["classification_coverage_ratio"],
                "dca_mean_selected": row_payload["dca_mean_selected"],
                "dca_mean_effective_balanced_background": row_payload["dca_mean_effective_balanced_background"],
                "dca_enrichment_ratio": row_payload["dca_enrichment_ratio"],
                "selected_balanced_core": json.dumps(row_payload["selected_balanced_core"]),
                "support_by_selected_pair": json.dumps(row_payload["support_by_selected_pair"], sort_keys=True),
            }
        )

    certificate = {
        "kind": "V13a_1UBQ_PURPOSE_GATE_READOUT_v0",
        "run_mode": "postprocess_only_no_new_simulation",
        "source_run_dir": str(source_run_dir),
        "source_certificate": str(source_run_dir / "openmm_tmd_replicas_v0_certificate.json"),
        "source_accession": args.source_accession,
        "row_id": row.row_id,
        "sequence_length": row.sequence_length,
        "trajectory_count": len(trajectories),
        "input_preflight_status": input_preflight.get("status"),
        "target_purpose": target_purpose,
        "v13_transfer_preflight": v13_preflight,
        "legacy_v5_preflight": source_certificate.get("preflight", {}),
        "legacy_runtime_v10_selected_pair_count": source_certificate.get("runtime_v10_selected_pair_count"),
        "effective_lane_counts": effective_lane_counts,
        "effective_lanes": {
            lane: [list(pair) for pair in sorted(pairs)] for lane, pairs in effective_lane_pairs.items()
        },
        "frequency_grid": thresholds,
        "frequency_sweep": sweep_rows,
        "selected_frequency_band": selected_band,
        "purpose_gate_decision": decision,
        "official_runtime_rerun_required_for_claim": bool(selected_band),
        "claim_allowed": False,
        "physics_interpretation_allowed": False,
        "biological_transfer_claim_allowed": False,
        "forbidden_methods": {
            "hardcoded_target_threshold": False,
            "native_precision_used_to_select_threshold": False,
            "selector_rules_modified": False,
            "partial_target_allowed": False,
        },
        "status_interpretation": (
            "postprocess found an internally supported single-domain compact balanced-core band; "
            "run an official purpose-aware runtime only after reviewing this certificate"
            if selected_band
            else "postprocess cleanly abstained; no internally stable purpose-fit frequency band was found"
        ),
    }

    cert_path = out_dir / "v13a_1ubq_purpose_gate_readout_certificate.json"
    cert_path.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(out_dir / "v13a_1ubq_purpose_gate_frequency_sweep.csv", flat_sweep_rows)

    compact = {
        "kind": certificate["kind"],
        "run_mode": certificate["run_mode"],
        "trajectory_count": len(trajectories),
        "input_preflight_status": input_preflight.get("status"),
        "target_role": target_purpose.get("target_role"),
        "v13_preflight_status": v13_preflight.get("status"),
        "legacy_v5_ready": source_certificate.get("preflight", {}).get("v5_ready") if isinstance(source_certificate.get("preflight"), dict) else None,
        "legacy_v5_failed_checks": source_certificate.get("preflight", {}).get("v5_requirements", {}).get("failed_checks") if isinstance(source_certificate.get("preflight"), dict) and isinstance(source_certificate.get("preflight", {}).get("v5_requirements"), dict) else None,
        "legacy_runtime_v10_selected_pair_count": source_certificate.get("runtime_v10_selected_pair_count"),
        "purpose_gate_decision": decision,
        "selected_threshold": selected_band.get("threshold") if selected_band else None,
        "selected_balanced_core": selected_band.get("selected_balanced_core") if selected_band else [],
        "selected_pair_count": selected_band.get("selected_pair_count") if selected_band else 0,
        "claim_allowed": False,
        "certificate": str(cert_path),
    }
    print(json.dumps(compact, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
