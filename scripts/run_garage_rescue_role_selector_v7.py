#!/usr/bin/env python3
from __future__ import annotations

"""Postprocess-only V7 rescue role selector repair.

Runs on an existing V5-like run directory and rebuilds per-replica stable pairs
using role-aware classification (local / medium / long-range) on top of lane-aware
admission checks, without re-running OpenMM trajectories.
"""

import argparse
import json
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Optional, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "scripts"))

from run_openmm_tmd_replicas_v0 import (  # type: ignore[import-not-found]
    DEFAULT_BENCHMARK_FILE,
    DEFAULT_EXTERNAL_COUPLING_FILE,
    _apply_degree_cap,
    _control_hits,
    _evaluate_audit_pair_trajectory_reachability,
    _extract_contacts,
    _lane_chemistry_config,
    _lane_dca_threshold_config,
    _load_anchor_classes,
    _load_dca_scores,
    _load_row,
    _parse_ca_trajectory,
    _parse_domain_boundaries,
    _resolve_pair_lane,
    _resolve_threshold,
    _topology_ok,
    _vote_pairs,
)
from pharmacotopology.folding_template_docking import chemical_score

Pair = tuple[int, int]


CORE_LONG_RANGE_SEPARATION = 24
MEDIUM_SUPPORT_MIN_SEPARATION = 8


def _safe_ratio(count: int, total: int) -> float:
    if total <= 0:
        return 0.0
    return round(count / total, 6)


def _coerce_pair(item: Sequence[Any]) -> Pair:
    left = int(item[0])
    right = int(item[1])
    return (left, right) if left < right else (right, left)


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        loaded = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return loaded if isinstance(loaded, dict) else {}


def _summarize_pair_reachability(
    reachability_by_replica: dict[str, object],
    selected: bool,
) -> dict[str, object]:
    states: Counter[str] = Counter()
    per_replica: dict[str, bool] = {}
    min_distances: list[float] = []
    tail_distances: list[float] = []
    tail_frequencies: list[float] = []
    tail_presence_count = 0

    for replica_name, details in sorted(reachability_by_replica.items()):
        if not isinstance(details, dict):
            continue
        state = str(details.get("trajectory_state", "missing_from_tail"))
        states[state] += 1

        tail_present = bool(details.get("tail_observed", False))
        per_replica[replica_name] = tail_present
        if tail_present:
            tail_presence_count += 1

        min_distance = details.get("min_distance")
        if isinstance(min_distance, (int, float)) and math.isfinite(float(min_distance)):
            min_distances.append(float(min_distance))

        min_tail_distance = details.get("min_tail_distance")
        if isinstance(min_tail_distance, (int, float)) and math.isfinite(float(min_tail_distance)):
            tail_distances.append(float(min_tail_distance))

        tail_frequency = details.get("tail_frequency")
        if isinstance(tail_frequency, (int, float)) and math.isfinite(float(tail_frequency)):
            tail_frequencies.append(float(tail_frequency))

    return {
        "tail_contact_frequency": round(sum(tail_frequencies) / len(tail_frequencies), 6)
        if tail_frequencies
        else 0.0,
        "min_distance_over_trajectory": min(min_distances) if min_distances else None,
        "mean_tail_distance": round(sum(tail_distances) / len(tail_distances), 6)
        if tail_distances
        else None,
        "per_replica_tail_presence": per_replica,
        "approached_then_lost": states["approach_then_disappear"] > 0,
        "never_approached": states["never_approached"] > 0,
        "selected_in_final_vote": bool(selected),
        "state_counts": dict(states),
        "tail_presence_count": tail_presence_count,
    }


def _aggregate_pair_status(values: Sequence[str]) -> str:
    if not values:
        return "missing_from_tail"
    priority = [
        "selected",
        "degree_cap",
        "control_overlap",
        "lane_min_separation",
        "topology",
        "dca",
        "chemical",
        "frequency",
        "approach_then_disappear",
        "never_approached",
        "missing_from_tail",
    ]
    for status in priority:
        if status in values:
            return status
    return str(values[0])


def _extract_target_pairs(
    source_payload: dict[str, Any],
    explicit_path: Optional[Path],
) -> list[Pair]:
    pairs: list[Pair] = []

    if explicit_path is not None and explicit_path.exists():
        try:
            raw_explicit = json.loads(explicit_path.read_text(encoding="utf-8"))
        except Exception:
            raw_explicit = []
        for item in raw_explicit or []:
            if isinstance(item, (list, tuple)) and len(item) == 2:
                pairs.append(_coerce_pair(item))

    if pairs:
        return sorted(set(pairs))

    pair_audit = source_payload.get("pair_audit")
    if isinstance(pair_audit, dict) and pair_audit:
        for raw_key in pair_audit:
            if not isinstance(raw_key, str):
                continue
            left_text, sep, right_text = raw_key.partition("-")
            if not sep:
                continue
            if left_text.isdigit() and right_text.isdigit():
                pairs.append(_coerce_pair((int(left_text), int(right_text))))

    if pairs:
        return sorted(set(pairs))

    preflight = source_payload.get("preflight")
    if isinstance(preflight, dict):
        effective = preflight.get("effective")
        if isinstance(effective, dict):
            rescue_pairs = effective.get("balanced_rescue", [])
            if isinstance(rescue_pairs, list):
                for item in rescue_pairs:
                    if isinstance(item, (list, tuple)) and len(item) == 2:
                        pairs.append(_coerce_pair(item))

    return sorted(set(pairs))


def _build_stable_pairs_with_lane_min_separation(
    trajectory: Path,
    *,
    row,
    sequence: str,
    tail_fraction: float,
    contact_cutoff_angstrom: float,
    min_separation_by_lane: dict[str, int],
    frequency_threshold: str | float,
    dca_threshold: str | float,
    chemical_threshold: float,
    topology_mode: str,
    domain_boundaries_raw: str,
    max_degree: int,
    control_count: int,
    max_control_overlap_fraction: float,
    coupling_file: Path,
    anchor_profile_file: Optional[Path],
    audit_pairs: set[Pair],
    forced_lane_pairs: Optional[dict[Pair, str]] = None,
    strict_dca_threshold: Optional[str | float] = None,
    balanced_strong_dca_threshold: Optional[str | float] = None,
    balanced_rescue_dca_threshold: Optional[str | float] = None,
    monitor_dca_threshold: Optional[str | float] = None,
    strict_chemical_threshold: Optional[float] = None,
    balanced_strong_chemical_threshold: Optional[float] = None,
    balanced_rescue_chemical_threshold: Optional[float] = None,
    monitor_chemical_threshold: Optional[float] = None,
    min_separation_inclusive: bool = True,
) -> tuple[
    dict[Pair, dict[str, float]],
    dict[str, object],
    dict[Pair, str],
    dict[Pair, str],
    dict[Pair, dict[str, object]],
]:
    frames = _parse_ca_trajectory(trajectory)
    if not frames:
        return {}, {
            "tail_count": 0,
            "considered_frames": 0,
            "lane_rejections": {},
        }, {}, {}, {}

    tail_count = max(1, math.ceil(len(frames) * tail_fraction))
    final_frames = frames[-tail_count:]
    required_frames = len(final_frames)

    # Global pass: extract raw contacts from tail with no min-separation filter.
    pair_counts: Counter[Pair] = Counter()
    for frame in final_frames:
        pair_counts.update(_extract_contacts(frame, contact_cutoff_angstrom, 0))

    reachability = _evaluate_audit_pair_trajectory_reachability(
        frames,
        audit_pairs=set(audit_pairs),
        contact_cutoff_angstrom=contact_cutoff_angstrom,
        tail_count=tail_count,
    )

    dca_scores, _dca_classes = _load_dca_scores(
        coupling_file=coupling_file,
        row=row,
        sequence_length=row.sequence_length,
    )
    raw_classes = _load_anchor_classes(
        coupling_file=anchor_profile_file or coupling_file,
        row=row,
        sequence_length=row.sequence_length,
    )
    if not raw_classes:
        raw_classes = _load_anchor_classes(
            coupling_file=coupling_file,
            row=row,
            sequence_length=row.sequence_length,
        )

    freq_values = [count / required_frames for count in pair_counts.values()] if pair_counts else []
    resolved_frequency_threshold = _resolve_threshold(
        frequency_threshold,
        freq_values,
        default=0.5,
    )

    lane_dca_thresholds = _lane_dca_threshold_config(
        dca_threshold,
        strict_dca_threshold,
        balanced_strong_dca_threshold,
        balanced_rescue_dca_threshold,
        monitor_dca_threshold,
    )
    lane_chem_thresholds = _lane_chemistry_config(
        chemical_threshold,
        strict_chemical_threshold,
        balanced_strong_chemical_threshold,
        balanced_rescue_chemical_threshold,
        monitor_chemical_threshold,
    )

    dca_values = [dca_scores.get(pair, 0.0) for pair in pair_counts]
    resolved_dca_threshold = _resolve_threshold(
        lane_dca_thresholds["unknown"],
        dca_values,
        default=0.0,
    )
    dca_thresholds = {
        lane: _resolve_threshold(value, dca_values, default=resolved_dca_threshold)
        for lane, value in lane_dca_thresholds.items()
    }
    resolved_chem_threshold = _resolve_threshold(
        chemical_threshold,
        list(pair_counts.values()) or [],
        default=0.0,
    )
    chemical_thresholds = {
        lane: _resolve_threshold(value, [], default=resolved_chem_threshold)
        for lane, value in lane_chem_thresholds.items()
    }

    boundaries = _parse_domain_boundaries(domain_boundaries_raw)
    lane_rejections: dict[str, Counter[str]] = {
        "lane_min_separation": Counter(),
        "frequency": Counter(),
        "chemical": Counter(),
        "dca": Counter(),
        "topology": Counter(),
        "control": Counter(),
    }
    boundary_24_rejections: Counter[str] = Counter()
    lane_degree_rejections: Counter[str] = Counter()
    lane_pre_counts: Counter[str] = Counter()

    pair_rejections: dict[Pair, str] = {pair: "missing_from_tail" for pair in audit_pairs}
    lane_lookup: dict[Pair, str] = {}
    pre_control: list[tuple[Pair, dict[str, float]]] = []

    # Assign every audited pair that appears in trajectory a first-stage rejection, then
    # overwrite as we pass each gate.
    for pair in sorted(set(pair_counts) & audit_pairs):
        pair_rejections[pair] = _evaluate_audit_pair_trajectory_reachability(
            [next(iter(final_frames))],
            audit_pairs={pair},
            contact_cutoff_angstrom=contact_cutoff_angstrom,
            tail_count=1,
        ).get(
            pair,
            {"trajectory_state": "missing_from_tail"},
        ).get("trajectory_state", "missing_from_tail")

    default_min_sep = max(int(v) for v in min_separation_by_lane.values()) if min_separation_by_lane else 0

    for (left, right), count in sorted(pair_counts.items()):
        if left > len(sequence) or right > len(sequence) or left >= right:
            continue
        seq_sep = right - left

        raw_pair_class = raw_classes.get((left, right), "")
        dca_score = dca_scores.get((left, right), 0.0)
        forced_lane = forced_lane_pairs.get((left, right)) if forced_lane_pairs else None
        if forced_lane is not None:
            lane = forced_lane
        else:
            lane = _resolve_pair_lane(
                raw_pair_class,
                dca_score,
                strict_threshold=dca_thresholds["strict"],
                balanced_strong_threshold=dca_thresholds["balanced"],
                balanced_rescue_threshold=dca_thresholds["balanced_rescue"],
            )
        lane_lookup[(left, right)] = lane
        lane_pre_counts[lane] += 1

        min_sep = min_separation_by_lane.get(lane, default_min_sep)
        if (seq_sep < min_sep) or (not min_separation_inclusive and seq_sep <= min_sep):
            lane_rejections["lane_min_separation"][lane] += 1
            if (left, right) in audit_pairs:
                pair_rejections[(left, right)] = "lane_min_separation"
                if min_sep == CORE_LONG_RANGE_SEPARATION and seq_sep == CORE_LONG_RANGE_SEPARATION:
                    boundary_24_rejections[lane] += 1
            continue

        freq = count / required_frames
        if freq < resolved_frequency_threshold:
            lane_rejections["frequency"][lane] += 1
            if (left, right) in audit_pairs:
                pair_rejections[(left, right)] = "frequency"
            continue

        chem = chemical_score(sequence[left - 1], sequence[right - 1])
        if chem < chemical_thresholds[lane]:
            lane_rejections["chemical"][lane] += 1
            if (left, right) in audit_pairs:
                pair_rejections[(left, right)] = "chemical"
            continue

        if dca_score < dca_thresholds[lane]:
            lane_rejections["dca"][lane] += 1
            if (left, right) in audit_pairs:
                pair_rejections[(left, right)] = "dca"
            continue

        if not _topology_ok(
            left,
            right,
            domain_boundaries=boundaries,
            topology_mode=topology_mode,
): 
            lane_rejections["topology"][lane] += 1
            if (left, right) in audit_pairs:
                pair_rejections[(left, right)] = "topology"
            continue

        pre_control.append(
            (
                (left, right),
                {
                    "count": float(count),
                    "frequency": freq,
                    "chemical_score": chem,
                    "dca_score": dca_score,
                    "pair_lane": lane,
                },
            )
        )

    pre_control_count = len(pre_control)
    auto_control_count = control_count
    if auto_control_count <= 0:
        auto_control_count = 0 if pre_control_count <= 0 else min(6, int(math.sqrt(pre_control_count)))

    control_hits = _control_hits(
        row=row,
        selected_pairs=[pair for pair, _ in pre_control],
        candidate_pairs=tuple(sorted(pair_counts)),
        control_count=auto_control_count,
    )

    surviving: list[tuple[Pair, dict[str, float]]] = []
    control_removed = 0
    for pair, payload in pre_control:
        lane = payload["pair_lane"]
        overlap = control_hits.get(pair, 0) / auto_control_count if auto_control_count else 0.0
        payload["control_overlap_fraction"] = overlap
        payload["control_gap"] = 1.0 - overlap
        if overlap > max_control_overlap_fraction:
            control_removed += 1
            if pair in audit_pairs:
                pair_rejections[pair] = "control_overlap"
            lane_rejections["control"][lane] += 1
            continue
        payload["survival_score"] = (
            0.50 * payload["frequency"]
            + 0.25 * payload["chemical_score"]
            + 0.25 * payload["dca_score"]
        )
        surviving.append((pair, payload))

    surviving.sort(
        key=lambda item: (
            -item[1]["survival_score"],
            -item[1]["frequency"],
            -item[1]["dca_score"],
            item[0][0],
            item[0][1],
        ),
    )

    before_degree = len(surviving)
    surviving = _apply_degree_cap(
        surviving,
        max_degree=max_degree,
        sequence_length=row.sequence_length,
    )
    degree_culled = before_degree - len(surviving)

    stable: dict[Pair, dict[str, float]] = {
        pair: payload for pair, payload in surviving
    }
    stable_lookup = {pair: payload for pair, payload in surviving}
    for pair, payload in [(pair, payload) for pair, payload in pre_control]:
        if pair not in stable_lookup:
            if pair in audit_pairs:
                pair_rejections[pair] = "degree_cap"
            lane_degree_rejections[payload["pair_lane"]] += 1

    for pair in audit_pairs:
        if pair in stable_lookup:
            pair_rejections[pair] = "selected"

    lane_pairs = {
        "strict": [list(pair) for pair, payload in stable.items() if payload["pair_lane"] == "strict"],
        "balanced": [list(pair) for pair, payload in stable.items() if payload["pair_lane"] == "balanced"],
        "balanced_rescue": [list(pair) for pair, payload in stable.items() if payload["pair_lane"] == "balanced_rescue"],
        "monitor": [list(pair) for pair, payload in stable.items() if payload["pair_lane"] == "monitor"],
        "unknown": [list(pair) for pair, payload in stable.items() if payload["pair_lane"] == "unknown"],
    }

    pair_reachability = {
        pair: details
        for pair, details in reachability.items()
    }

    metadata = {
        "tail_count": tail_count,
        "considered_frames": required_frames,
        "raw_contact_pair_count": len(pair_counts),
        "resolved_frequency_threshold": resolved_frequency_threshold,
        "resolved_dca_threshold": resolved_dca_threshold,
        "resolved_dca_threshold_by_lane": dca_thresholds,
        "resolved_chemical_threshold_by_lane": chemical_thresholds,
        "pre_control_pair_count": pre_control_count,
        "control_removed_pair_count": control_removed,
        "degree_culled_pair_count": degree_culled,
        "control_count": auto_control_count,
        "max_control_overlap_fraction": max_control_overlap_fraction,
        "topology_mode": topology_mode,
        "domain_boundaries": domain_boundaries_raw,
        "max_degree": max_degree,
        "min_separation_by_lane": dict(min_separation_by_lane),
        "lane_candidate_counts": {name: value for name, value in lane_pre_counts.items()},
        "seq_sep_24_boundary_failures_by_lane": {lane: int(count) for lane, count in boundary_24_rejections.items()},
        "lane_rejections": {
            reason: {lane: int(count) for lane, count in by_lane.items()}
            for reason, by_lane in lane_rejections.items()
        },
        "lane_degree_rejections": dict(lane_degree_rejections),
        "lane_pairs": {name: sorted(items) for name, items in lane_pairs.items()},
        "pair_trajectory_reachability": {
            f"{left}-{right}": details for (left, right), details in pair_reachability.items()
        },
        "pair_rejection_reasons": [
            [left, right, reason] for (left, right), reason in sorted(pair_rejections.items())
        ],
    }

    return stable, metadata, pair_rejections, lane_lookup, pair_reachability


def _classify_contact_role(lane: str, sequence_sep: int) -> str:
    if sequence_sep >= CORE_LONG_RANGE_SEPARATION:
        if lane == "balanced_rescue":
            return "border_long_range_rescue"
        if lane in {"strict", "balanced"}:
            return "true_long_range_core"
        return "diagnostic_shell"
    if sequence_sep >= MEDIUM_SUPPORT_MIN_SEPARATION:
        return "medium_support"
    return "local_support"


def _is_long_range_role(role: str) -> bool:
    return role == "true_long_range_core"


def _is_border_rescue_role(role: str) -> bool:
    return role == "border_long_range_rescue"


def _is_rescue_support_role(role: str) -> bool:
    return role in {"border_long_range_rescue", "medium_support", "local_support"}


def _collect_replica_map(items: list[dict[Pair, str]], pair: Pair) -> dict[str, str]:
    values: dict[str, str] = {}
    for idx, payload in enumerate(items):
        if pair in payload:
            values[f"replica_{idx + 1:02d}"] = payload[pair]
    return values


def _collect_replica_lanes(items: list[dict[Pair, str]], pair: Pair) -> dict[str, str]:
    values: dict[str, str] = {}
    for idx, payload in enumerate(items):
        if pair in payload:
            values[f"replica_{idx + 1:02d}"] = payload[pair]
    return values


def _collect_replica_reachability(
    items: list[dict[Pair, dict[str, object]]],
    pair: Pair,
) -> dict[str, dict[str, object]]:
    values: dict[str, dict[str, object]] = {}
    for idx, payload in enumerate(items):
            if pair in payload:
                values[f"replica_{idx + 1:02d}"] = payload[pair]
    return values


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run V7 rescue role selector postprocess on existing V5 trajectories.",
    )
    parser.add_argument(
        "--source-run-dir",
        default=str(
            REPO_ROOT
            / "first_contact_clean_pharmacotopology_layer_run"
            / "GARAGE_TAIL_REACHABILITY_WITH_STRICT_GUARD_V5"
        ),
    )
    parser.add_argument(
        "--source-certificate",
        default="",
        help="Optional explicit source certificate JSON path.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(
            REPO_ROOT
            / "first_contact_clean_pharmacotopology_layer_run"
            / "GARAGE_RESCUE_ROLE_SELECTOR_V7"
        ),
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_EXTERNAL_COUPLING_FILE))
    parser.add_argument("--readout-coupling-file", default="")
    parser.add_argument("--anchor-profile-file", default="")
    parser.add_argument("--audit-pairs-json", default="")

    parser.add_argument("--strict-min-separation", type=int, default=24)
    parser.add_argument("--balanced-min-separation", type=int, default=24)
    parser.add_argument("--rescue-min-separation", type=int, default=12)
    parser.add_argument("--monitor-min-separation", type=int, default=24)
    parser.add_argument("--long-range-min-separation", type=int, default=CORE_LONG_RANGE_SEPARATION)
    parser.add_argument(
        "--strict-seq-boundary",
        action="store_true",
        default=False,
        help="Use strict (<=) sequence-separation min-separation boundary checks.",
    )

    parser.add_argument("--write-json", default="")

    args = parser.parse_args()

    source_run = Path(args.source_run_dir)
    if not source_run.exists():
        raise SystemExit(f"source run directory missing: {source_run}")

    if args.source_certificate:
        certificate_path = Path(args.source_certificate).expanduser()
        if not certificate_path.exists():
            certificate_path = source_run / "openmm_tmd_replicas_v0_certificate.json"
    else:
        certificate_path = source_run / "openmm_tmd_replicas_v0_certificate.json"
    source_payload = _read_json(certificate_path)

    source_accession = source_payload.get("source_accession", "4AKE:A")
    row = _load_row(Path(args.benchmark_file), source_accession)

    topology_mode = str(source_payload.get("topology_mode", "interdomain"))
    domain_boundaries = str(source_payload.get("domain_boundaries", "1-60,61-110,111-170,171-214"))
    tail_fraction = float(source_payload.get("tail_fraction", 0.20))
    contact_cutoff = float(source_payload.get("contact_cutoff_angstrom", 7.0))

    readout_coupling = Path(args.readout_coupling_file) if args.readout_coupling_file else Path(args.external_coupling_file)
    if not readout_coupling.exists():
        readout_coupling = Path(args.external_coupling_file)
    anchor_profile = Path(args.anchor_profile_file) if args.anchor_profile_file else None

    dca_threshold = source_payload.get("dca_threshold", "auto")
    frequency_threshold = source_payload.get("frequency_threshold", "auto")
    strict_dca_threshold = source_payload.get("strict_dca_threshold")
    balanced_strong_dca_threshold = source_payload.get("balanced_strong_dca_threshold", 0.80)
    balanced_rescue_dca_threshold = source_payload.get("balanced_rescue_dca_threshold", 0.70)
    monitor_dca_threshold = source_payload.get("monitor_dca_threshold")

    chemical_threshold = float(source_payload.get("chemical_threshold", 0.5))
    strict_chemical_threshold = source_payload.get("strict_chemical_threshold")
    balanced_strong_chemical_threshold = source_payload.get("balanced_strong_chemical_threshold")
    balanced_rescue_chemical_threshold = source_payload.get("balanced_rescue_chemical_threshold")
    monitor_chemical_threshold = source_payload.get("monitor_chemical_threshold")

    max_degree = int(source_payload.get("max_degree", 4))
    control_count = int(source_payload.get("control_count", 0))
    max_control_overlap_fraction = float(source_payload.get("max_control_overlap_fraction", 0.25))

    vote_threshold = int(source_payload.get("vote_threshold", 3))
    vote_min_average_chem = float(source_payload.get("vote_min_average_chemical", 0.5))
    strict_vote_threshold = source_payload.get("strict_vote_threshold")
    balanced_vote_threshold = source_payload.get("balanced_vote_threshold")
    rescue_vote_threshold = source_payload.get("balanced_rescue_vote_threshold")
    monitor_vote_threshold = source_payload.get("monitor_vote_threshold")
    strict_vote_min_average_chem = source_payload.get("strict_vote_min_average_chemical")
    balanced_vote_min_average_chem = source_payload.get("balanced_vote_min_average_chemical")
    rescue_vote_min_average_chem = source_payload.get("balanced_rescue_vote_min_average_chemical")
    monitor_vote_min_average_chem = source_payload.get("monitor_vote_min_average_chemical")

    target_pairs = _extract_target_pairs(
        source_payload,
        Path(args.audit_pairs_json) if args.audit_pairs_json else None,
    )
    if not target_pairs:
        raise SystemExit("target rescue pairs are empty; provide --audit-pairs-json or source pair_audit")
    target_pairs_set = set(target_pairs)

    forced_lane_pairs: dict[Pair, str] = {}
    preflight = source_payload.get("preflight")
    if isinstance(preflight, dict):
        preflight_effective = preflight.get("effective")
        if isinstance(preflight_effective, dict):
            for lane_name, raw_pairs in preflight_effective.items():
                if not isinstance(raw_pairs, list):
                    continue
                for item in raw_pairs:
                    if isinstance(item, (list, tuple)) and len(item) == 2:
                        forced_lane_pairs[_coerce_pair(item)] = str(lane_name)

    min_sep_before = {
        "strict": int(args.strict_min_separation),
        "balanced": int(args.balanced_min_separation),
        "balanced_rescue": int(args.rescue_min_separation),
        "monitor": int(args.monitor_min_separation),
        "unknown": int(args.long_range_min_separation),
    }
    min_sep_after = {
        "strict": int(args.long_range_min_separation),
        "balanced": int(args.long_range_min_separation),
        "balanced_rescue": int(args.rescue_min_separation),
        "monitor": int(args.monitor_min_separation),
        "unknown": int(args.long_range_min_separation),
    }

    replica_dirs = sorted(
        [
            item
            for item in source_run.glob("replica_*")
            if item.is_dir() and (item / "openmm_dca_restrained_trajectory.pdb").exists()
        ]
    )
    if not replica_dirs:
        raise SystemExit(f"no replica trajectories found in {source_run}")

    before_payloads: list[dict[Pair, dict[str, float]]] = []
    after_payloads: list[dict[Pair, dict[str, float]]] = []
    before_metadata: list[dict[str, object]] = []
    after_metadata: list[dict[str, object]] = []
    before_reasons: list[dict[Pair, str]] = []
    after_reasons: list[dict[Pair, str]] = []
    before_lanes: list[dict[Pair, str]] = []
    after_lanes: list[dict[Pair, str]] = []
    before_reachability: list[dict[Pair, dict[str, object]]] = []
    after_reachability: list[dict[Pair, dict[str, object]]] = []

    for replica_dir in replica_dirs:
        trajectory = replica_dir / "openmm_dca_restrained_trajectory.pdb"
        if not trajectory.exists():
            continue

        before_stable, before_meta, before_reason_map, before_lane_map, before_reach = _build_stable_pairs_with_lane_min_separation(
            trajectory=trajectory,
            row=row,
            sequence=row.sequence,
            tail_fraction=tail_fraction,
            contact_cutoff_angstrom=contact_cutoff,
            min_separation_by_lane=min_sep_before,
            frequency_threshold=frequency_threshold,
            dca_threshold=dca_threshold,
            chemical_threshold=chemical_threshold,
            topology_mode=topology_mode,
            domain_boundaries_raw=domain_boundaries,
            max_degree=max_degree,
            control_count=control_count,
            max_control_overlap_fraction=max_control_overlap_fraction,
            coupling_file=readout_coupling,
            anchor_profile_file=anchor_profile,
            audit_pairs=target_pairs_set,
            forced_lane_pairs=forced_lane_pairs,
            strict_dca_threshold=strict_dca_threshold,
            balanced_strong_dca_threshold=balanced_strong_dca_threshold,
            balanced_rescue_dca_threshold=balanced_rescue_dca_threshold,
            monitor_dca_threshold=monitor_dca_threshold,
            strict_chemical_threshold=strict_chemical_threshold,
            balanced_strong_chemical_threshold=balanced_strong_chemical_threshold,
            balanced_rescue_chemical_threshold=balanced_rescue_chemical_threshold,
            monitor_chemical_threshold=monitor_chemical_threshold,
            min_separation_inclusive=not args.strict_seq_boundary,
        )

        after_stable, after_meta, after_reason_map, after_lane_map, after_reach = _build_stable_pairs_with_lane_min_separation(
            trajectory=trajectory,
            row=row,
            sequence=row.sequence,
            tail_fraction=tail_fraction,
            contact_cutoff_angstrom=contact_cutoff,
            min_separation_by_lane=min_sep_after,
            frequency_threshold=frequency_threshold,
            dca_threshold=dca_threshold,
            chemical_threshold=chemical_threshold,
            topology_mode=topology_mode,
            domain_boundaries_raw=domain_boundaries,
            max_degree=max_degree,
            control_count=control_count,
            max_control_overlap_fraction=max_control_overlap_fraction,
            coupling_file=readout_coupling,
            anchor_profile_file=anchor_profile,
            audit_pairs=target_pairs_set,
            forced_lane_pairs=forced_lane_pairs,
            strict_dca_threshold=strict_dca_threshold,
            balanced_strong_dca_threshold=balanced_strong_dca_threshold,
            balanced_rescue_dca_threshold=balanced_rescue_dca_threshold,
            monitor_dca_threshold=monitor_dca_threshold,
            strict_chemical_threshold=strict_chemical_threshold,
            balanced_strong_chemical_threshold=balanced_strong_chemical_threshold,
            balanced_rescue_chemical_threshold=balanced_rescue_chemical_threshold,
            monitor_chemical_threshold=monitor_chemical_threshold,
            min_separation_inclusive=not args.strict_seq_boundary,
        )

        before_payloads.append(before_stable)
        after_payloads.append(after_stable)
        before_metadata.append(before_meta)
        after_metadata.append(after_meta)
        before_reasons.append(before_reason_map)
        after_reasons.append(after_reason_map)
        before_lanes.append(before_lane_map)
        after_lanes.append(after_lane_map)
        before_reachability.append(before_reach)
        after_reachability.append(after_reach)

    if not before_payloads or not after_payloads:
        raise SystemExit("no successful trajectories could be postprocessed")

    before_voted, _, before_lane_support, before_lane_selected, _, before_fail_reasons = _vote_pairs(
        before_payloads,
        vote_threshold=vote_threshold,
        min_average_chem=vote_min_average_chem,
        strict_vote_threshold=strict_vote_threshold,
        balanced_vote_threshold=balanced_vote_threshold,
        balanced_rescue_vote_threshold=rescue_vote_threshold,
        monitor_vote_threshold=monitor_vote_threshold,
        strict_min_average_chem=strict_vote_min_average_chem,
        balanced_strong_min_average_chem=balanced_vote_min_average_chem,
        balanced_rescue_min_average_chem=rescue_vote_min_average_chem,
        monitor_min_average_chem=monitor_vote_min_average_chem,
    )
    after_voted, _, after_lane_support, after_lane_selected, _, after_fail_reasons = _vote_pairs(
        after_payloads,
        vote_threshold=vote_threshold,
        min_average_chem=vote_min_average_chem,
        strict_vote_threshold=strict_vote_threshold,
        balanced_vote_threshold=balanced_vote_threshold,
        balanced_rescue_vote_threshold=rescue_vote_threshold,
        monitor_vote_threshold=monitor_vote_threshold,
        strict_min_average_chem=strict_vote_min_average_chem,
        balanced_strong_min_average_chem=balanced_vote_min_average_chem,
        balanced_rescue_min_average_chem=rescue_vote_min_average_chem,
        monitor_min_average_chem=monitor_vote_min_average_chem,
    )

    before_selected = {(left, right) for left, right, _, _, _ in before_voted}
    after_selected = {(left, right) for left, right, _, _, _ in after_voted}

    before_selected_by_lane: dict[str, set[Pair]] = {
        "strict": set(),
        "balanced": set(),
        "balanced_rescue": set(),
        "monitor": set(),
        "unknown": set(),
    }
    after_selected_by_lane: dict[str, set[Pair]] = {
        "strict": set(),
        "balanced": set(),
        "balanced_rescue": set(),
        "monitor": set(),
        "unknown": set(),
    }
    for left, right, _, _, lane in before_voted:
        before_selected_by_lane.setdefault(lane, set()).add((left, right))
    for left, right, _, _, lane in after_voted:
        after_selected_by_lane.setdefault(lane, set()).add((left, right))

    before_support: Counter[Pair] = Counter()
    after_support: Counter[Pair] = Counter()
    for payload in before_payloads:
        for pair in payload:
            before_support[pair] += 1
    for payload in after_payloads:
        for pair in payload:
            after_support[pair] += 1

    before_selected_by_role: dict[str, set[Pair]] = {}
    after_selected_by_role: dict[str, set[Pair]] = {}

    before_support_failures_by_role: dict[str, int] = Counter()
    after_support_failures_by_role: dict[str, int] = Counter()
    boundary_failures_before: dict[str, int] = Counter()
    boundary_failures_after: dict[str, int] = Counter()
    tail_presence_sum_by_role: dict[str, int] = Counter()
    tail_presence_pair_count_by_role: dict[str, int] = Counter()
    pair_reports: list[dict[str, Any]] = []

    contact_role_counts_before: Counter[str] = Counter()
    contact_role_counts_after: Counter[str] = Counter()

    for pair in sorted(target_pairs_set):
        left, right = pair
        sequence_sep = right - left

        before_rep_status = _collect_replica_map(before_reasons, pair)
        after_rep_status = _collect_replica_map(after_reasons, pair)
        before_lane_rep = _collect_replica_lanes(before_lanes, pair)
        after_lane_rep = _collect_replica_lanes(after_lanes, pair)

        before_status = _aggregate_pair_status(list(before_rep_status.values()))
        after_status = _aggregate_pair_status(list(after_rep_status.values()))

        before_lanes_seen = sorted(set(before_lane_rep.values()))
        after_lanes_seen = sorted(set(after_lane_rep.values()))
        before_lane = before_lanes_seen[0] if before_lanes_seen else "monitor"
        after_lane = after_lanes_seen[0] if after_lanes_seen else "monitor"

        before_role = _classify_contact_role(before_lane, sequence_sep)
        after_role = _classify_contact_role(after_lane, sequence_sep)

        contact_role_counts_before[before_role] += 1
        contact_role_counts_after[after_role] += 1

        after_pair_reach = _collect_replica_reachability(after_reachability, pair)
        reach_summary = _summarize_pair_reachability(
            after_pair_reach,
            selected=(pair in after_selected),
        )
        before_fail_reason = before_fail_reasons.get(pair, "missing_from_tail")
        after_fail_reason = after_fail_reasons.get(pair, "missing_from_tail")

        before_support_failures_by_role[before_role] += 1 if before_fail_reason != "selected" else 0
        after_support_failures_by_role[after_role] += 1 if after_fail_reason != "selected" else 0
        if before_status == "lane_min_separation" and sequence_sep == CORE_LONG_RANGE_SEPARATION:
            boundary_failures_before[before_lane] += 1
        if after_status == "lane_min_separation" and sequence_sep == CORE_LONG_RANGE_SEPARATION:
            boundary_failures_after[after_lane] += 1

        if pair in before_selected:
            before_selected_by_role.setdefault(before_role, set()).add(pair)
        if pair in after_selected:
            after_selected_by_role.setdefault(after_role, set()).add(pair)

        pair_reports.append(
            {
                "pair": [left, right],
                "sequence_separation": sequence_sep,
                "lane_admission_status_before": before_status,
                "lane_admission_status_after": after_status,
                "before_lane": before_lane,
                "after_lane": after_lane,
                "before_role_classification": before_role,
                "after_role_classification": after_role,
                "long_range_ready_by_role": after_status != "lane_min_separation"
                and _is_long_range_role(after_role),
                "border_rescue_ready_by_role": after_status != "lane_min_separation"
                and _is_border_rescue_role(after_role),
                "blocked_by_global_min_sep_before": before_status == "lane_min_separation",
                "blocked_by_global_min_sep_after": after_status == "lane_min_separation",
                "allowed_by_seq_sep_boundary_rule": sequence_sep >= args.long_range_min_separation,
                "allowed_by_rescue_lane_min_sep": (
                    after_status != "lane_min_separation"
                    and after_lane == "balanced_rescue"
                    and sequence_sep >= args.rescue_min_separation
                ),
                "tail_contact_frequency": reach_summary["tail_contact_frequency"],
                "vote_count_after_lane_admission": int(after_support.get(pair, 0)),
                "selected_after_rescue_lane": pair in after_selected_by_lane.get("balanced_rescue", set()),
                "selected_after_lane": pair in after_selected,
                "selected_after_role": pair in after_selected_by_role.get(after_role, set()),
                "selected_after_in_repair_lane": pair in after_selected_by_lane.get("balanced_rescue", set()),
                "selected_by_rescue_ontology": pair in after_selected_by_role.get("border_long_range_rescue", set()),
                "role_classification": after_role,
                "evidence_class": after_role,
                "vote_failure_before": before_fail_reason,
                "vote_failure_after": after_fail_reason,
                "reason_by_replica_before": before_rep_status,
                "reason_by_replica_after": after_rep_status,
                "lane_by_replica_after": after_lane_rep,
                "tail_reachability_by_replica": after_pair_reach,
                "reachability_summary": reach_summary,
            }
        )

        tail_presence_pair_count_by_role[after_role] += 1
        tail_presence_sum_by_role[after_role] += int(reach_summary["tail_presence_count"])

    rescue_pairs_by_category = {
        "local_support_candidates": 0,
        "medium_support_candidates": 0,
        "border_long_range_rescue_candidates": 0,
        "true_long_range_core_candidates": 0,
        "diagnostic_shell_candidates": 0,
    }
    for report in pair_reports:
        role = str(report["before_role_classification"])
        if role == "local_support":
            rescue_pairs_by_category["local_support_candidates"] += 1
        elif role == "medium_support":
            rescue_pairs_by_category["medium_support_candidates"] += 1
        elif role == "border_long_range_rescue":
            rescue_pairs_by_category["border_long_range_rescue_candidates"] += 1
        elif role == "true_long_range_core":
            rescue_pairs_by_category["true_long_range_core_candidates"] += 1
        else:
            rescue_pairs_by_category["diagnostic_shell_candidates"] += 1

    def _role_selected_count(role: str) -> int:
        return len(after_selected_by_role.get(role, set()))

    rescue_candidates_before = (
        rescue_pairs_by_category["local_support_candidates"]
        + rescue_pairs_by_category["medium_support_candidates"]
        + rescue_pairs_by_category["border_long_range_rescue_candidates"]
    )
    rescue_candidates_after = (
        sum(
            1
            for report in pair_reports
            if report["after_role_classification"] in {"local_support", "medium_support", "border_long_range_rescue"}
            and report["lane_admission_status_after"] != "lane_min_separation"
        )
    )

    local_support_selected = _role_selected_count("local_support")
    medium_support_selected = _role_selected_count("medium_support")
    border_long_range_rescue_selected = _role_selected_count("border_long_range_rescue")
    true_long_range_core_selected = _role_selected_count("true_long_range_core")
    diagnostic_shell_selected = _role_selected_count("diagnostic_shell")

    vote_failure_by_role = {
        role: {"before_failures": int(before_support_failures_by_role.get(role, 0)),
               "after_failures": int(after_support_failures_by_role.get(role, 0))}
        for role in sorted(set(before_support_failures_by_role) | set(after_support_failures_by_role))
    }
    all_roles = {
        "local_support",
        "medium_support",
        "border_long_range_rescue",
        "true_long_range_core",
        "diagnostic_shell",
    }
    tail_presence_by_role = {
        role: _safe_ratio(tail_presence_sum_by_role.get(role, 0), tail_presence_pair_count_by_role.get(role, 0))
        for role in all_roles
    }
    rescue_long_range_success = border_long_range_rescue_selected > 0
    rescue_local_success = local_support_selected > 0
    rescue_medium_success = medium_support_selected > 0
    rescue_core_success = true_long_range_core_selected > 0
    rescue_local_not_core = rescue_local_success != rescue_core_success
    rescue_medium_not_core = rescue_medium_success != rescue_core_success

    before_long_range_count = len(before_selected_by_role.get("true_long_range_core", set()))
    after_long_range_count = len(after_selected_by_role.get("true_long_range_core", set()))

    effective_anchor_pairs: set[Pair] = set()
    preflight = source_payload.get("preflight", {})
    if isinstance(preflight, dict):
        effective = preflight.get("effective")
        if isinstance(effective, dict):
            for lane in ("strict", "balanced", "balanced_rescue"):
                raw = effective.get(lane, [])
                if isinstance(raw, list):
                    for item in raw:
                        if isinstance(item, (list, tuple)) and len(item) == 2:
                            effective_anchor_pairs.add(_coerce_pair(item))

    after_pair_in_blueprint = _safe_ratio(
        len(after_selected & effective_anchor_pairs),
        len(after_selected),
    )

    summary: dict[str, object] = {
        "run_type": "GARAGE_RESCUE_ROLE_SELECTOR_V7",
        "run_mode": "postprocess_no_new_simulation",
        "claim_allowed": False,
        "source_run_dir": str(source_run),
        "source_certificate": str(certificate_path),
        "before_min_separation_by_lane": min_sep_before,
        "after_min_separation_by_lane": min_sep_after,
        "long_range_min_separation": args.long_range_min_separation,
        "strict_seq_boundary": args.strict_seq_boundary,
        "rescue_pairs": pair_reports,
        "rescue_candidates_before": rescue_candidates_before,
        "rescue_candidates_after": rescue_candidates_after,
        "rescue_selected_after": border_long_range_rescue_selected + medium_support_selected + local_support_selected,
        "selected_by_role": {
            "local_support_selected": sorted([
                item["pair"] for item in pair_reports if item.get("selected_after_role")
                and item.get("after_role_classification") == "local_support"
            ]),
            "medium_support_selected": sorted([
                item["pair"] for item in pair_reports if item.get("selected_after_role")
                and item.get("after_role_classification") == "medium_support"
            ]),
            "border_long_range_rescue_selected": sorted([
                item["pair"] for item in pair_reports if item.get("selected_after_role")
                and item.get("after_role_classification") == "border_long_range_rescue"
            ]),
            "true_long_range_core_selected": sorted([
                item["pair"] for item in pair_reports if item.get("selected_after_role")
                and item.get("after_role_classification") == "true_long_range_core"
            ]),
        },
        "contact_role_counts": {
            "before": dict(contact_role_counts_before),
            "after": dict(contact_role_counts_after),
        },
        "rescue_pairs_by_category": dict(rescue_pairs_by_category),
        "local_support_selected": local_support_selected,
        "medium_support_selected": medium_support_selected,
        "border_long_range_rescue_selected": border_long_range_rescue_selected,
        "true_long_range_core_selected": true_long_range_core_selected,
        "diagnostic_shell_selected": diagnostic_shell_selected,
        "tail_presence_by_role": tail_presence_by_role,
        "vote_failure_by_role": vote_failure_by_role,
        "rescue_selected_by_role": {
            "local_support_selected": bool(rescue_pairs_by_category["local_support_candidates"] > 0 and local_support_selected > 0),
            "medium_support_selected": bool(rescue_pairs_by_category["medium_support_candidates"] > 0 and medium_support_selected > 0),
            "border_long_range_rescue_selected": bool(
                rescue_pairs_by_category["border_long_range_rescue_candidates"] > 0 and border_long_range_rescue_selected > 0
            ),
        },
        "rescue_local_success": rescue_local_success,
        "rescue_medium_success": rescue_medium_success,
        "rescue_long_range_success": rescue_long_range_success,
        "rescue_local_success_differs_from_folding_core": rescue_local_not_core,
        "rescue_medium_success_differs_from_folding_core": rescue_medium_not_core,
        "min_separation_boundary_failures": {
            "before": dict(boundary_failures_before),
            "after": dict(boundary_failures_after),
            "seq_sep_24_boundary_handling": {
                "value": CORE_LONG_RANGE_SEPARATION,
                "inclusive": not args.strict_seq_boundary,
                "before_24_blocked_pairs": int(sum(boundary_failures_before.values())),
                "after_24_blocked_pairs": int(sum(boundary_failures_after.values())),
                "pairs_unblocked_by_inclusive_boundary": int(
                    sum(boundary_failures_before.values()) - sum(boundary_failures_after.values())
                ),
            },
        },
        "seq_sep_24_boundary_handling": {
            "required": f">={CORE_LONG_RANGE_SEPARATION}",
            "strict_mode": bool(args.strict_seq_boundary),
            "configured": int(args.long_range_min_separation),
        },
        "long_range_selected_before": before_long_range_count,
        "long_range_selected_after": after_long_range_count,
        "long_range_evidence_polluted": (before_long_range_count != after_long_range_count),
        "noise_added": after_lane_selected.get("monitor", 0) + after_lane_selected.get("unknown", 0),
        "pair_in_blueprint": after_pair_in_blueprint,
        "selected_pair_count": len(after_selected),
        "selected_pairs": sorted([list(pair) for pair in sorted(after_selected)]),
        "selected_by_lane_after": {
            "strict": sorted([list(pair) for pair in sorted(after_selected_by_lane["strict"])]),
            "balanced": sorted([list(pair) for pair in sorted(after_selected_by_lane["balanced"])]),
            "balanced_rescue": sorted([list(pair) for pair in sorted(after_selected_by_lane["balanced_rescue"])]),
            "monitor": sorted([list(pair) for pair in sorted(after_selected_by_lane["monitor"])]),
            "unknown": sorted([list(pair) for pair in sorted(after_selected_by_lane["unknown"])]),
        },
        "selected_by_lane_before": {
            "strict": sorted([list(pair) for pair in sorted(before_selected_by_lane["strict"])]),
            "balanced": sorted([list(pair) for pair in sorted(before_selected_by_lane["balanced"])]),
            "balanced_rescue": sorted([list(pair) for pair in sorted(before_selected_by_lane["balanced_rescue"])]),
            "monitor": sorted([list(pair) for pair in sorted(before_selected_by_lane["monitor"])]),
            "unknown": sorted([list(pair) for pair in sorted(before_selected_by_lane["unknown"])]),
        },
        "selected_by_role_before": {
            role: sorted([list(pair) for pair in sorted(before_selected_by_role.get(role, set()))])
            for role in sorted(set(before_selected_by_role))
        },
        "selected_by_role_after": {
            role: sorted([list(pair) for pair in sorted(after_selected_by_role.get(role, set()))])
            for role in sorted(set(after_selected_by_role))
        },
        "lane_vote_support_after": dict(after_lane_support),
        "lane_vote_support_before": dict(before_lane_support),
        "preflight": {
            "source_preflight_ready": bool(source_payload.get("preflight", {}).get("v5_ready", False)),
            "v5_requirements": source_payload.get("preflight", {}).get("v5_requirements", {}),
            "failure_types": source_payload.get("preflight", {}).get("failure_type")
            or source_payload.get("preflight", {}).get("notes", []),
        },
    }

    out_root = Path(args.out_dir)
    out_root.mkdir(parents=True, exist_ok=True)
    out_file = out_root / "garbage_road_selector_v7_postprocess.json"
    if args.write_json:
        out_file = Path(args.write_json)
        out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
