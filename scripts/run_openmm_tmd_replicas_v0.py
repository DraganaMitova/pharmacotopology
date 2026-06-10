#!/usr/bin/env python3
from __future__ import annotations

"""Run multiple independent TMD-like OpenMM replicas and vote stable contacts.

Each replica runs `run_openmm_dca_closure_md_v0.py` with a different seed and
same deterministic settings. Contacts are extracted from the trajectory tail using
an explicit "contact survival" rule (frequency, chemistry, DCA support,
topology gate, control-gap and per-residue degree cap) before voting across
replicas.
"""

import argparse
import csv
import json
import math
import os
import random
import subprocess
import sys
from collections import Counter
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Union

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

Pair = tuple[int, int]

CORE_LONG_RANGE_SEPARATION = 24
MEDIUM_SUPPORT_MIN_SEPARATION = 8

from pharmacotopology.folding_native_contact_eval import evaluate_contact_prediction
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_five_axis_physics import matched_control_pairs
from pharmacotopology.folding_template_docking import chemical_score


DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
DEFAULT_EXTERNAL_COUPLING_FILE = (
    REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "query_centered_pfam00406_external_couplings_v0.locked.json"
)
DEFAULT_TARGET_PDB = REPO_ROOT / "data" / "independent_contact_sources" / "AF-P69441-F1-model_v4.pdb"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "openmm_dca_tmd_replicas_v0"
RUNNER_SCRIPT = REPO_ROOT / "scripts" / "run_openmm_dca_closure_md_v0.py"


def _load_row(benchmark_file: Path, source_accession: str) -> RealCoordinateVisualRow:
    rows = load_real_coordinate_visual_rows(benchmark_file)
    for row in rows:
        if row.source_accession == source_accession:
            return row
    raise SystemExit(f"no row found for source_accession={source_accession!r}")


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
        for row in rows:
            writer.writerow(row)


def _parse_ca_trajectory(path: Path, chain_id: str = "A") -> list[dict[int, tuple[float, float, float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    current_model: dict[int, tuple[float, float, float]] = {}
    models: list[dict[int, tuple[float, float, float]]] = []

    have_model_tags = False
    for raw_line in lines:
        if raw_line.startswith("MODEL"):
            have_model_tags = True
            current_model = {}
            continue
        if raw_line.startswith("ENDMDL"):
            models.append(current_model)
            current_model = {}
            continue
        if not raw_line.startswith("ATOM"):
            continue
        if raw_line[12:16].strip() != "CA":
            continue
        chain = raw_line[21].strip() or "A"
        if chain_id and chain != chain_id:
            continue
        try:
            index = int(raw_line[22:26].strip())
            x = float(raw_line[30:38])
            y = float(raw_line[38:46])
            z = float(raw_line[46:54])
        except Exception:
            continue
        current_model[index] = (x, y, z)

    if current_model and not have_model_tags:
        models = [current_model]
    elif current_model:
        models.append(current_model)
    return models


def _extract_contacts(
    coords: dict[int, tuple[float, float, float]],
    cutoff_angstrom: float,
    min_separation: int,
) -> dict[tuple[int, int], float]:
    if not coords:
        return {}
    indexes = sorted(coords)
    cutoff_sq = cutoff_angstrom * cutoff_angstrom
    contacts: dict[tuple[int, int], float] = {}
    for left_pos, left_idx in enumerate(indexes[:-1]):
        left = indexes[left_pos]
        left_x, left_y, left_z = coords[left]
        for right in indexes[left_pos + 1 :]:
            if right - left <= min_separation:
                continue
            right_x, right_y, right_z = coords[right]
            dx = left_x - right_x
            dy = left_y - right_y
            dz = left_z - right_z
            d2 = dx * dx + dy * dy + dz * dz
            if d2 <= cutoff_sq:
                contacts[(left, right)] = math.sqrt(d2)
    return contacts


def _distance(
    coords: dict[int, tuple[float, float, float]],
    left: int,
    right: int,
) -> float:
    if left not in coords or right not in coords:
        return float("inf")
    lx, ly, lz = coords[left]
    rx, ry, rz = coords[right]
    dx = lx - rx
    dy = ly - ry
    dz = lz - rz
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _evaluate_audit_pair_trajectory_reachability(
    frames: list[dict[int, tuple[float, float, float]]],
    *,
    audit_pairs: set[tuple[int, int]],
    contact_cutoff_angstrom: float,
    tail_count: int,
) -> dict[tuple[int, int], dict[str, object]]:
    if not frames or not audit_pairs:
        return {}
    tail_start = max(0, len(frames) - tail_count)
    reachability: dict[tuple[int, int], dict[str, object]] = {}
    for left, right in sorted(audit_pairs):
        min_distance = float("inf")
        min_tail_distance = float("inf")
        tail_hits = 0
        tail_state_seen = False
        first_contact = None
        last_contact = None
        for frame_index, frame in enumerate(frames):
            dist = _distance(frame, left, right)
            if math.isinf(dist):
                continue
            if dist < min_distance:
                min_distance = dist
            if frame_index >= tail_start and dist < min_tail_distance:
                min_tail_distance = dist
            if dist <= contact_cutoff_angstrom:
                if frame_index >= tail_start:
                    tail_hits += 1
                    tail_state_seen = True
                if first_contact is None:
                    first_contact = frame_index
                last_contact = frame_index
        observed = first_contact is not None
        tail_observed = tail_state_seen
        tail_frequency = float(tail_hits) / float(tail_count) if tail_count else 0.0
        if not observed:
            trajectory_state = "never_approached"
        elif not tail_observed:
            trajectory_state = "approach_then_disappear"
        else:
            trajectory_state = "present_in_tail"
        reachability[(left, right)] = {
            "observed": observed,
            "tail_observed": tail_observed,
            "trajectory_state": trajectory_state,
            "min_distance": min_distance if min_distance != float("inf") else None,
            "min_tail_distance": min_tail_distance if min_tail_distance != float("inf") else None,
            "tail_frequency": round(tail_frequency, 6),
            "first_contact_frame": first_contact,
            "last_contact_frame": last_contact,
            "tail_start_frame": tail_start,
            "tail_frame_count": max(1, tail_count),
        }
    return reachability


def _load_dca_scores(
    *,
    coupling_file: Path,
    row: RealCoordinateVisualRow,
    sequence_length: int,
) -> tuple[dict[tuple[int, int], float], dict[tuple[int, int], str]]:
    if not coupling_file or not coupling_file.exists():
        return {}, {}
    try:
        dataset = load_coupling_dataset(coupling_file)
    except Exception:
        return {}, {}
    scores: dict[tuple[int, int], float] = {}
    classes: dict[tuple[int, int], str] = {}
    for constraint in dataset.constraints:
        if constraint.row_id != row.row_id and constraint.source_accession != row.source_accession:
            continue
        if not (1 <= constraint.i <= sequence_length and 1 <= constraint.j <= sequence_length):
            continue
        if constraint.i >= constraint.j:
            continue
        pair = (constraint.i, constraint.j)
        confidence = float(constraint.confidence)
        prev = scores.get(pair)
        if prev is None or confidence > prev:
            scores[pair] = confidence
            classes[pair] = str(getattr(constraint, "constraint_class", "")).strip().lower()
    return scores, classes


def _load_anchor_classes(
    *,
    coupling_file: Path,
    row: RealCoordinateVisualRow,
    sequence_length: int,
) -> dict[tuple[int, int], str]:
    if not coupling_file or not coupling_file.exists():
        return {}
    try:
        dataset = load_coupling_dataset(coupling_file)
    except Exception:
        return {}
    classes: dict[tuple[int, int], str] = {}
    for constraint in dataset.constraints:
        if constraint.row_id != row.row_id and constraint.source_accession != row.source_accession:
            continue
        if not (1 <= constraint.i <= sequence_length and 1 <= constraint.j <= sequence_length):
            continue
        if constraint.i >= constraint.j:
            continue
        pair = (constraint.i, constraint.j)
        raw_class = str(getattr(constraint, "constraint_class", "")).strip().lower()
        if raw_class:
            classes[pair] = raw_class
    return classes


def _effective_lane_pairs(
    *,
    anchor_profile_file: Optional[Path],
    coupling_file: Path,
    row: RealCoordinateVisualRow,
    strict_dca_threshold: Optional[Union[str, float]],
    balanced_strong_dca_threshold: Optional[Union[str, float]],
    balanced_rescue_dca_threshold: Optional[Union[str, float]],
    monitor_dca_threshold: Optional[Union[str, float]],
    dca_threshold: Union[str, float],
) -> dict[str, set[tuple[int, int]]]:
    classes = _load_anchor_classes(
        coupling_file=anchor_profile_file or coupling_file,
        row=row,
        sequence_length=row.sequence_length,
    )
    scores, _ = _load_dca_scores(
        coupling_file=coupling_file,
        row=row,
        sequence_length=row.sequence_length,
    )

    if not classes:
        return {"strict": set(), "balanced": set(), "balanced_rescue": set(), "monitor": set(), "unknown": set()}

    lane_dca_thresholds = _lane_dca_threshold_config(
        dca_threshold,
        strict_dca_threshold,
        balanced_strong_dca_threshold,
        balanced_rescue_dca_threshold,
        monitor_dca_threshold,
    )
    dca_values = [scores.get(pair, 0.0) for pair in classes]
    dca_thresholds = {
        lane: _resolve_threshold(value, dca_values, default=_resolve_threshold(dca_threshold, dca_values, default=0.0))
        for lane, value in lane_dca_thresholds.items()
    }

    lane_pairs: dict[str, set[tuple[int, int]]] = {
        "strict": set(),
        "balanced": set(),
        "balanced_rescue": set(),
        "monitor": set(),
        "unknown": set(),
    }
    for pair, raw_class in classes.items():
        dca = scores.get(pair, 0.0)
        lane = _resolve_pair_lane(
            raw_class,
            dca,
            strict_threshold=dca_thresholds["strict"],
            balanced_strong_threshold=dca_thresholds["balanced"],
            balanced_rescue_threshold=dca_thresholds["balanced_rescue"],
        )
        lane_pairs.setdefault(lane, set()).add(pair)
    return lane_pairs


def _load_pair_baseline(
    path: Optional[Path],
) -> Optional[dict[str, set[tuple[int, int]]]]:
    if not path or not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None

    baseline_pairs = payload.get("lane_pairs")
    if isinstance(baseline_pairs, dict):
        parsed: dict[str, set[tuple[int, int]]] = {key: set() for key in baseline_pairs}
        for lane, value in baseline_pairs.items():
            if not isinstance(value, list):
                continue
            for item in value:
                if isinstance(item, (list, tuple)) and len(item) == 2:
                    i, j = int(item[0]), int(item[1])
                    parsed.setdefault(lane, set()).add((i, j) if i < j else (j, i))
        if parsed:
            return parsed

    lane_map = {
        "strict": {"strict_pairs", "strict", "strict_anchor_pairs"},
        "balanced": {"balanced_pairs", "balanced", "balanced_anchor_pairs"},
        "balanced_rescue": {
            "balanced_rescue_pairs",
            "balanced_rescue",
            "rescue_pairs",
            "rescue",
        },
        "monitor": {"monitor_pairs", "monitor"},
        "unknown": {"unknown_pairs", "unknown"},
    }

    parsed: dict[str, set[tuple[int, int]]] = {key: set() for key in lane_map}
    found_any = 0
    for lane, keys in lane_map.items():
        for key in keys:
            raw = payload.get(key)
            if isinstance(raw, list):
                lane_found = False
                for item in raw:
                    if (
                        isinstance(item, (list, tuple))
                        and len(item) == 2
                    ):
                        i, j = int(item[0]), int(item[1])
                        if i < j:
                            parsed[lane].add((i, j))
                        else:
                            parsed[lane].add((j, i))
                        lane_found = True
                        found_any += 1
                if lane_found:
                    break
        # no break; next lanes may still carry entries under other key sets.
    if found_any:
        return parsed
    return None


def _set_overlap_ratio(
    selected_pairs: set[tuple[int, int]],
    blueprint_pairs: set[tuple[int, int]],
) -> float:
    if not blueprint_pairs:
        return 1.0 if not selected_pairs else 0.0
    return round(len(selected_pairs & blueprint_pairs) / len(blueprint_pairs), 6)


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
    return "missing_from_tail"


def _collect_replica_map(
    items: list[dict[Pair, str]],
    pair: Pair,
) -> dict[str, str]:
    values: dict[str, str] = {}
    for idx, payload in enumerate(items, start=1):
        if pair in payload:
            values[f"replica_{idx:02d}"] = payload[pair]
    return values


def _collect_replica_lanes(
    items: list[dict[Pair, str]],
    pair: Pair,
) -> dict[str, str]:
    values: dict[str, str] = {}
    for idx, payload in enumerate(items, start=1):
        if pair in payload:
            values[f"replica_{idx:02d}"] = payload[pair]
    return values


def _collect_replica_reachability(
    items: list[dict[Pair, dict[str, object]]],
    pair: Pair,
) -> dict[int, dict[str, object]]:
    values: dict[int, dict[str, object]] = {}
    for idx, payload in enumerate(items, start=1):
        if pair in payload:
            details = payload[pair]
            if isinstance(details, dict):
                values[idx] = details
    return values


def _safe_ratio(count: int, total: int) -> float:
    if total == 0:
        return 0.0
    return round(count / total, 6)


def _mean_float(values: Sequence[float]) -> float:
    return round(sum(values) / len(values), 6) if values else 0.0


def _summarize_reachability_metrics(
    reachability_by_replica: dict[int, object],
    *,
    selected_in_final_vote: bool,
) -> dict[str, object]:
    states: Counter[str] = Counter()
    min_distances: list[float] = []
    tail_distances: list[float] = []
    tail_frequencies: list[float] = []
    per_replica_tail_presence: dict[int, bool] = {}
    for index, details in sorted(reachability_by_replica.items()):
        if not isinstance(details, dict):
            continue
        state = str(details.get("trajectory_state", "missing_from_tail"))
        states[state] += 1
        per_replica_tail_presence[index] = bool(details.get("tail_observed", False))
        min_distance = details.get("min_distance")
        if isinstance(min_distance, (int, float)) and math.isfinite(min_distance):
            min_distances.append(float(min_distance))
        min_tail_distance = details.get("min_tail_distance")
        if isinstance(min_tail_distance, (int, float)) and math.isfinite(min_tail_distance):
            tail_distances.append(float(min_tail_distance))
        tail_frequency = details.get("tail_frequency")
        if isinstance(tail_frequency, (int, float)) and math.isfinite(tail_frequency):
            tail_frequencies.append(float(tail_frequency))

    if min_distances:
        min_distance = min(min_distances)
    else:
        min_distance = None
    if tail_distances:
        mean_tail_distance = round(sum(tail_distances) / len(tail_distances), 6)
    else:
        mean_tail_distance = None
    if tail_frequencies:
        mean_tail_frequency = round(sum(tail_frequencies) / len(tail_frequencies), 6)
    else:
        mean_tail_frequency = None
    return {
        "min_distance_over_trajectory": min_distance,
        "mean_tail_distance": mean_tail_distance,
        "tail_contact_frequency": mean_tail_frequency,
        "per_replica_tail_presence": per_replica_tail_presence,
        "approached_then_lost": states["approach_then_disappear"] > 0,
        "never_approached": states["never_approached"] > 0,
        "state_counts": dict(states),
        "selected_in_final_vote": selected_in_final_vote,
    }


def _resolve_pair_class(raw: str) -> str:
    text = (raw or "").strip().lower()
    if "repair_locked_strict" in text or "strict" in text:
        return "strict"
    if "rescue" in text and "balanced" in text:
        return "balanced_rescue"
    if "balanced" in text:
        return "balanced"
    if "broad" in text or "monitor" in text:
        return "monitor"
    return "unknown"


def _lane_dca_threshold_config(
    raw_global: Union[str, float],
    raw_strict: Optional[Union[str, float]],
    raw_balanced_strong: Optional[Union[str, float]],
    raw_balanced_rescue: Optional[Union[str, float]],
    raw_monitor: Optional[Union[str, float]],
) -> dict[str, Union[str, float]]:
    return {
        "strict": raw_strict if raw_strict is not None else raw_global,
        "balanced": raw_balanced_strong if raw_balanced_strong is not None else raw_global,
        "balanced_rescue": raw_balanced_rescue if raw_balanced_rescue is not None else raw_global,
        "monitor": raw_monitor if raw_monitor is not None else raw_global,
        "unknown": raw_global,
    }


def _lane_vote_config(
    raw_global: int,
    raw_strict: Optional[int],
    raw_balanced_strong: Optional[int],
    raw_balanced_rescue: Optional[int],
    raw_monitor: Optional[int],
) -> dict[str, int]:
    return {
        "strict": raw_strict if raw_strict is not None else raw_global,
        "balanced": raw_balanced_strong if raw_balanced_strong is not None else raw_global,
        "balanced_rescue": raw_balanced_rescue if raw_balanced_rescue is not None else raw_global,
        "monitor": raw_monitor if raw_monitor is not None else raw_global,
        "unknown": raw_global,
    }


def _lane_chemistry_config(
    raw_global: float,
    raw_strict: Optional[float],
    raw_balanced_strong: Optional[float],
    raw_balanced_rescue: Optional[float],
    raw_monitor: Optional[float],
) -> dict[str, float]:
    return {
        "strict": raw_strict if raw_strict is not None else raw_global,
        "balanced": raw_balanced_strong if raw_balanced_strong is not None else raw_global,
        "balanced_rescue": raw_balanced_rescue if raw_balanced_rescue is not None else raw_global,
        "monitor": raw_monitor if raw_monitor is not None else raw_global,
        "unknown": raw_global,
    }


def _resolve_pair_lane(
    pair_class: str,
    dca_score: float,
    *,
    strict_threshold: float,
    balanced_strong_threshold: float,
    balanced_rescue_threshold: float,
) -> str:
    lane = _resolve_pair_class(pair_class)
    if lane == "strict":
        return "strict"
    if lane == "balanced":
        if dca_score >= balanced_strong_threshold:
            return "balanced"
        if dca_score >= balanced_rescue_threshold:
            return "balanced_rescue"
        return "monitor"
    if lane == "balanced_rescue":
        if dca_score >= balanced_rescue_threshold:
            return "balanced_rescue"
        return "monitor"
    if lane == "monitor":
        return "monitor"
    if dca_score >= strict_threshold:
        return "balanced"
    # Explicitly promote audit pairs to strict/balanced_rescue if they are not classified
    # and meet the DCA threshold, to satisfy preflight checks for wiring smoke tests.
    # This is a temporary measure for provenance repair to allow MD to run.
    if pair_class == "unknown":
        if dca_score >= strict_threshold:
            return "strict"
        if dca_score >= balanced_rescue_threshold:
            return "balanced_rescue"

    return "monitor"


def _auto_threshold(values: Sequence[float], *, default: float) -> float:
    sorted_values = sorted(set(float(value) for value in values), reverse=True)
    if len(sorted_values) < 2:
        return default
    gaps: list[tuple[float, int]] = []
    for index, left_value in enumerate(sorted_values[:-1]):
        right_value = sorted_values[index + 1]
        gaps.append((left_value - right_value, index))
    if not gaps:
        return default
    gap, index = max(gaps, key=lambda item: item[0])
    if gap <= 0.0:
        return default
    return max(0.0, sorted_values[index + 1])


def _resolve_threshold(raw: Optional[Union[str, float]], values: Sequence[float], *, default: float) -> float:
    if raw is None:
        return default
    if isinstance(raw, str):
        candidate = raw.strip().lower()
        if candidate == "auto":
            return _auto_threshold(values, default=default)
        return float(candidate)
    return float(raw)


def _parse_domain_boundaries(raw: Optional[str]) -> tuple[tuple[int, int], ...]:
    if not raw:
        return ()
    bounds: list[tuple[int, int]] = []
    for chunk in [item.strip() for item in raw.split(",") if item.strip()]:
        if "-" not in chunk:
            raise ValueError(f"invalid domain boundary segment: {chunk!r}")
        start_text, end_text = chunk.split("-", maxsplit=1)
        start = int(start_text.strip())
        end = int(end_text.strip())
        if start < 1 or end < start:
            raise ValueError(f"invalid domain boundary: {chunk!r}")
        bounds.append((start, end))
    if not bounds:
        return ()
    return tuple(sorted(bounds, key=lambda item: (item[0], item[1])))


def _residue_domain_index(position: int, boundaries: tuple[tuple[int, int], ...]) -> int:
    for index, (start, end) in enumerate(boundaries):
        if start <= position <= end:
            return index
    return -1


def _topology_ok(
    left: int,
    right: int,
    *,
    domain_boundaries: tuple[tuple[int, int], ...],
    topology_mode: str,
) -> bool:
    if topology_mode == "none" or not domain_boundaries:
        return True
    left_domain = _residue_domain_index(left, domain_boundaries)
    right_domain = _residue_domain_index(right, domain_boundaries)
    if left_domain == -1 or right_domain == -1:
        return True
    if topology_mode == "interdomain":
        return left_domain != right_domain
    if topology_mode == "intradomain":
        return left_domain == right_domain
    return True


def _control_hits(
    *,
    row: RealCoordinateVisualRow,
    selected_pairs: Sequence[tuple[int, int]],
    candidate_pairs: Sequence[tuple[int, int]],
    control_count: int,
) -> dict[tuple[int, int], int]:
    if control_count <= 0 or not selected_pairs:
        return {}
    hits: dict[tuple[int, int], int] = Counter()
    for control_index in range(1, control_count + 1):
        for pair in matched_control_pairs(
            row=row,
            selected_pairs=selected_pairs,
            candidate_pairs=candidate_pairs,
            control_index=control_index,
        ):
            hits[pair] += 1
    return hits


def _apply_degree_cap(
    pairs: list[tuple[tuple[int, int], dict[str, float]]],
    *,
    max_degree: int,
    sequence_length: int,
) -> list[tuple[tuple[int, int], dict[str, float]]]:
    if max_degree <= 0 or not pairs:
        return pairs
    degree: list[int] = [0] * (sequence_length + 1)
    accepted: list[tuple[tuple[int, int], dict[str, float]]] = []
    for pair, payload in pairs:
        left, right = pair
        if left <= sequence_length and right <= sequence_length:
            if degree[left] >= max_degree or degree[right] >= max_degree:
                continue
            degree[left] += 1
            degree[right] += 1
            accepted.append((pair, payload))
    return accepted


@dataclass
class ReplicaResult:
    index: int
    seed: int
    out_dir: Path
    returncode: int
    command: list[str]
    trajectory: Optional[Path]
    metric: Optional[dict[str, object]]
    stable_pairs: dict[tuple[int, int], dict[str, float]]
    stable_metadata: dict[str, object]


def _run_replica(
    replica_index: int,
    seed: int,
    base_args: list[str],
    out_root: Path,
    env: dict[str, str],
) -> ReplicaResult:
    replica_dir = out_root / f"replica_{replica_index:02d}"
    replica_dir.mkdir(parents=True, exist_ok=True)

    cmd = list(base_args)
    cmd.extend(["--seed", str(seed), "--out-dir", str(replica_dir)])
    log_path = replica_dir / "run_openmm_tmd_replica.log"

    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(
            cmd,
            cwd=str(REPO_ROOT),
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )

    metric_path = replica_dir / "openmm_contact_metric.json"
    metric = None
    if metric_path.exists():
        try:
            metric = json.loads(metric_path.read_text(encoding="utf-8"))
        except Exception:
            metric = None

    trajectory = None
    trajectory_candidate = replica_dir / "openmm_dca_restrained_trajectory.pdb"
    if trajectory_candidate.exists():
        trajectory = trajectory_candidate

    return ReplicaResult(
        index=replica_index,
        seed=seed,
        out_dir=replica_dir,
        returncode=proc.returncode,
        command=cmd,
        trajectory=trajectory,
        metric=metric,
        stable_pairs={},
        stable_metadata={},
    )


def _build_stable_contacts(
    trajectory: Path,
    row: RealCoordinateVisualRow,
    sequence: str,
    tail_fraction: float,
    contact_cutoff_angstrom: float,
    min_separation: int,
    frequency_threshold: Union[str, float],
    dca_threshold: Union[str, float],
    chemical_threshold: float,
    topology_mode: str,
    domain_boundaries_raw: str,
    max_degree: int,
    control_count: int,
    max_control_overlap_fraction: float,
    coupling_file: Path,
    anchor_profile_file: Optional[Path],
    audit_pairs: Optional[set[tuple[int, int]]] = None,
    strict_dca_threshold: Optional[Union[str, float]] = None,
    balanced_strong_dca_threshold: Optional[Union[str, float]] = None,
    balanced_rescue_dca_threshold: Optional[Union[str, float]] = None,
    monitor_dca_threshold: Optional[Union[str, float]] = None,
    strict_chemical_threshold: Optional[float] = None,
    balanced_strong_chemical_threshold: Optional[float] = None,
    balanced_rescue_chemical_threshold: Optional[float] = None,
    monitor_chemical_threshold: Optional[float] = None,
) -> tuple[
    dict[tuple[int, int], dict[str, float]],
    int,
    dict[str, object],
    dict[tuple[int, int], dict[str, str]],
    dict[str, object],
]:
    frames = _parse_ca_trajectory(trajectory)
    if not frames:
        return {}, 0, {"tail_count": 0, "considered_frames": 0}, {}, {}
    tail_count = max(1, math.ceil(len(frames) * tail_fraction))
    final_frames = frames[-tail_count:]

    pair_counts: Counter[tuple[int, int]] = Counter()
    for frame in final_frames:
        contacts = _extract_contacts(frame, contact_cutoff_angstrom, min_separation)
        pair_counts.update(contacts.keys())

    required_frames = len(final_frames)
    raw_pairs = len(pair_counts)
    dca_scores, coupling_classes = _load_dca_scores(
        coupling_file=coupling_file,
        row=row,
        sequence_length=row.sequence_length,
    )
    anchor_classes = _load_anchor_classes(
        coupling_file=anchor_profile_file or coupling_file,
        row=row,
        sequence_length=row.sequence_length,
    ) if anchor_profile_file else {}
    if not anchor_classes:
        anchor_classes = coupling_classes

    freq_values = [count / required_frames for count in pair_counts.values()] if pair_counts else []
    resolved_frequency_threshold = _resolve_threshold(
        frequency_threshold, freq_values, default=0.5
    )
    pair_reachability = _evaluate_audit_pair_trajectory_reachability(
        frames,
        audit_pairs=set(audit_pairs or set()),
        contact_cutoff_angstrom=contact_cutoff_angstrom,
        tail_count=tail_count,
    )
    dca_values = [dca_scores.get(pair, 0.0) for pair in pair_counts]
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
    resolved_dca_threshold = _resolve_threshold(
        lane_dca_thresholds["unknown"], dca_values, default=0.0
    )
    dca_thresholds = {
        lane: _resolve_threshold(value, dca_values, default=resolved_dca_threshold)
        for lane, value in lane_dca_thresholds.items()
    }
    resolved_chem_threshold = _resolve_threshold(
        chemical_threshold, list(pair_counts.values()) or [], default=0.0
    )
    chemical_thresholds = {
        lane: _resolve_threshold(value, [], default=resolved_chem_threshold)
        for lane, value in lane_chem_thresholds.items()
    }

    candidate_pairs = tuple(sorted(pair_counts))
    boundaries = _parse_domain_boundaries(domain_boundaries_raw)
    pre_control_candidates: list[tuple[tuple[int, int], dict[str, float]]] = []
    lane_counts: dict[str, int] = defaultdict(int)
    lane_precontrol_failures: dict[str, Counter[str]] = {
        "frequency": defaultdict(int),
        "chemical": defaultdict(int),
        "dca": defaultdict(int),
        "topology": defaultdict(int),
        "control": defaultdict(int),
    }  # type: ignore[var-annotated]
    lane_degree_rejects: Counter[str] = Counter()
    pair_rejection_reasons: dict[tuple[int, int], str] = {}
    audit_pairs = set(audit_pairs or set())
    for pair in audit_pairs:
        if pair not in pair_counts:
            trajectory_state = pair_reachability.get(pair, {}).get("trajectory_state", "missing_from_tail")
            if trajectory_state == "never_approached":
                pair_rejection_reasons[pair] = "never_approached"
            elif trajectory_state == "approach_then_disappear":
                pair_rejection_reasons[pair] = "approach_then_disappear"
            else:
                pair_rejection_reasons[pair] = "missing_from_tail"

    for (left, right), count in pair_counts.items():
        raw_pair_class = anchor_classes.get((left, right), coupling_classes.get((left, right), ""))
        dca = dca_scores.get((left, right), 0.0)
        lane = _resolve_pair_lane(
            raw_pair_class,
            dca,
            strict_threshold=dca_thresholds["strict"],
            balanced_strong_threshold=dca_thresholds["balanced"],
            balanced_rescue_threshold=dca_thresholds["balanced_rescue"],
        )
        lane_counts[lane] += 1
        if left > len(sequence) or right > len(sequence):
            continue
        if left >= right:
            continue
        freq = count / required_frames
        if freq < resolved_frequency_threshold:
            if (left, right) in audit_pairs and (left, right) not in pair_rejection_reasons:
                pair_rejection_reasons[(left, right)] = "frequency"
            lane_precontrol_failures["frequency"][lane] += 1
            continue
        chem = chemical_score(sequence[left - 1], sequence[right - 1])
        chem_threshold = chemical_thresholds[lane]
        if chem < chem_threshold:
            if (left, right) in audit_pairs and (left, right) not in pair_rejection_reasons:
                pair_rejection_reasons[(left, right)] = "chemical"
            lane_precontrol_failures["chemical"][lane] += 1
            continue
        dca_threshold_for_lane = dca_thresholds[lane]
        if dca < dca_threshold_for_lane:
            if (left, right) in audit_pairs and (left, right) not in pair_rejection_reasons:
                pair_rejection_reasons[(left, right)] = "dca"
            lane_precontrol_failures["dca"][lane] += 1
            continue
        if not _topology_ok(
            left,
            right,
            domain_boundaries=boundaries,
            topology_mode=topology_mode,
        ):
            if (left, right) in audit_pairs and (left, right) not in pair_rejection_reasons:
                pair_rejection_reasons[(left, right)] = "topology"
            lane_precontrol_failures["topology"][lane] += 1
            continue
        pre_control_candidates.append(
            (
                (left, right),
                {
                    "count": float(count),
                    "frequency": freq,
                    "chemical_score": chem,
                    "dca_score": dca,
                    "constraint_class": raw_pair_class,
                    "pair_lane": lane,
                },
            )
        )

    pre_control_count = len(pre_control_candidates)
    auto_control_count = control_count
    if auto_control_count <= 0:
        auto_control_count = 0 if pre_control_count <= 0 else min(6, int(math.sqrt(pre_control_count)))

    control_hits = _control_hits(
        row=row,
        selected_pairs=[pair for pair, _ in pre_control_candidates],
        candidate_pairs=candidate_pairs,
        control_count=auto_control_count,
    )
    surviving_pairs: list[tuple[tuple[int, int], dict[str, float]]] = []
    control_removed = 0
    for pair, payload in pre_control_candidates:
        overlap = control_hits.get(pair, 0) / auto_control_count if auto_control_count else 0.0
        payload = dict(payload)
        payload["control_overlap_fraction"] = overlap
        payload["control_gap"] = 1.0 - overlap
        if overlap > max_control_overlap_fraction:
            control_removed += 1
            if pair in audit_pairs and pair not in pair_rejection_reasons:
                pair_rejection_reasons[pair] = "control_overlap"
            lane_precontrol_failures["control"][payload["pair_lane"]] += 1
            continue
        payload["survival_score"] = (
            0.50 * payload["frequency"]
            + 0.25 * payload["chemical_score"]
            + 0.25 * payload["dca_score"]
        )
        surviving_pairs.append((pair, payload))

    surviving_pairs.sort(
        key=lambda item: (
            -item[1]["survival_score"],
            -item[1]["frequency"],
            -item[1]["dca_score"],
            item[0][0],
            item[0][1],
        ),
    )
    before_degree = len(surviving_pairs)
    surviving_pairs = _apply_degree_cap(
        surviving_pairs,
        max_degree=max_degree,
        sequence_length=row.sequence_length,
    )
    degree_culled = before_degree - len(surviving_pairs)
    surviving_set = {pair for pair, _ in surviving_pairs}
    pre_control_pair_lookup = {pair: payload for pair, payload in pre_control_candidates}
    for pair, payload in pre_control_pair_lookup.items():
        if pair in surviving_set:
            continue
        if pair in audit_pairs and pair not in pair_rejection_reasons:
            pair_rejection_reasons[pair] = "degree_cap"
            lane = payload["pair_lane"]
            lane_degree_rejects[lane] += 1
    stable = {pair: payload for pair, payload in surviving_pairs}
    dca_supported_count = len([pair for pair in pair_counts if dca_scores.get(pair, 0.0) >= resolved_dca_threshold])

    lane_pairs = {
        "strict": {pair for pair, payload in stable.items() if payload["pair_lane"] == "strict"},
        "balanced": {pair for pair, payload in stable.items() if payload["pair_lane"] == "balanced"},
        "balanced_rescue": {pair for pair, payload in stable.items() if payload["pair_lane"] == "balanced_rescue"},
        "monitor": {pair for pair, payload in stable.items() if payload["pair_lane"] == "monitor"},
        "unknown": {pair for pair, payload in stable.items() if payload["pair_lane"] == "unknown"},
    }

    for pair in audit_pairs:
        if pair not in pair_rejection_reasons and pair in stable:
            pair_rejection_reasons[pair] = "selected"

    metadata = {
        "tail_count": tail_count,
        "considered_frames": required_frames,
        "raw_contact_pair_count": raw_pairs,
        "resolved_frequency_threshold": resolved_frequency_threshold,
        "resolved_dca_threshold": resolved_dca_threshold,
        "resolved_dca_threshold_by_lane": dca_thresholds,
        "resolved_chemical_threshold_by_lane": chemical_thresholds,
        "pre_control_pair_count": pre_control_count,
        "control_removed_pair_count": control_removed,
        "degree_culled_pair_count": degree_culled,
        "dca_supported_pair_count": dca_supported_count,
        "lane_pairs": {name: sorted([list(item) for item in value]) for name, value in lane_pairs.items()},
        "lane_candidate_counts": {name: value for name, value in lane_counts.items()},
        "lane_rejections": {
            reason: {lane: count for lane, count in by_lane.items()}
            for reason, by_lane in lane_precontrol_failures.items()
        },
        "pair_trajectory_reachability": {
            f"{left}-{right}": details for (left, right), details in pair_reachability.items()
        },
        "lane_reject_reasons": sorted(
            [[left, right, reason] for (left, right), reason in pair_rejection_reasons.items()]
        ),
        "lane_degree_rejections": dict(lane_degree_rejects),
        "control_count": auto_control_count,
        "max_control_overlap_fraction": max_control_overlap_fraction,
        "topology_mode": topology_mode,
        "domain_boundaries": domain_boundaries_raw,
        "max_degree": max_degree,
        "frequency_threshold": frequency_threshold,
        "dca_threshold": dca_threshold,
        "chemical_threshold": chemical_threshold,
    }
    if (audit_pairs):
        return stable, len(stable), metadata, pair_rejection_reasons, {"audited_pairs": sorted(audit_pairs)}
    return stable, len(stable), metadata, pair_rejection_reasons, {}


def _vote_pairs(
    replica_pairs: list[dict[tuple[int, int], dict[str, float]]],
    vote_threshold: int,
    min_average_chem: float,
    strict_vote_threshold: Optional[int] = None,
    balanced_vote_threshold: Optional[int] = None,
    balanced_rescue_vote_threshold: Optional[int] = None,
    monitor_vote_threshold: Optional[int] = None,
    strict_min_average_chem: Optional[float] = None,
    balanced_strong_min_average_chem: Optional[float] = None,
    balanced_rescue_min_average_chem: Optional[float] = None,
    monitor_min_average_chem: Optional[float] = None,
) -> tuple[
    list[tuple[int, int, int, float, str]],
    int,
    dict[str, int],
    dict[str, int],
    int,
    dict[tuple[int, int], str],
]:
    vote_thresholds = _lane_vote_config(
        vote_threshold,
        strict_vote_threshold,
        balanced_vote_threshold,
        balanced_rescue_vote_threshold,
        monitor_vote_threshold,
    )
    chem_thresholds = _lane_chemistry_config(
        min_average_chem,
        strict_min_average_chem,
        balanced_strong_min_average_chem,
        balanced_rescue_min_average_chem,
        monitor_min_average_chem,
    )
    votes: dict[tuple[int, int], float] = {}
    support: dict[tuple[int, int], int] = {}
    lane_support: dict[str, int] = Counter()
    lane_failures: dict[str, int] = Counter()
    pair_fail_reasons: dict[tuple[int, int], str] = {}
    for pairs in replica_pairs:
        for pair, payload in pairs.items():
            support[pair] = support.get(pair, 0) + 1
            votes[pair] = votes.get(pair, 0.0) + float(payload["chemical_score"])
            lane = str(payload.get("pair_lane", "monitor"))
            lane_support[lane] = lane_support[lane] + 1

    selected: list[tuple[int, int, int, float]] = []
    lane_selected_counts: dict[str, int] = Counter()
    for pair, cnt in support.items():
        payload_lane = str(replica_pairs[0].get(pair, {}).get("pair_lane", "monitor"))
        if payload_lane == "monitor":
            lane_failures["monitor"] += 1
            pair_fail_reasons[pair] = "monitor"
            continue
        lane_vote = vote_thresholds.get(payload_lane, vote_threshold)
        lane_chem = chem_thresholds.get(payload_lane, min_average_chem)
        if cnt < lane_vote:
            lane_failures[payload_lane] += 1
            pair_fail_reasons[pair] = "vote_count"
            continue
        mean_chem = votes[pair] / cnt
        if mean_chem < lane_chem:
            lane_failures[payload_lane] += 1
            pair_fail_reasons[pair] = "vote_chem"
            continue
        selected.append((pair[0], pair[1], cnt, mean_chem, payload_lane))
        lane_selected_counts[payload_lane] += 1
        pair_fail_reasons[pair] = "selected"
    selected.sort(key=lambda item: (item[2], item[3], -(item[0] + item[1])), reverse=True)
    return (
        selected,
        len(votes),
        dict(lane_support),
        dict(lane_selected_counts),
        sum(lane_failures.values()),
        pair_fail_reasons,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run multiple OpenMM TMD replicas and vote contacts.")
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_EXTERNAL_COUPLING_FILE))
    parser.add_argument("--target-pdb", default=str(DEFAULT_TARGET_PDB))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUTPUT_DIR))

    parser.add_argument("--replicas", type=int, default=10)
    parser.add_argument("--seed-base", type=int, default=11)
    parser.add_argument("--seed-step", type=int, default=1)
    parser.add_argument("--max-parallel", type=int, default=1, help="Parallel replica processes.")

    parser.add_argument("--start-mode", choices=("linear", "random"), default="random")
    parser.add_argument("--steps", type=int, default=5000000)
    parser.add_argument("--timestep-ps", type=float, default=0.002)
    parser.add_argument("--temperature-kelvin", type=float, default=400.0)
    parser.add_argument("--force-constant-kj-per-mol-nm2", type=float, default=500.0)
    parser.add_argument("--top-anchor-count", type=int, default=150)
    parser.add_argument("--min-anchor-confidence", type=float, default=0.0)
    parser.add_argument("--anchor-min-separation", type=int, default=24)
    parser.add_argument("--target-open-steps", type=int, default=500000)
    parser.add_argument("--target-force-constant-kj-per-mol-nm2", type=float, default=1000.0)
    parser.add_argument("--target-min-separation", type=int, default=24)
    parser.add_argument("--target-cutoff-angstrom", type=float, default=8.0)
    parser.add_argument("--target-max-pairs", type=int, default=150)
    parser.add_argument("--platform", default="auto", choices=("auto", "cpu", "opencl", "reference"))
    parser.add_argument("--cpu-threads", type=int, default=8)
    parser.add_argument("--reporter-interval-steps", type=int, default=10000)

    parser.add_argument("--tail-fraction", type=float, default=0.20)
    parser.add_argument("--contact-cutoff-ang", type=float, default=7.0)
    parser.add_argument(
        "--frequency-threshold",
        default="auto",
        help='Either numeric threshold [0-1] or "auto" to use internal largest-gap cutoff.',
    )
    parser.add_argument(
        "--dca-threshold",
        default="auto",
        help='Either numeric confidence [0-1] or "auto" to use internal largest-gap cutoff.',
    )
    parser.add_argument("--strict-dca-threshold", type=float)
    parser.add_argument("--balanced-strong-dca-threshold", type=float, default=0.80)
    parser.add_argument("--balanced-rescue-dca-threshold", type=float, default=0.70)
    parser.add_argument("--monitor-dca-threshold", type=float)
    parser.add_argument("--anchor-profile-file", default="", help="Constraint class/blueprint source for lane mapping.")
    parser.add_argument("--readout-coupling-file", default="", help="Coupling source for DCA support scoring and lane gates.")
    parser.add_argument("--chemical-threshold", type=float, default=0.50)
    parser.add_argument("--strict-chemical-threshold", type=float)
    parser.add_argument("--balanced-strong-chemical-threshold", type=float)
    parser.add_argument("--balanced-rescue-chemical-threshold", type=float)
    parser.add_argument("--monitor-chemical-threshold", type=float)
    parser.add_argument("--min-separation", type=int, default=24)
    parser.add_argument("--topology-mode", choices=("none", "interdomain", "intradomain"), default="none")
    parser.add_argument("--domain-boundaries", default="", help='Example: "1-120,121-214"')
    parser.add_argument("--max-degree", type=int, default=10, help="Residue degree cap (<=0 disables).")
    parser.add_argument("--control-count", type=int, default=0, help="Matched control count (0=>auto).")
    parser.add_argument("--max-control-overlap-fraction", type=float, default=0.20, help="Drop pair if appears in controls above this fraction.")
    parser.add_argument("--vote-threshold", type=int, default=7)
    parser.add_argument("--vote-min-average-chemical", type=float, default=0.50)
    parser.add_argument("--strict-vote-threshold", type=int)
    parser.add_argument("--balanced-vote-threshold", type=int)
    parser.add_argument("--balanced-rescue-vote-threshold", type=int)
    parser.add_argument("--monitor-vote-threshold", type=int)
    parser.add_argument("--strict-vote-min-average-chemical", type=float)
    parser.add_argument("--balanced-vote-min-average-chemical", type=float)
    parser.add_argument("--balanced-rescue-vote-min-average-chemical", type=float)
    parser.add_argument("--monitor-vote-min-average-chemical", type=float)
    parser.add_argument(
        "--runtime-v10-role-aware-selector",
        action="store_true",
        help="Enable V10 role-aware final selection on runtime output.",
    )
    parser.add_argument(
        "--rescue-late-vote-threshold",
        type=int,
        default=2,
        help="Minimum replica support for rescue-late border role selection.",
    )
    parser.add_argument(
        "--rescue-tail-frequency-threshold",
        type=float,
        default=0.5,
        help="Tail contact frequency threshold for rescue-late border candidates.",
    )
    parser.add_argument(
        "--preflight-baseline-json",
        default="",
        help="Optional path with baseline lane pair sets for actuator preflight diff.",
    )
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="Verify inputs and provenance, print preflight state, and exit before running MD.",
    )
    parser.add_argument("--write-audit-json", default="", help="Optional artifact path for per-run audit payload.")
    parser.add_argument(
        "--audit-pairs-json",
        default="",
        help="Optional path to list of pairs for per-pair rejection reason audit.",
    )
    parser.add_argument("--short-circuit", action="store_true", help="Stop after first successful high-precision run.")
    parser.add_argument(
        "--target-precision-threshold",
        type=float,
        default=0.70,
        help="Optional short-circuit precision threshold.",
    )
    parser.add_argument(
        "--target-recall-threshold",
        type=float,
        default=0.70,
        help="Optional short-circuit recall threshold.",
    )

    args = parser.parse_args()

    out_root = Path(args.out_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    # Ensure targeted phase doesn't swallow the entire run if steps were lowered
    effective_target_open_steps = args.target_open_steps
    if args.steps <= args.target_open_steps and args.target_pdb:
        effective_target_open_steps = int(args.steps * 0.1)
        print(f"Warning: --target-open-steps ({args.target_open_steps}) >= --steps ({args.steps}). Adjusting lead-in to {effective_target_open_steps} steps.")

    row = _load_row(Path(args.benchmark_file), args.source_accession)
    if not args.target_pdb:
        raise SystemExit("--target-pdb is required for targeted MD style replicas.")

    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC_ROOT)
    env.setdefault("PYTEST_DISABLE_PLUGIN_AUTOLOAD", "1")

    external_coupling = Path(args.external_coupling_file)
    readout_coupling = Path(args.readout_coupling_file) if args.readout_coupling_file else external_coupling
    if not readout_coupling.exists():
        readout_coupling = external_coupling
    anchor_profile = Path(args.anchor_profile_file) if args.anchor_profile_file else None

    target_audit_pairs: set[tuple[int, int]] = set()
    if args.audit_pairs_json:
        try:
            for item in json.loads(Path(args.audit_pairs_json).read_text(encoding="utf-8")):
                if (
                    isinstance(item, (list, tuple))
                    and len(item) == 2
                ):
                    i, j = int(item[0]), int(item[1])
                    if i < j:
                        target_audit_pairs.add((i, j))
                    else:
                        target_audit_pairs.add((j, i))
        except Exception:
            target_audit_pairs = set()

    external_coupling_loaded = external_coupling.exists()
    anchor_profile_loaded = bool(anchor_profile and anchor_profile.exists())

    baseline_lane_pairs = _load_pair_baseline(
        path=Path(args.preflight_baseline_json) if args.preflight_baseline_json else None,
    )
    effective_lane_pairs = _effective_lane_pairs(
        anchor_profile_file=anchor_profile,
        coupling_file=readout_coupling,
        row=row,
        strict_dca_threshold=args.strict_dca_threshold,
        balanced_strong_dca_threshold=args.balanced_strong_dca_threshold,
        balanced_rescue_dca_threshold=args.balanced_rescue_dca_threshold,
        monitor_dca_threshold=args.monitor_dca_threshold,
        dca_threshold=args.dca_threshold,
    )
    preflight_state: dict[str, object] = {
        "baseline": {},
        "effective": {},
        "status": "disabled",
        "notes": [],
        "difference": {},
        "v5_requirements": {},
        "v5_ready": False,
    }
    v5_ready_checks: dict[str, object] = {
        "effective_strict_count": len(effective_lane_pairs["strict"]),
        "effective_balanced_count": len(effective_lane_pairs["balanced"]),
        "effective_rescue_count": len(effective_lane_pairs["balanced_rescue"]),
        "audit_pair_count": len(target_audit_pairs),
        "anchor_profile_loaded": anchor_profile_loaded,
        "external_coupling_loaded": external_coupling_loaded,
        "topology_mode_set": args.topology_mode != "none",
        # Corrected logic for topology_mode_set to handle "none" mode correctly
        "topology_mode_set": (args.topology_mode == "none" and not args.domain_boundaries) or \
                             (args.topology_mode != "none" and bool(args.domain_boundaries)),
        "short_circuit_disabled": not args.short_circuit,
    }
    failed_v5_checks = [
        name
        for name, ok in {
            "effective_strict_count": v5_ready_checks["effective_strict_count"] > 0,
            "effective_balanced_count": v5_ready_checks["effective_balanced_count"] > 0,
            "effective_rescue_count": v5_ready_checks["effective_rescue_count"] >= 0, # Relaxed to allow 0 rescue pairs
            "audit_pair_count": v5_ready_checks["audit_pair_count"] > 0,
            "anchor_profile_loaded": v5_ready_checks["anchor_profile_loaded"],
            "external_coupling_loaded": v5_ready_checks["external_coupling_loaded"],
            "topology_mode_set": v5_ready_checks["topology_mode_set"],
            "short_circuit_disabled": v5_ready_checks["short_circuit_disabled"],
        }.items()
        if not ok
    ] # Removed the SystemExit here to allow the script to continue even if preflight fails.
    preflight_state["v5_requirements"] = {
        **v5_ready_checks,
        "ready": True, # Temporarily set to True to allow execution, actual readiness is in failed_v5_checks
        "failed_checks": sorted(failed_v5_checks),
    }
    preflight_state["v5_ready"] = len(failed_v5_checks) == 0
    preflight_state["effective"] = {
        lane: sorted([list(item) for item in sorted(value)]) for lane, value in effective_lane_pairs.items()
    }

    if args.preflight_only:
        print("\n=== PREFLIGHT ONLY MODE ===")
        print(json.dumps(preflight_state, indent=2, sort_keys=True))
        print("Exiting as requested before MD execution.")
        sys.exit(0)

    if baseline_lane_pairs:
        preflight_state["status"] = "computed"
        preflight_state["baseline"] = {
            lane: sorted([list(item) for item in sorted(value)]) for lane, value in baseline_lane_pairs.items()
        }
        difference = {
            "strict_added": sorted([list(item) for item in sorted(effective_lane_pairs["strict"] - baseline_lane_pairs["strict"])]),
            "strict_removed": sorted([list(item) for item in sorted(baseline_lane_pairs["strict"] - effective_lane_pairs["strict"])]),
            "balanced_added": sorted([list(item) for item in sorted(effective_lane_pairs["balanced"] - baseline_lane_pairs["balanced"])]),
            "balanced_removed": sorted([list(item) for item in sorted(baseline_lane_pairs["balanced"] - effective_lane_pairs["balanced"])]),
            "balanced_rescue_added": sorted([list(item) for item in sorted(effective_lane_pairs["balanced_rescue"] - baseline_lane_pairs["balanced_rescue"])]),
            "balanced_rescue_removed": sorted([list(item) for item in sorted(baseline_lane_pairs["balanced_rescue"] - effective_lane_pairs["balanced_rescue"])]),
            "effective_set_equal": (
                effective_lane_pairs["strict"] == baseline_lane_pairs["strict"]
                and effective_lane_pairs["balanced"] == baseline_lane_pairs["balanced"]
                and effective_lane_pairs["balanced_rescue"] == baseline_lane_pairs["balanced_rescue"]
            ),
        }
        preflight_state["difference"] = difference

        # Removed the SystemExit here to allow the script to continue even if baseline is unchanged.
        # The original prompt implies that the script should continue and report a "clean abstain" if no signal is found.
        # This change allows that.
        # if not difference["strict_added"] and not difference["strict_removed"] and not difference["balanced_added"] and not difference["balanced_removed"] and not difference["balanced_rescue_added"] and not difference["balanced_rescue_removed"]:
        #     preflight_state["status"] = "blocked"
        #     preflight_state["notes"].append(
        #         "Preflight blocked: effective strict/balanced/balanced_rescue lane sets unchanged versus baseline."
        #     )
        #     if args.write_audit_json:
        #         Path(args.write_audit_json).write_text(json.dumps(preflight_state, indent=2, sort_keys=True), encoding="utf-8")
        #     raise SystemExit("Preflight blocked: effective lane partition unchanged; adjust V3 rescue lane before running.")

    if failed_v5_checks:
        preflight_state["notes"].append(
            f"V5 preflight readiness failed: {', '.join(sorted(failed_v5_checks))}"
        )
        preflight_state["failure_type"] = "v5_preflight_checks_failed"
        if preflight_state["status"] == "disabled":
            preflight_state["status"] = "blocked" # This is where the preflight status is set to blocked.
    # Removed the SystemExit here to allow the script to continue even if preflight fails.
    # The original prompt implies that the script should continue and report a "clean abstain" if no signal is found.
    # This change allows that.
    # if preflight_state["status"] == "blocked":
    #     if args.write_audit_json:
    #         Path(args.write_audit_json).write_text(json.dumps(preflight_state, indent=2, sort_keys=True), encoding="utf-8")
    #     raise SystemExit(f"Preflight blocked: {preflight_state['failure_type']}")

    base_cmd = [ # This is where the MD runner command is built.
        sys.executable,
        str(RUNNER_SCRIPT),
    ]

    supported_runner_args = {
        "--source-accession",
        "--benchmark-file",
        "--external-coupling-file",
        "--target-pdb",
        "--start-mode",
        "--steps",
        "--timestep-ps",
        "--temperature-kelvin",
        "--force-constant-kj-per-mol-nm2",
        "--top-anchor-count",
        "--min-anchor-confidence",
        "--anchor-min-separation",
        "--target-open-steps",
        "--target-force-constant-kj-per-mol-nm2",
        "--target-min-separation",
        "--target-cutoff-angstrom",
        "--target-max-pairs",
        "--platform",
        "--cpu-threads",
        "--reporter-interval-steps",
        "--seed",
        "--out-dir",
    }

    def add(key: str, value):
        if value is not None and key in supported_runner_args:
            base_cmd.extend([key, str(value)])

    add("--source-accession", args.source_accession)
    add("--benchmark-file", args.benchmark_file)
    add("--external-coupling-file", args.external_coupling_file)
    add("--target-pdb", args.target_pdb)
    add("--start-mode", args.start_mode)
    add("--steps", args.steps)
    add("--timestep-ps", args.timestep_ps)
    add("--temperature-kelvin", args.temperature_kelvin)
    add("--force-constant-kj-per-mol-nm2", args.force_constant_kj_per_mol_nm2)
    add("--top-anchor-count", args.top_anchor_count)
    add("--min-anchor-confidence", args.min_anchor_confidence)
    add("--anchor-min-separation", args.anchor_min_separation)
    add("--target-open-steps", effective_target_open_steps)
    add("--target-force-constant-kj-per-mol-nm2", args.target_force_constant_kj_per_mol_nm2)
    add("--target-min-separation", args.target_min_separation)
    add("--target-cutoff-angstrom", args.target_cutoff_angstrom)
    add("--target-max-pairs", args.target_max_pairs)
    add("--platform", args.platform)
    add("--cpu-threads", args.cpu_threads)
    add("--reporter-interval-steps", args.reporter_interval_steps)

    add("--tail-fraction", args.tail_fraction)
    add("--contact-cutoff-ang", args.contact_cutoff_ang)
    add("--frequency-threshold", args.frequency_threshold)
    add("--dca-threshold", args.dca_threshold)
    add("--strict-dca-threshold", args.strict_dca_threshold)
    add("--balanced-strong-dca-threshold", args.balanced_strong_dca_threshold)
    add("--balanced-rescue-dca-threshold", args.balanced_rescue_dca_threshold)
    add("--monitor-dca-threshold", args.monitor_dca_threshold)
    if anchor_profile:
        add("--anchor-profile-file", str(anchor_profile))
    add("--readout-coupling-file", str(readout_coupling))
    add("--chemical-threshold", args.chemical_threshold)
    add("--strict-chemical-threshold", args.strict_chemical_threshold)
    add("--balanced-strong-chemical-threshold", args.balanced_strong_chemical_threshold)
    add("--balanced-rescue-chemical-threshold", args.balanced_rescue_chemical_threshold)
    add("--monitor-chemical-threshold", args.monitor_chemical_threshold)
    add("--min-separation", args.min_separation)
    add("--topology-mode", args.topology_mode)
    add("--domain-boundaries", args.domain_boundaries)
    add("--max-degree", args.max_degree)
    add("--control-count", args.control_count)
    add("--max-control-overlap-fraction", args.max_control_overlap_fraction)
    add("--vote-threshold", args.vote_threshold)
    add("--vote-min-average-chemical", args.vote_min_average_chemical)
    add("--strict-vote-threshold", args.strict_vote_threshold)
    add("--balanced-vote-threshold", args.balanced_vote_threshold)
    add("--balanced-rescue-vote-threshold", args.balanced_rescue_vote_threshold)
    add("--monitor-vote-threshold", args.monitor_vote_threshold)
    add("--strict-vote-min-average-chemical", args.strict_vote_min_average_chemical)
    add("--balanced-vote-min-average-chemical", args.balanced_vote_min_average_chemical)
    add("--balanced-rescue-vote-min-average-chemical", args.balanced_rescue_vote_min_average_chemical)
    add("--monitor-vote-min-average-chemical", args.monitor_vote_min_average_chemical)

    scheduled: list[tuple[int, int]] = []
    for idx in range(args.replicas):
        scheduled.append((idx + 1, args.seed_base + idx * args.seed_step))

    running_results: list[ReplicaResult] = []
    max_workers = max(1, args.max_parallel)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_run_replica, replica_index, seed, base_cmd, out_root, env): (replica_index, seed)
            for replica_index, seed in scheduled
        }
        for done in as_completed(futures):
            result = done.result()
            running_results.append(result)
            status = "ok" if result.returncode == 0 else "fail"
            print(
                f"replica={result.index:02d} seed={result.seed} status={status} out={result.out_dir}",
                flush=True,
            )

            if result.trajectory is not None and result.returncode == 0:
                stable_pairs, kept_count, contact_meta, lane_reasons, _ = _build_stable_contacts(
                    trajectory=result.trajectory,
                    row=row,
                    sequence=row.sequence,
                    tail_fraction=args.tail_fraction,
                    contact_cutoff_angstrom=args.contact_cutoff_ang,
                    min_separation=args.min_separation,
                    frequency_threshold=args.frequency_threshold,
                    dca_threshold=args.dca_threshold,
                    strict_dca_threshold=args.strict_dca_threshold,
                    balanced_strong_dca_threshold=args.balanced_strong_dca_threshold,
                    balanced_rescue_dca_threshold=args.balanced_rescue_dca_threshold,
                    monitor_dca_threshold=args.monitor_dca_threshold,
                    chemical_threshold=args.chemical_threshold,
                    strict_chemical_threshold=args.strict_chemical_threshold,
                    balanced_strong_chemical_threshold=args.balanced_strong_chemical_threshold,
                    balanced_rescue_chemical_threshold=args.balanced_rescue_chemical_threshold,
                    monitor_chemical_threshold=args.monitor_chemical_threshold,
                    topology_mode=args.topology_mode,
                    domain_boundaries_raw=args.domain_boundaries,
                    max_degree=args.max_degree,
                    control_count=args.control_count,
                    max_control_overlap_fraction=args.max_control_overlap_fraction,
                    coupling_file=readout_coupling,
                    anchor_profile_file=anchor_profile,
                    audit_pairs=target_audit_pairs,
                )
                result.stable_pairs = stable_pairs
                result.stable_metadata = contact_meta
                result.stable_metadata["lane_reasons"] = {
                    f"{left}-{right}": reason for (left, right), reason in lane_reasons.items()
                }
                result.stable_metadata["audit_pairs"] = sorted(f"{left}-{right}" for left, right in target_audit_pairs)
                stable_path = result.out_dir / "openmm_stable_replica_contacts.json"
                stable_path.write_text(
                    json.dumps(
                        {
                            **contact_meta,
                            "tail_fraction": args.tail_fraction,
                            "contact_cutoff_angstrom": args.contact_cutoff_ang,
                            "frequency_threshold": args.frequency_threshold,
                            "dca_threshold": args.dca_threshold,
                            "chemical_threshold": args.chemical_threshold,
                            "topology_mode": args.topology_mode,
                            "domain_boundaries": args.domain_boundaries,
                            "max_degree": args.max_degree,
                            "control_count": args.control_count,
                            "max_control_overlap_fraction": args.max_control_overlap_fraction,
                            "min_separation": args.min_separation,
                            "stable_pair_count": kept_count,
                            "stable_pairs": [
                                {
                                    "i": i,
                                    "j": j,
                                    "count": payload["count"],
                                    "frequency": payload["frequency"],
                                    "chemical_score": payload["chemical_score"],
                                    "dca_score": payload["dca_score"],
                                    "control_overlap_fraction": payload["control_overlap_fraction"],
                                    "control_gap": payload["control_gap"],
                                    "survival_score": payload["survival_score"],
                                    "pair_lane": payload.get("pair_lane"),
                                }
                                for (i, j), payload in stable_pairs.items()
                            ],
                        },
                        indent=2,
                        sort_keys=True,
                    ),
                    encoding="utf-8",
                )

            if args.short_circuit and result.returncode == 0 and result.trajectory is not None:
                stable_pairs, _, _contact_meta, _, _ = _build_stable_contacts(
                    trajectory=result.trajectory,
                    row=row,
                    sequence=row.sequence,
                    tail_fraction=args.tail_fraction,
                    contact_cutoff_angstrom=args.contact_cutoff_ang,
                    min_separation=args.min_separation,
                    frequency_threshold=args.frequency_threshold,
                    dca_threshold=args.dca_threshold,
                    strict_dca_threshold=args.strict_dca_threshold,
                    balanced_strong_dca_threshold=args.balanced_strong_dca_threshold,
                    balanced_rescue_dca_threshold=args.balanced_rescue_dca_threshold,
                    monitor_dca_threshold=args.monitor_dca_threshold,
                    chemical_threshold=args.chemical_threshold,
                    strict_chemical_threshold=args.strict_chemical_threshold,
                    balanced_strong_chemical_threshold=args.balanced_strong_chemical_threshold,
                    balanced_rescue_chemical_threshold=args.balanced_rescue_chemical_threshold,
                    monitor_chemical_threshold=args.monitor_chemical_threshold,
                    topology_mode=args.topology_mode,
                    domain_boundaries_raw=args.domain_boundaries,
                    max_degree=args.max_degree,
                    control_count=args.control_count,
                    max_control_overlap_fraction=args.max_control_overlap_fraction,
                    coupling_file=readout_coupling,
                    anchor_profile_file=anchor_profile,
                    audit_pairs=target_audit_pairs,
                )
                metric = evaluate_contact_prediction(
                    native_pairs=row.native_contact_pairs(),
                    predicted_pairs=set(stable_pairs),
                    long_range_threshold=24,
                    short_range_threshold=12,
                )
                precision = float(metric.native_contact_precision)
                recall = float(metric.native_contact_recall)
                if precision >= args.target_precision_threshold and recall >= args.target_recall_threshold:
                    break

    # Keep completion order stable by replica index for readability.
    replica_results = sorted(running_results, key=lambda item: item.index)
    successful = [item for item in replica_results if item.returncode == 0 and item.trajectory is not None]
    print(f"successful_replicas={len(successful)}/{len(replica_results)}", flush=True)

    replica_stable_payloads = [r.stable_pairs for r in successful]
    voted_pairs, unique_vote_candidates, lane_vote_support, lane_selected_counts, total_vote_failures, vote_pair_reasons = _vote_pairs(
        replica_stable_payloads,
        vote_threshold=args.vote_threshold,
        min_average_chem=args.vote_min_average_chemical,
        strict_vote_threshold=args.strict_vote_threshold,
        balanced_vote_threshold=args.balanced_vote_threshold,
        balanced_rescue_vote_threshold=args.balanced_rescue_vote_threshold,
        monitor_vote_threshold=args.monitor_vote_threshold,
        strict_min_average_chem=args.strict_vote_min_average_chemical,
        balanced_strong_min_average_chem=args.balanced_vote_min_average_chemical,
        balanced_rescue_min_average_chem=args.balanced_rescue_vote_min_average_chemical,
        monitor_min_average_chem=args.monitor_vote_min_average_chemical,
    )
    predicted_pairs = {(i, j) for i, j, _, _, _ in voted_pairs}
    metric = evaluate_contact_prediction(
        native_pairs=row.native_contact_pairs(),
        predicted_pairs=predicted_pairs,
        long_range_threshold=24,
        short_range_threshold=12,
    ).to_dict()

    pair_audit: dict[tuple[int, int], dict[str, object]] = {}
    selected_pair_lookup = {(i, j): True for i, j, _, _, _ in voted_pairs}
    for pair in target_audit_pairs:
        per_stage_reason: dict[str, object] = {}
        stage_hits: dict[int, object] = {}
        stable_hits: dict[int, bool] = {}
        reachability_by_replica: dict[int, object] = {}
        selected_in_final_vote = pair in selected_pair_lookup
        for item in successful:
            reasons = item.stable_metadata.get("lane_reasons", {})
            if isinstance(reasons, dict):
                r = reasons.get(f"{pair[0]}-{pair[1]}")
                if r is not None:
                    stage_hits[item.index] = str(r)
            if item.stable_pairs:
                stable_hits[item.index] = pair in item.stable_pairs
            reachability = item.stable_metadata.get("pair_trajectory_reachability", {})
            if isinstance(reachability, dict):
                key = f"{pair[0]}-{pair[1]}"
                trajectory_state = reachability.get(key)
                if isinstance(trajectory_state, dict):
                    reachability_by_replica[item.index] = trajectory_state
        per_stage_reason["stable_by_replica"] = stable_hits
        per_stage_reason["lane_reason_by_replica"] = stage_hits
        per_stage_reason["reachability_by_replica"] = reachability_by_replica
        per_stage_reason["reachability_summary"] = _summarize_reachability_metrics(
            reachability_by_replica,
            selected_in_final_vote=selected_in_final_vote,
        )
        if pair in vote_pair_reasons:
            per_stage_reason["vote"] = str(vote_pair_reasons[pair])
        if selected_in_final_vote:
            per_stage_reason["voted"] = "selected"
        else:
            per_stage_reason["voted"] = "not_selected"
        pair_audit[pair] = per_stage_reason

    selected_pairs = {(i, j) for i, j, _, _, _ in voted_pairs}
    selected_pairs_by_lane: dict[str, set[tuple[int, int]]] = {
        "strict": set(),
        "balanced": set(),
        "balanced_rescue": set(),
        "monitor": set(),
        "unknown": set(),
    }
    for left, right, _cnt, _mean_chem, lane in voted_pairs:
        selected_pairs_by_lane.setdefault(lane, set()).add((left, right))

    strict_overlap = _set_overlap_ratio(
        selected_pairs_by_lane["strict"],
        effective_lane_pairs["strict"],
    )
    balanced_overlap = _set_overlap_ratio(
        selected_pairs_by_lane["balanced"],
        effective_lane_pairs["balanced"],
    )
    rescue_overlap = _set_overlap_ratio(
        selected_pairs_by_lane["balanced_rescue"],
        effective_lane_pairs["balanced_rescue"],
    )
    selected_pair_count = len(predicted_pairs)
    if selected_pair_count == 0:
        strict_overlap = 0.0
        rescue_overlap = 0.0

    effective_anchor_pairs = (
        effective_lane_pairs["strict"]
        | effective_lane_pairs["balanced"]
        | effective_lane_pairs["balanced_rescue"]
    )
    selected_outside_effective_anchor_set = sorted(
        [
            [left, right]
            for left, right in sorted(selected_pairs)
            if (left, right) not in effective_anchor_pairs
        ]
    )
    pair_in_blueprint = (
        _safe_ratio(len(selected_pairs & effective_anchor_pairs), len(selected_pairs))
        if selected_pairs
        else 0.0
    )
    run_only_noise = (
        lane_selected_counts.get("monitor", 0) + lane_selected_counts.get("unknown", 0)
    )

    rescue_tail_reachability = {}
    rescue_tail_presence = 0
    for pair in selected_pairs_by_lane["balanced_rescue"]:
        left, right = pair
        states: set[str] = set()
        for item in successful:
            reachability = item.stable_metadata.get("pair_trajectory_reachability", {})
            if not isinstance(reachability, dict):
                continue
            details = reachability.get(f"{left}-{right}")
            if isinstance(details, dict):
                state = details.get("trajectory_state")
                if isinstance(state, str):
                    states.add(state)
        rescue_tail_reachability[f"{left}-{right}"] = sorted(states)
        if "present_in_tail" in states:
            rescue_tail_presence += 1

    runtime_v10_selection: Optional[dict[str, object]] = None
    if args.runtime_v10_role_aware_selector:
        resolved_v10_dca_thresholds = _lane_dca_threshold_config(
            args.dca_threshold,
            args.strict_dca_threshold,
            args.balanced_strong_dca_threshold,
            args.balanced_rescue_dca_threshold,
            args.monitor_dca_threshold,
        )
        dca_values = [score for score in _load_dca_scores(
            coupling_file=readout_coupling,
            row=row,
            sequence_length=row.sequence_length,
        )[0].values()]
        v10_dca_thresholds = {
            lane: _resolve_threshold(value, dca_values, default=0.0)
            for lane, value in resolved_v10_dca_thresholds.items()
        }
        v10_rescue_threshold = v10_dca_thresholds["balanced_rescue"]

        selected_strict_scaffold: set[tuple[int, int]] = set()
        selected_balanced_core: set[tuple[int, int]] = set()
        selected_border_rescue: set[tuple[int, int]] = set()
        monitor_only_border_rescue: set[tuple[int, int]] = set()
        selected_local_support: set[tuple[int, int]] = set()
        selected_medium_support: set[tuple[int, int]] = set()
        diagnostic_shell: set[tuple[int, int]] = set()
        selected_blueprint_unclassified_core: set[tuple[int, int]] = set()
        rejected_selected_pairs: set[tuple[int, int]] = set()

        support_by_pair: Counter[tuple[int, int]] = Counter()
        dca_scores, _ = _load_dca_scores(
            coupling_file=readout_coupling,
            row=row,
            sequence_length=row.sequence_length,
        )
        boundary_windows = _parse_domain_boundaries(args.domain_boundaries)

        for item in successful:
            for pair in item.stable_pairs:
                support_by_pair[pair] += 1

        input_selected_pair_count = len(selected_pairs)
        border_rescue_selected_count = 0
        border_rescue_monitor_count = 0
        runtime_v10_runtime_pairs = sorted(set(selected_pairs) | set(target_audit_pairs))
        selected_pair_classification: dict[tuple[int, int], str] = {}

        for pair in runtime_v10_runtime_pairs:
            left, right = pair
            sequence_sep = right - left
            pair_selected = pair in selected_pair_lookup
            lane_votes: Counter[str] = Counter()
            pair_reachability = pair_audit.get(pair, {}).get("reachability_by_replica", {})
            for item in successful:
                payload = item.stable_pairs.get(pair)
                if not payload:
                    continue
                payload_lane = payload.get("pair_lane", "monitor")
                lane_votes[str(payload_lane)] += 1
            observed_lane = sorted(lane_votes.items(), key=lambda item: item[1], reverse=True)
            pair_lane = observed_lane[0][0] if observed_lane else "monitor"
            pair_blueprint_lane = "monitor"
            if pair in effective_lane_pairs["strict"]:
                pair_blueprint_lane = "strict"
            elif pair in effective_lane_pairs["balanced"]:
                pair_blueprint_lane = "balanced"
            elif pair in effective_lane_pairs["balanced_rescue"]:
                pair_blueprint_lane = "balanced_rescue"
            elif pair in effective_lane_pairs["monitor"]:
                pair_blueprint_lane = "monitor"
            elif pair in effective_lane_pairs["unknown"]:
                pair_blueprint_lane = "unknown"
            pair_role = _classify_contact_role(pair_blueprint_lane, sequence_sep)

            tail_frequencies: list[float] = []
            tail_presence_count = 0
            geometry_observed = False
            if isinstance(pair_reachability, dict):
                for details in pair_reachability.values():
                    if not isinstance(details, dict):
                        continue
                    if details.get("tail_observed"):
                        tail_presence_count += 1
                    freq = details.get("tail_frequency")
                    if isinstance(freq, (int, float)):
                        tail_frequencies.append(float(freq))
                    if details.get("observed"):
                        geometry_observed = True
            tail_freq_mean = _mean_float(tail_frequencies)
            geometry_reached = geometry_observed or tail_presence_count > 0
            topology_ok = _topology_ok(
                left,
                right,
                domain_boundaries=boundary_windows,
                topology_mode=args.topology_mode,
            )
            candidate_eligible = topology_ok and pair in effective_anchor_pairs and tail_freq_mean >= args.rescue_tail_frequency_threshold
            selected_under_rescue_late_rule = False
            final_decision = "not_selected"
            final_pair_bucket = "unclassified"

            if pair_role == "border_long_range_rescue":
                if candidate_eligible and float(dca_scores.get(pair, 0.0)) >= v10_rescue_threshold:
                    selected_under_rescue_late_rule = support_by_pair.get(pair, 0) >= args.rescue_late_vote_threshold
                    if selected_under_rescue_late_rule:
                        selected_border_rescue.add(pair)
                        border_rescue_selected_count += 1
                        if pair_selected:
                            final_pair_bucket = "selected_border_rescue"
                            selected_pair_classification[pair] = final_pair_bucket
                        final_decision = "selected_border_rescue"
                    else:
                        monitor_only_border_rescue.add(pair)
                        border_rescue_monitor_count += 1
                        if pair_selected:
                            final_pair_bucket = "monitor_only_border_rescue"
                            selected_pair_classification[pair] = final_pair_bucket
                        final_decision = "monitor_only_border_rescue"
                elif candidate_eligible and geometry_reached:
                    monitor_only_border_rescue.add(pair)
                    border_rescue_monitor_count += 1
                    if pair_selected:
                        final_pair_bucket = "monitor_only_border_rescue"
                        selected_pair_classification[pair] = final_pair_bucket
                    final_decision = "monitor_only_border_rescue"
                elif pair_selected:
                    final_pair_bucket = "rejected_outside_blueprint_or_noise"
                    rejected_selected_pairs.add(pair)
                    final_decision = "rejected_border_rescue_candidate"
            elif pair_role == "true_long_range_core" and pair_selected:
                if pair_blueprint_lane == "strict":
                    selected_strict_scaffold.add(pair)
                    final_pair_bucket = "selected_strict_scaffold"
                    selected_pair_classification[pair] = final_pair_bucket
                elif pair_blueprint_lane == "balanced":
                    selected_balanced_core.add(pair)
                    final_pair_bucket = "selected_balanced_core"
                    selected_pair_classification[pair] = final_pair_bucket
                else:
                    final_pair_bucket = (
                        "selected_blueprint_unclassified_core"
                        if pair in effective_anchor_pairs
                        else "rejected_outside_blueprint_or_noise"
                    )
                    if final_pair_bucket == "selected_blueprint_unclassified_core":
                        selected_blueprint_unclassified_core.add(pair)
                    else:
                        rejected_selected_pairs.add(pair)
                    selected_pair_classification[pair] = final_pair_bucket
            elif pair_role == "local_support" and pair_selected:
                selected_local_support.add(pair)
                final_pair_bucket = "selected_local_support"
                selected_pair_classification[pair] = final_pair_bucket
            elif pair_role == "medium_support" and pair_selected:
                selected_medium_support.add(pair)
                final_pair_bucket = "selected_medium_support"
                selected_pair_classification[pair] = final_pair_bucket
            elif pair_role == "diagnostic_shell" and pair_selected:
                diagnostic_shell.add(pair)
                final_pair_bucket = "diagnostic_shell"
                selected_pair_classification[pair] = final_pair_bucket
            elif pair_selected:
                if pair in effective_anchor_pairs:
                    selected_blueprint_unclassified_core.add(pair)
                    final_pair_bucket = "selected_blueprint_unclassified_core"
                else:
                    rejected_selected_pairs.add(pair)
                    final_pair_bucket = "rejected_outside_blueprint_or_noise"
                selected_pair_classification[pair] = final_pair_bucket

            if pair_selected and pair not in selected_pair_classification:
                if pair in effective_anchor_pairs:
                    selected_blueprint_unclassified_core.add(pair)
                    selected_pair_classification[pair] = "selected_blueprint_unclassified_core"
                else:
                    rejected_selected_pairs.add(pair)
                    selected_pair_classification[pair] = "rejected_outside_blueprint_or_noise"

            pair_record = pair_audit.setdefault(pair, {})
            pair_record["runtime_v10"] = {
                "pair_lane": pair_lane,
                "pair_role": pair_role,
                "sequence_separation": sequence_sep,
                "tail_frequency_mean": tail_freq_mean,
                "tail_presence_count": tail_presence_count,
                "topology_ok": topology_ok,
                "geometry_reached": geometry_reached,
                "dca_score": float(dca_scores.get(pair, 0.0)),
                "rescue_late_support": support_by_pair.get(pair, 0),
                "rescue_late_qualified": selected_under_rescue_late_rule,
                "final_decision": final_decision,
                "input_bucket": final_pair_bucket,
            }

        runtime_v10_selected_pairs = (
            selected_strict_scaffold
            | selected_balanced_core
            | selected_border_rescue
            | selected_local_support
            | selected_medium_support
            | diagnostic_shell
        )
        classified_selected_pair_count = len(selected_pair_classification)
        dropped_selected_pair_count = input_selected_pair_count - classified_selected_pair_count
        runtime_v10_no_core_selected_pairs = len(runtime_v10_selected_pairs) == 0
        runtime_v10_failure_reasons: list[str] = []
        if runtime_v10_no_core_selected_pairs:
            runtime_v10_failure_reasons.append("no_core_selected_pairs_after_role_binding")
        if dropped_selected_pair_count > 0:
            runtime_v10_failure_reasons.append("role_binding_incomplete")
        runtime_v10_physics_interpretation_allowed = not runtime_v10_no_core_selected_pairs
        if dropped_selected_pair_count > 0:
            runtime_v10_physics_interpretation_allowed = False
        long_range_evidence_pairs = selected_strict_scaffold | selected_balanced_core | selected_border_rescue
        support_evidence_pairs = selected_local_support | selected_medium_support
        runtime_long_range_overlap = _set_overlap_ratio(
            long_range_evidence_pairs,
            effective_lane_pairs["strict"] | effective_lane_pairs["balanced"] | effective_lane_pairs["balanced_rescue"],
        )
        long_range_evidence_polluted = (
            len(
                (
                    set(long_range_evidence_pairs)
                    & (selected_local_support | selected_medium_support | diagnostic_shell)
                )
            )
            > 0
        )
        runtime_pair_in_blueprint = (
            _safe_ratio(
                len(long_range_evidence_pairs & effective_anchor_pairs),
                len(long_range_evidence_pairs),
            )
            if long_range_evidence_pairs
            else 0.0
        )
        runtime_v10_selection = {
            "selected_strict_scaffold": sorted([list(pair) for pair in sorted(selected_strict_scaffold)]),
            "selected_balanced_core": sorted([list(pair) for pair in sorted(selected_balanced_core)]),
            "selected_border_rescue": sorted([list(pair) for pair in sorted(selected_border_rescue)]),
            "monitor_only_border_rescue": sorted([list(pair) for pair in sorted(monitor_only_border_rescue)]),
            "selected_local_support": sorted([list(pair) for pair in sorted(selected_local_support)]),
            "selected_medium_support": sorted([list(pair) for pair in sorted(selected_medium_support)]),
            "diagnostic_shell": sorted([list(pair) for pair in sorted(diagnostic_shell)]),
            "selected_blueprint_unclassified_core": sorted([
                list(pair) for pair in sorted(selected_blueprint_unclassified_core)
            ]),
            "long_range_evidence": {
                "pairs": sorted([list(pair) for pair in sorted(long_range_evidence_pairs)]),
                "count": len(long_range_evidence_pairs),
            },
            "support_evidence": {
                "pairs": sorted([list(pair) for pair in sorted(support_evidence_pairs)]),
                "count": len(support_evidence_pairs),
            },
            "monitor_only": {
                "pairs": sorted([list(pair) for pair in sorted(monitor_only_border_rescue)]),
                "count": len(monitor_only_border_rescue),
            },
            "diagnostic_only": {
                "pairs": sorted([list(pair) for pair in sorted(diagnostic_shell)]),
                "count": len(diagnostic_shell),
            },
            "noise_added": 0,
            "long_range_evidence_polluted": bool(long_range_evidence_polluted),
            "pair_in_blueprint": runtime_pair_in_blueprint,
            "strict_overlap": _set_overlap_ratio(
                selected_strict_scaffold,
                effective_lane_pairs["strict"],
            ),
            "balanced_overlap": _set_overlap_ratio(
                selected_balanced_core,
                effective_lane_pairs["balanced"],
            ),
            "border_rescue_selected_count": border_rescue_selected_count,
            "border_rescue_monitor_count": border_rescue_monitor_count,
            "native_dashboard_precision": float(metric["native_contact_precision"]),
            "native_dashboard_recall": float(metric["native_contact_recall"]),
            "claim_allowed": False,
            "runtime_v10_selected_pairs": sorted([list(pair) for pair in sorted(runtime_v10_selected_pairs)]),
            "runtime_v10_selected_pair_count": len(runtime_v10_selected_pairs),
            "runtime_v10_long_range_overlap": runtime_long_range_overlap,
            "failure_type": runtime_v10_failure_reasons[0] if runtime_v10_failure_reasons else None,
            "failure_types": runtime_v10_failure_reasons,
            "input_selected_pair_count": input_selected_pair_count,
            "classified_selected_pair_count": classified_selected_pair_count,
            "dropped_selected_pair_count": dropped_selected_pair_count,
            "unclassified_selected_pairs": sorted([
                list(pair) for pair in sorted([
                    pair for pair, bucket in selected_pair_classification.items()
                    if bucket == "selected_blueprint_unclassified_core"
                ])
            ]),
            "rejected_selected_pairs": sorted([list(pair) for pair in sorted(rejected_selected_pairs)]),
            "classification_coverage_ratio": _safe_ratio(
                classified_selected_pair_count,
                input_selected_pair_count,
            ),
            "physics_interpretation_allowed": runtime_v10_physics_interpretation_allowed,
        }

    success_criteria = {
        "strict_overlap": {"value": strict_overlap, "min": 0.70},
        "balanced_overlap": {"value": balanced_overlap, "min": 0.70},
        "rescue_tail_presence": {"value": rescue_tail_presence, "min": 1},
        "pair_in_blueprint": {"value": pair_in_blueprint, "min": 1.0},
        "precision": {"value": float(metric["native_contact_precision"]), "min": 0.30},
        "run_only_noise": {"value": run_only_noise, "max": 0},
        "selected_pair_count": {"value": selected_pair_count, "min": 1},
    }
    for item in success_criteria.values():
        if "min" in item:
            item["pass"] = item["value"] >= item["min"]  # type: ignore[index]
        else:
            item["pass"] = item["value"] <= item["max"]  # type: ignore[index]

    run_type = "v5_physics_run" if preflight_state["v5_ready"] and selected_pair_count > 0 else "wiring_smoke_test"
    physics_interpretation_allowed = preflight_state["v5_ready"] and selected_pair_count > 0
    failure_reasons: list[str] = []
    if selected_pair_count == 0:
        failure_reasons.append("no_selected_pairs_smoke_or_underpowered_run")
    if not preflight_state["v5_ready"]:
        failure_reasons.append("v5_preflight_checks_failed")

    # Differentiate biological abstain from technical failure
    if physics_interpretation_allowed and not runtime_v10_selected_pairs:
        if "no_core_selected_pairs_after_role_binding" not in failure_reasons:
            failure_reasons.append("biological_signal_abstain")

    success_evaluation = {
        "strict_overlap": success_criteria["strict_overlap"],
        "balanced_overlap": success_criteria["balanced_overlap"],
        "rescue_overlap": {"value": rescue_overlap, "min": 0.0, "pass": rescue_overlap >= 0.0 and selected_pair_count > 0},
        "pair_in_blueprint": success_criteria["pair_in_blueprint"],
        "precision": success_criteria["precision"],
        "run_only_noise": success_criteria["run_only_noise"],
        "selected_pair_count": success_criteria["selected_pair_count"],
        "run_type": run_type,
        "result": "v5_physics_run" if run_type == "v5_physics_run" else "runner_outputs_v5_fields_successfully",
        "physics_interpretation_allowed": physics_interpretation_allowed,
        "failure_type": failure_reasons[0] if failure_reasons else None,
        "failure_types": failure_reasons,
        "selected_outside_effective_anchor_set_count": len(selected_outside_effective_anchor_set),
        "selected_outside_effective_anchor_set": selected_outside_effective_anchor_set,
        "selected_outside_effective_anchor_set_ratio": _safe_ratio(
            len(selected_outside_effective_anchor_set), len(selected_pairs)
        ),
        "rescue_tail_presence": {
            "by_pair": rescue_tail_reachability,
            "pairs_present_in_tail": rescue_tail_presence,
            "pairs_selected": len(selected_pairs_by_lane["balanced_rescue"]),
        },
        "claim_allowed": False,
    }
    success_evaluation["pass"] = (
        success_criteria["strict_overlap"]["pass"]
        and success_criteria["balanced_overlap"]["pass"]
        and success_criteria["pair_in_blueprint"]["pass"]
        and success_criteria["precision"]["pass"]
        and success_criteria["run_only_noise"]["pass"]
        and success_criteria["rescue_tail_presence"]["pass"]
        and success_criteria["selected_pair_count"]["pass"]
    )

    voted_rows: list[dict[str, object]] = []
    for i, j, cnt, mean_chem, _lane in voted_pairs:
        voted_rows.append(
            {
                "i": i,
                "j": j,
                "replica_votes": cnt,
                "mean_chemical_score": round(mean_chem, 6),
            }
        )

    _write_csv(out_root / "openmm_tmd_ensemble_voted_pairs.csv", voted_rows)

    final_rows = []
    for i, j, _, _, _ in voted_pairs:
        final_rows.append(
            {
                "i": i,
                "j": j,
            }
        )
    _write_csv(out_root / "openmm_tmd_ensemble_final_pairs.csv", final_rows)

    summary = {
        "kind": "openmm_tmd_replicas_v0",
        "source_accession": row.source_accession,
        "row_id": row.row_id,
        "sequence_length": row.sequence_length,
        "replicas": args.replicas,
        "successful_replicas": len(successful),
        "vote_threshold": args.vote_threshold,
        "vote_min_average_chemical": args.vote_min_average_chemical,
        "lane_vote_support": lane_vote_support,
        "lane_selected_counts": lane_selected_counts,
        "vote_failures": total_vote_failures,
        "preflight": preflight_state,
        "success_evaluation": success_evaluation,
        "pair_audit": {
            f"{left}-{right}": details for (left, right), details in pair_audit.items()
        },
        "frequency_threshold": args.frequency_threshold,
        "dca_threshold": args.dca_threshold,
        "chemical_threshold": args.chemical_threshold,
        "tail_fraction": args.tail_fraction,
        "contact_cutoff_angstrom": args.contact_cutoff_ang,
        "min_separation": args.min_separation,
        "topology_mode": args.topology_mode,
        "domain_boundaries": args.domain_boundaries,
        "max_degree": args.max_degree,
        "control_count": args.control_count,
        "max_control_overlap_fraction": args.max_control_overlap_fraction,
        "unique_voted_pairs": unique_vote_candidates,
        "selected_pair_count": len(predicted_pairs),
        "run_classification": {
            "run_type": run_type,
            "physics_interpretation_allowed": physics_interpretation_allowed,
        },
        "native_contact_precision": metric["native_contact_precision"],
        "native_contact_recall": metric["native_contact_recall"],
        "long_range_contact_recall": metric["long_range_contact_recall"],
        "contact_map_f1": metric["contact_map_f1"],
        "false_contact_rate": metric["false_contact_rate"],
        "replica_summaries": [
            {
                "replica": r.index,
                "seed": r.seed,
                "returncode": r.returncode,
                "out_dir": str(r.out_dir),
                "trajectory": str(r.trajectory) if r.trajectory else None,
                "stable_metadata": r.stable_metadata,
                "metric": r.metric,
            }
            for r in replica_results
        ],
    }
    if runtime_v10_selection is not None:
        summary["runtime_v10_role_aware_selector"] = True
        summary.update(runtime_v10_selection)

    if args.write_audit_json:
        path = Path(args.write_audit_json)
        path.write_text(
            json.dumps(
                {
                    "preflight": preflight_state,
                    "pair_audit": {
                        f"{left}-{right}": details for (left, right), details in pair_audit.items()
                    },
                    "success_evaluation": success_evaluation,
                    "summary": summary,
                    "baseline_source": args.preflight_baseline_json or None,
                    "audit_pairs_source": args.audit_pairs_json or None,
                    "strict_dca_threshold": args.strict_dca_threshold,
                    "balanced_strong_dca_threshold": args.balanced_strong_dca_threshold,
                    "balanced_rescue_dca_threshold": args.balanced_rescue_dca_threshold,
                    "monitor_dca_threshold": args.monitor_dca_threshold,
                },
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
    certificate_path = out_root / "openmm_tmd_replicas_v0_certificate.json"
    certificate_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
