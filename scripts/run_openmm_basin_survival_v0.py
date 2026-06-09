#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics as stats
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_five_axis_physics import matched_control_pairs
from pharmacotopology.folding_native_contact_eval import evaluate_contact_prediction
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)
from pharmacotopology.folding_template_docking import chemical_score


DEFAULT_BENCHMARK = Path("data/folding_real_coordinate_visual_8.locked.json")
DEFAULT_COUPLING_FILE = (
    Path("first_contact_clean_pharmacotopology_layer_run")
    / "query_centered_pfam00406_external_couplings_v0.locked.json"
)


FREQ_BUDGET = {"strict": 0.20, "balanced": 0.50, "broad": 1.00}
FREQ_PROFILE_WEIGHTS = {
    "strict": (0.72, 0.20, 0.08),
    "balanced": (0.38, 0.48, 0.14),
    "broad": (0.22, 0.60, 0.18),
}
SURVIVAL_SCORE_WEIGHTS = {
    "strict": (0.50, 0.25, 0.25),
    "balanced": (0.52, 0.30, 0.18),
    "broad": (0.65, 0.30, 0.05),
}
PLATEAU_PROFILE_BUDGET_DIVISORS = {"strict": 5, "balanced": 2, "broad": 1}
PLATEAU_PROFILE_NAMES = ("strict", "balanced", "broad", "diagnostic")
PLATEAU_CORE_PROFILE_NAMES = ("strict", "balanced", "broad")
PLATEAU_STABILITY_LOW_THRESHOLD = 0.25
PLATEAU_TAIL_SENSITIVITY_THRESHOLD = 0.20
PLATEAU_CUTOFF_SENSITIVITY_THRESHOLD = 0.20
PLATEAU_BASIN_STABILITY_THRESHOLD = 0.30
PLATEAU_STABLE_CORE_SUPPORT = 0.65
GATE_SEQUENCE = ("frequency", "chemical", "dca", "topology", "control", "survival", "degree", "budget", "collapse")
ABLATION_CUMULATIVE_MAP = {
    "frequency_only": {"frequency"},
    "budget_only": {"frequency", "budget"},
    "dca_only": {"frequency", "dca"},
    "topology_only": {"frequency", "chemical", "topology"},
    "control_gap_only": {"frequency", "chemical", "dca", "topology", "control"},
    "collapse_only": set(),
    "frequency_plus_budget": {"frequency", "budget"},
    "frequency_plus_budget_plus_collapse": {"frequency", "budget", "collapse"},
    "frequency_plus_budget_plus_dca": {"frequency", "budget", "dca"},
    "full": {"frequency", "chemical", "dca", "topology", "control", "survival", "degree", "budget"},
}
ABLATION_VALID_CASES = set(ABLATION_CUMULATIVE_MAP.keys())


@dataclass
class BasinPrediction:
    basin_id: int
    frame_count: int
    raw_pair_count: int
    predicted_count: int
    metric: dict
    contacts: list[tuple[int, int]]
    frequencies: dict[tuple[int, int], int]
    confidence: float
    collapse_score: float
    feature_centroid: dict[str, float]
    selected_frequency: float
    survival_mean_frequency: float
    survival_mean_score: float
    internal_quality: float
    dca_supported_count: int
    dca_met_ratio: float
    gate_counts: dict[str, int]
    ablation_counts: dict[str, int]
    ablation_profiles: dict[str, dict[str, object]]


def _parse_ca_trajectory(path: Path, chain_id: str = "A") -> list[dict[int, tuple[float, float, float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    current_model: dict[int, tuple[float, float, float]] = {}
    models: list[dict[int, tuple[float, float, float]]] = []

    in_multi = False
    for raw in lines:
        if raw.startswith("MODEL"):
            in_multi = True
            current_model = {}
            continue
        if raw.startswith("ENDMDL"):
            models.append(current_model)
            current_model = {}
            continue
        if not raw.startswith("ATOM"):
            continue
        if raw[12:16].strip() != "CA":
            continue
        chain = raw[21].strip() or "A"
        if chain != chain_id:
            continue
        try:
            idx = int(raw[22:26])
            x = float(raw[30:38])
            y = float(raw[38:46])
            z = float(raw[46:54])
        except Exception:
            continue
        current_model[idx] = (x, y, z)

    if current_model and not in_multi:
        models = [current_model]
    elif current_model:
        models.append(current_model)
    return models


def _parse_domain_boundaries(raw: str | None) -> tuple[tuple[int, int], ...]:
    if not raw:
        return ()
    out: list[tuple[int, int]] = []
    for segment in [piece.strip() for piece in raw.split(",") if piece.strip()]:
        if "-" not in segment:
            raise ValueError(f"Invalid domain boundary {segment!r}")
        start, end = segment.split("-", 1)
        s = int(start)
        e = int(end)
        if s <= 0 or e < s:
            raise ValueError(f"Invalid domain boundary {segment!r}")
        out.append((s, e))
    return tuple(sorted(out))


def _residue_domain_index(position: int, boundaries: tuple[tuple[int, int], ...]) -> int:
    for idx, (start, end) in enumerate(boundaries):
        if start <= position <= end:
            return idx
    return -1


def _topology_ok(left: int, right: int, *, boundaries: tuple[tuple[int, int], ...], topology_mode: str) -> bool:
    if topology_mode == "none" or not boundaries:
        return True
    left_domain = _residue_domain_index(left, boundaries)
    right_domain = _residue_domain_index(right, boundaries)
    if left_domain == -1 or right_domain == -1:
        return True
    if topology_mode == "interdomain":
        return left_domain != right_domain
    if topology_mode == "intradomain":
        return left_domain == right_domain
    return True


def _load_dca_scores(*, coupling_file: Path | None, row: RealCoordinateVisualRow, sequence_length: int) -> dict[tuple[int, int], float]:
    if not coupling_file or not coupling_file.exists():
        return {}
    try:
        dataset = load_coupling_dataset(coupling_file)
    except Exception:
        return {}
    scores = {}
    for c in dataset.constraints:
        if c.row_id != row.row_id and c.source_accession != row.source_accession:
            continue
        if not (1 <= c.i <= sequence_length and 1 <= c.j <= sequence_length):
            continue
        if c.i >= c.j:
            continue
        pair = (c.i, c.j)
        value = float(c.confidence)
        prev = scores.get(pair, -1.0)
        if value > prev:
            scores[pair] = value
    return scores


def _extract_contacts(coords: dict[int, tuple[float, float, float]], cutoff: float, min_separation: int = 2) -> dict[tuple[int, int], float]:
    if not coords:
        return {}
    keys = sorted(coords)
    cutoff_sq = cutoff * cutoff
    out: dict[tuple[int, int], float] = {}
    for left_idx, left in enumerate(keys[:-1]):
        lx, ly, lz = coords[left]
        for right in keys[left_idx + 1 :]:
            if right - left <= min_separation:
                continue
            rx, ry, rz = coords[right]
            dx = lx - rx
            dy = ly - ry
            dz = lz - rz
            d2 = dx * dx + dy * dy + dz * dz
            if d2 <= cutoff_sq:
                out[(left, right)] = math.sqrt(d2)
    return out


def _distance(coords: dict[int, tuple[float, float, float]], left: int, right: int) -> float:
    if left not in coords or right not in coords:
        return float("inf")
    lx, ly, lz = coords[left]
    rx, ry, rz = coords[right]
    dx = lx - rx
    dy = ly - ry
    dz = lz - rz
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _radius_of_gyration(coords: dict[int, tuple[float, float, float]]) -> float:
    if not coords:
        return 0.0
    xs = [p[0] for p in coords.values()]
    ys = [p[1] for p in coords.values()]
    zs = [p[2] for p in coords.values()]
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)
    total = 0.0
    for x, y, z in coords.values():
        total += (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2
    return math.sqrt(total / len(coords))


def _domain_centroid_distance(coords: dict[int, tuple[float, float, float]], boundaries: tuple[tuple[int, int], ...], seq_len: int) -> float:
    if not coords or not boundaries:
        split = seq_len // 2
        d1 = [r for r in range(1, split + 1) if r in coords]
        d2 = [r for r in range(split + 1, seq_len + 1) if r in coords]
    else:
        d1 = [r for r in range(boundaries[0][0], boundaries[0][1] + 1) if r in coords]
        d2 = [r for r in range(boundaries[1][0], boundaries[1][1] + 1) if r in coords] if len(boundaries) > 1 else []

    if not d1 or not d2:
        return 0.0
    c1 = [coords[r] for r in d1]
    c2 = [coords[r] for r in d2]
    cx1 = sum(p[0] for p in c1) / len(c1)
    cy1 = sum(p[1] for p in c1) / len(c1)
    cz1 = sum(p[2] for p in c1) / len(c1)
    cx2 = sum(p[0] for p in c2) / len(c2)
    cy2 = sum(p[1] for p in c2) / len(c2)
    cz2 = sum(p[2] for p in c2) / len(c2)
    return math.sqrt((cx1 - cx2) ** 2 + (cy1 - cy2) ** 2 + (cz1 - cz2) ** 2)


def _control_gap(row: RealCoordinateVisualRow, selected_pairs: Sequence[tuple[int, int]], candidate_pairs: Sequence[tuple[int, int]], control_count: int) -> dict[tuple[int, int], int]:
    if control_count <= 0 or not selected_pairs:
        return {}
    hits = defaultdict(int)
    for idx in range(1, control_count + 1):
        control_pairs = matched_control_pairs(row=row, selected_pairs=selected_pairs, candidate_pairs=candidate_pairs, control_index=idx)
        for pair in control_pairs:
            hits[pair] += 1
    return hits


def _apply_degree_cap(
    pairs: list[tuple[tuple[int, int], dict[str, float]]],
    sequence_length: int,
    max_degree: int,
) -> list[tuple[tuple[int, int], dict[str, float]]]:
    if max_degree <= 0 or not pairs:
        return pairs
    degree = [0] * (sequence_length + 1)
    out = []
    for pair, payload in pairs:
        left, right = pair
        if degree[left] >= max_degree or degree[right] >= max_degree:
            continue
        degree[left] += 1
        degree[right] += 1
        out.append((pair, payload))
    return out


def _resolve_auto_threshold(values: Sequence[float], default: float) -> float:
    sorted_values = sorted(set(float(v) for v in values), reverse=True)
    if len(sorted_values) < 2:
        return default
    best_gap = 0.0
    best_index = 0
    for idx, left in enumerate(sorted_values[:-1]):
        right = sorted_values[idx + 1]
        gap = left - right
        if gap > best_gap:
            best_gap = gap
            best_index = idx
    if best_gap <= 0:
        return default
    return max(0.0, sorted_values[best_index + 1])


def _percentile_95(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    if len(values) == 1:
        return float(values[0])
    sorted_values = sorted(values)
    idx = max(0, min(len(sorted_values) - 1, int(0.95 * (len(sorted_values) - 1))))
    return float(sorted_values[idx])


def _dca_mode_for_profile(profile: str) -> str:
    return "strict" if profile == "strict" else "soft"


def _resolve_cumulative_ablation_gates(case_name: str) -> set[str]:
    return set(ABLATION_CUMULATIVE_MAP.get(case_name, {case_name}) )


def _build_stable_pairs_for_frequency_threshold(
    pair_counts: Counter[tuple[int, int]],
    required_frames: int,
    row: RealCoordinateVisualRow,
    sequence: str,
    dca_scores: dict[tuple[int, int], float],
    dca_threshold: float,
    chemical_threshold: float,
    topology_mode: str,
    domain_boundaries: tuple[tuple[int, int], ...],
    control_count: int,
    max_control_overlap_fraction: float,
    max_degree: int,
    frequency_threshold: float,
    profile: str,
    dca_mode: str,
    *,
    enabled_gates: set[str] | None = None,
) -> tuple[list[tuple[tuple[int, int], dict[str, float]]], dict[str, int]]:
    if not pair_counts:
        return [], {
            "raw_pair_count": 0,
            "pre_control_pair_count": 0,
            "control_removed_pair_count": 0,
            "control_count": 0,
            "frequency_gate_pair_count": 0,
            "chemical_gate_pair_count": 0,
            "dca_gate_pair_count": 0,
            "topology_gate_pair_count": 0,
            "degree_culled_pair_count": 0,
            "dca_supported_pair_count": 0,
            "survival_scored_pair_count": 0,
            "budgeted_pair_count": 0,
        }

    gates = enabled_gates or {"frequency", "chemical", "dca", "topology", "control", "survival", "degree"}
    freq_enabled = "frequency" in gates
    chem_enabled = "chemical" in gates
    dca_enabled = "dca" in gates
    topo_enabled = "topology" in gates
    control_enabled = "control" in gates
    degree_enabled = "degree" in gates
    surv_enabled = "survival" in gates

    raw_pairs = 0
    freq_count = 0
    chem_count = 0
    dca_count = 0
    topo_count = 0

    pre_candidates = []
    for (left, right), count in pair_counts.items():
        raw_pairs += 1
        if left >= right or right > len(sequence):
            continue
        freq = count / required_frames
        if freq_enabled and freq < frequency_threshold:
            continue
        freq_count += 1
        chem = chemical_score(sequence[left - 1], sequence[right - 1])
        if chem_enabled and chem < chemical_threshold:
            continue
        chem_count += 1
        dca = dca_scores.get((left, right), 0.0)
        if dca_enabled and dca_mode == "strict" and dca < dca_threshold:
            continue
        dca_count += 1
        if topo_enabled and not _topology_ok(left, right, boundaries=domain_boundaries, topology_mode=topology_mode):
            continue
        topo_count += 1
        pre_candidates.append(
            ((left, right), {
                "count": count,
                "frequency": freq,
                "chemical_score": chem,
                "dca_score": dca,
            })
        )

    pre_control_count = len(pre_candidates)
    auto_control_count = control_count
    if auto_control_count <= 0:
        auto_control_count = 0
    elif auto_control_count > pre_control_count:
        auto_control_count = pre_control_count

    hits = _control_gap(row=row, selected_pairs=[pair for pair, _ in pre_candidates], candidate_pairs=tuple(sorted(pair_counts)), control_count=auto_control_count) if control_enabled else {}

    surviving = []
    control_removed = 0
    for pair, payload in pre_candidates:
        overlap = hits.get(pair, 0) / auto_control_count if auto_control_count else 0.0
        payload = dict(payload)
        payload["control_overlap_fraction"] = overlap
        payload["control_gap"] = 1.0 - overlap
        if control_enabled and overlap > max_control_overlap_fraction:
            control_removed += 1
            continue
        freq_weight, chem_weight, dca_weight = SURVIVAL_SCORE_WEIGHTS.get(profile, SURVIVAL_SCORE_WEIGHTS["balanced"])
        if not surv_enabled:
            payload["survival_score"] = payload["frequency"]
        else:
            payload["survival_score"] = (
                freq_weight * payload["frequency"]
                + chem_weight * payload["chemical_score"]
                + dca_weight * payload["dca_score"]
            )
        surviving.append((pair, payload))

    surviving.sort(key=lambda item: (-item[1]["survival_score"], -item[1]["frequency"], -item[1]["dca_score"], item[0][0], item[0][1]))
    before = len(surviving)
    surviving = _apply_degree_cap(surviving, sequence_length=len(sequence), max_degree=max_degree) if degree_enabled else surviving
    degree_culled = before - len(surviving)
    degree_scored_count = len(surviving)

    metadata = {
        "raw_pair_count": raw_pairs,
        "pre_control_pair_count": pre_control_count,
        "control_removed_pair_count": control_removed,
        "control_count": auto_control_count,
        "frequency_gate_pair_count": freq_count,
        "chemical_gate_pair_count": chem_count,
        "dca_gate_pair_count": dca_count,
        "topology_gate_pair_count": topo_count,
        "degree_culled_pair_count": degree_culled,
        "dca_supported_pair_count": len([pair for pair in pair_counts if dca_scores.get(pair, 0.0) >= dca_threshold]),
        "survival_scored_pair_count": len(surviving),
        "degree_scored_pair_count": degree_scored_count,
        "cumulative_after_frequency_pair_count": freq_count,
        "cumulative_after_chemical_pair_count": chem_count,
        "cumulative_after_dca_pair_count": dca_count,
        "cumulative_after_topology_pair_count": topo_count,
        "cumulative_after_control_pair_count": max(0, pre_control_count - control_removed),
        "cumulative_after_survival_pair_count": len(surviving),
        "cumulative_after_degree_pair_count": degree_scored_count,
        "budgeted_pair_count": 0,
    }
    return surviving, metadata


def _auto_basin_threshold(
    pair_counts: Counter[tuple[int, int]],
    required_frames: int,
    row: RealCoordinateVisualRow,
    sequence: str,
    dca_scores: dict[tuple[int, int], float],
    dca_threshold: float,
    chemical_threshold: float,
    topology_mode: str,
    domain_boundaries: tuple[tuple[int, int], ...],
    control_count: int,
    max_control_overlap_fraction: float,
    max_degree: int,
    profile: str,
    *,
    enabled_gates: set[str] | None = None,
) -> tuple[float, list[tuple[tuple[int, int], dict[str, float]]], dict[str, object]]:
    if not pair_counts:
        return 0.0, [], {
            "auto_profile": profile,
            "selected_frequency_threshold": 0.0,
            "frequency_budget_target": 0,
            "frequency_auto_selection_reason": "no_frequency_distribution",
            "frequency_auto_candidates": 0,
            "frequency_auto_valid_candidates": 0,
            "best_coverage_to_budget": 0.0,
            "best_mean_frequency": 0.0,
            "best_mean_survival_score": 0.0,
            "best_largest_count_drop": 0.0,
            "raw_pair_count": 0,
            "survival_scored_pair_count": 0,
            "budgeted_pair_count": 0,
        }

    frequencies = sorted(set(count / required_frames for count in pair_counts.values()), reverse=True)
    if not frequencies:
        return 0.0, [], {
            "auto_profile": profile,
            "selected_frequency_threshold": 0.0,
            "frequency_budget_target": 0,
            "frequency_auto_selection_reason": "no_frequency_distribution",
            "frequency_auto_candidates": 0,
            "frequency_auto_valid_candidates": 0,
            "best_coverage_to_budget": 0.0,
            "best_mean_frequency": 0.0,
            "best_mean_survival_score": 0.0,
            "best_largest_count_drop": 0.0,
            "raw_pair_count": 0,
            "survival_scored_pair_count": 0,
            "budgeted_pair_count": 0,
        }

    if frequencies[-1] > 0.0:
        frequencies.append(0.0)

    freq_weight, coverage_weight, surv_weight = FREQ_PROFILE_WEIGHTS.get(profile, FREQ_PROFILE_WEIGHTS["balanced"])
    budget = max(1, int(len(sequence) * FREQ_BUDGET.get(profile, 0.50)))

    best_score = -1.0
    best_threshold = frequencies[-1]
    best_payload: dict[str, object] = {}
    best_selected: list[tuple[tuple[int, int], dict[str, float]]] = []
    last_count = None
    valid_candidates = 0
    best_drop = 0.0

    for th in frequencies:
        stable_pairs, metadata = _build_stable_pairs_for_frequency_threshold(
            pair_counts,
            required_frames=required_frames,
            row=row,
            sequence=sequence,
            dca_scores=dca_scores,
            dca_threshold=dca_threshold,
            chemical_threshold=chemical_threshold,
            topology_mode=topology_mode,
            domain_boundaries=domain_boundaries,
            control_count=control_count,
            max_control_overlap_fraction=max_control_overlap_fraction,
            max_degree=max_degree,
            profile=profile,
            dca_mode=_dca_mode_for_profile(profile),
            frequency_threshold=th,
            enabled_gates=enabled_gates,
        )

        selected = stable_pairs[:budget]
        raw_selected_count = len(stable_pairs)
        selected_count = len(selected)
        if selected_count <= 0:
            continue
        valid_candidates += 1
        coverage = selected_count / budget
        mean_freq = sum(item[1]["frequency"] for item in selected) / selected_count
        mean_survival = sum(item[1]["survival_score"] for item in selected) / selected_count

        current_drop = 0.0
        if last_count is not None and last_count > raw_selected_count:
            current_drop = (last_count - raw_selected_count) / max(1, last_count)
        last_count = raw_selected_count
        best_drop = max(best_drop, current_drop)

        utility = freq_weight * mean_freq + coverage_weight * coverage + surv_weight * mean_survival
        candidate = {
            "frequency_threshold": th,
            "raw_selected_count": raw_selected_count,
            "selected_count": selected_count,
            "coverage": coverage,
            "mean_frequency": mean_freq,
            "mean_survival_score": mean_survival,
            "largest_count_drop": current_drop,
            **metadata,
        }
        candidate["budgeted_pair_count"] = selected_count
        if utility > best_score:
            best_score = utility
            best_threshold = th
            best_payload = candidate
            best_selected = selected

    if not best_payload:
        fallback = frequencies[-1]
        stable_pairs, metadata = _build_stable_pairs_for_frequency_threshold(
            pair_counts,
            required_frames=required_frames,
            row=row,
            sequence=sequence,
            dca_scores=dca_scores,
            dca_threshold=dca_threshold,
            chemical_threshold=chemical_threshold,
            topology_mode=topology_mode,
            domain_boundaries=domain_boundaries,
            control_count=control_count,
            max_control_overlap_fraction=max_control_overlap_fraction,
            max_degree=max_degree,
            profile=profile,
            dca_mode=_dca_mode_for_profile(profile),
            frequency_threshold=fallback,
            enabled_gates=enabled_gates,
        )
        metadata["budgeted_pair_count"] = len(stable_pairs[:budget])
        return fallback, stable_pairs[:budget], {
            **metadata,
            "auto_profile": profile,
            "selected_frequency_threshold": fallback,
            "frequency_budget_target": budget,
            "frequency_auto_selection_reason": "fallback_due_to_zero_utility",
            "frequency_auto_candidates": len(frequencies),
            "frequency_auto_valid_candidates": valid_candidates,
            "best_coverage_to_budget": 0.0,
            "best_mean_frequency": 0.0,
            "best_mean_survival_score": 0.0,
            "best_largest_count_drop": best_drop,
        }

    payload = dict(best_payload)
    payload.update(
        {
            "auto_profile": profile,
            "selected_frequency_threshold": best_threshold,
            "frequency_budget_target": budget,
            "frequency_auto_candidates": len(frequencies),
            "frequency_auto_valid_candidates": valid_candidates,
            "best_coverage_to_budget": payload["coverage"],
            "best_mean_frequency": payload["mean_frequency"],
            "best_mean_survival_score": payload["mean_survival_score"],
            "best_largest_count_drop": best_drop,
            "frequency_auto_selection_reason": "best_auto_utility",
        }
    )
    return float(best_threshold), best_selected, payload


def _load_trajectory_metadata(replica_dir: Path) -> dict:
    out = {}
    for key in ["openmm_dca_closure_md_certificate_v0.json", "openmm_contact_metric.json"]:
        cert = replica_dir / key
        if cert.exists():
            try:
                out[key] = json.loads(cert.read_text(encoding="utf-8"))
            except Exception:
                pass
    return out


def _load_anchor_csv(path: Path | None) -> list[tuple[int, int, float]]:
    if not path or not path.exists():
        return []
    anchors = []
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                i = int(row["i"])
                j = int(row["j"])
                d = float(row.get("distance_angstrom", row.get("distance", "6.0")))
            except Exception:
                continue
            if i >= j:
                i, j = j, i
            anchors.append((i, j, d))
    return anchors


def _frame_features(
    frames: list[dict[int, tuple[float, float, float]]],
    *,
    contact_cutoff: float,
    domain_boundaries: tuple[tuple[int, int], ...],
    anchor_pairs: list[tuple[int, int, float]],
    min_separation: int,
    seq_len: int,
) -> tuple[list[list[float]], list[dict[str, float]], list[dict[tuple[int, int], float]]]:
    feature_vectors: list[list[float]] = []
    per_frame: list[dict[str, float]] = []
    frame_contacts: list[dict[tuple[int, int], float]] = []
    pair_usage: Counter[tuple[int, int]] = Counter()

    for frame in frames:
        contacts = _extract_contacts(frame, cutoff=contact_cutoff, min_separation=min_separation)
        pair_usage.update(contacts)
        frame_contacts.append(contacts)

        if seq_len > 1:
            max_possible = (seq_len * (seq_len - 1)) / 2 - (seq_len - 1) * min_separation
            max_possible = max(1.0, max_possible)
        else:
            max_possible = 1.0

        density = len(contacts) / max_possible

        sat_count = 0
        for left, right, target_distance in anchor_pairs:
            d = _distance(frame, left, right)
            if d <= max(1e-6, target_distance):
                sat_count += 1
        sat_ratio = sat_count / max(1, len(anchor_pairs))

        # maximum node degree in this frame's contact map
        deg = defaultdict(int)
        for left, right in contacts:
            deg[left] += 1
            deg[right] += 1
        max_deg = max(deg.values()) if deg else 0
        if deg:
            high_degree_threshold = max(1, max(1, seq_len // 25))
            high_degree_fraction = sum(1 for value in deg.values() if value >= high_degree_threshold) / len(deg)
        else:
            high_degree_fraction = 0.0

        radius = _radius_of_gyration(frame)
        interdom = _domain_centroid_distance(frame, boundaries=domain_boundaries, seq_len=seq_len)

        feats = [radius, interdom, density, sat_ratio, float(max_deg), high_degree_fraction, 0.0]
        feature_vectors.append(feats)
        per_frame.append(
            {
                "radius_gyration": radius,
                "domain_gap": interdom,
                "contact_density": density,
                "anchor_satisfaction": sat_ratio,
                "max_degree": float(max_deg),
                "high_degree_fraction": float(high_degree_fraction),
                "contact_count": float(len(contacts)),
            }
        )

    if pair_usage:
        support_threshold = max(1, math.ceil(len(frames) * 0.45))
        core_pairs = {pair for pair, count in pair_usage.items() if count >= support_threshold}
    else:
        core_pairs = set()

    for feats, features, contacts in zip(feature_vectors, per_frame, frame_contacts):
        if core_pairs:
            union = len(contacts) + len(core_pairs) - len(set(contacts).intersection(core_pairs))
            similarity = 0.0 if union <= 0 else len(set(contacts).intersection(core_pairs)) / union
            feats[-1] = similarity
            features["contact_map_similarity"] = similarity

    return feature_vectors, per_frame, frame_contacts


def _cluster_basin_labels(feature_vectors: list[list[float]], basin_count: int) -> list[int]:
    if not feature_vectors:
        return []
    if len(feature_vectors) <= basin_count:
        return list(range(len(feature_vectors)))
    scaled = StandardScaler().fit_transform(feature_vectors)
    model = KMeans(n_clusters=basin_count, random_state=0, n_init=10)
    labels = model.fit_predict(scaled)
    return [int(x) for x in labels]


def _collapse_from_feature_stats(all_features: list[dict[str, float]], basin_indices: list[int], labels: list[int], basin: int) -> tuple[float, dict[str, float]]:
    feat_by_name = {"contact_density": [], "max_degree": [], "anchor_satisfaction": [], "high_degree_fraction": []}
    for i, label in enumerate(labels):
        if label != basin:
            continue
        f = all_features[i]
        for k in feat_by_name:
            feat_by_name[k].append(f[k])

    if not feat_by_name["contact_density"]:
        return 1.0, {"contact_density": 0.0, "max_degree": 0.0, "anchor_satisfaction": 0.0}

    d_mean = stats.mean(feat_by_name["contact_density"])
    md_mean = stats.mean(feat_by_name["max_degree"])
    d95 = _percentile_95([f["contact_density"] for f in all_features]) if all_features else d_mean
    md95 = _percentile_95([f["max_degree"] for f in all_features]) if all_features else md_mean
    anchor = stats.mean(feat_by_name["anchor_satisfaction"])

    # scale penalties to [0,1]
    density_penalty = 0.0 if d95 == 0 else min(1.0, d_mean / (d95 * 1.2))
    degree_penalty = 0.0 if md95 == 0 else min(1.0, md_mean / (md95 * 1.1))

    score = 0.7 * density_penalty + 0.3 * degree_penalty + 0.0 * (1.0 - anchor)
    return score, {
        "contact_density": d_mean,
        "max_degree": md_mean,
        "anchor_satisfaction": anchor,
        "density_quantile_ref": d95,
        "degree_quantile_ref": md95,
        "high_degree_fraction": stats.mean(feat_by_name["high_degree_fraction"]),
    }


def _select_contacts_for_frames(
    frames: list[dict[int, tuple[float, float, float]]],
    *,
    row: RealCoordinateVisualRow,
    dca_scores: dict[tuple[int, int], float],
    dca_threshold: float,
    contact_cutoff: float,
    chemical_threshold: float,
    topology_mode: str,
    domain_boundaries: tuple[tuple[int, int], ...],
    min_separation: int,
    profile: str,
    control_count: int,
    max_control_overlap_fraction: float,
    max_degree: int,
    enabled_gates: set[str] | None = None,
    ablation_gates: set[str] | None = None,
) -> BasinPrediction | None:
    if not frames:
        return None

    pair_counts: Counter[tuple[int, int]] = Counter()
    for frame in frames:
        for pair in _extract_contacts(frame, cutoff=contact_cutoff, min_separation=min_separation):
            pair_counts[pair] += 1

    raw_pair_count = len(pair_counts)
    if not pair_counts:
        return BasinPrediction(
            0,
            len(frames),
            0,
            0,
            {"native_contact_precision": 0.0, "native_contact_recall": 0.0, "long_range_contact_recall": 0.0, "contact_map_f1": 0.0, "false_contact_rate": 0.0},
            [],
            {},
            0.0,
            1.0,
            {},
            0.0,
            0.0,
            0.0,
            0.0,
            0,
            0.0,
            {},
            {},
            {},
        )

    dca_thr = dca_threshold
    if dca_thr == 0.0 and dca_scores:
        dca_thr = _resolve_auto_threshold(dca_scores.values(), default=0.0)

    threshold, selected_pairs, metadata = _auto_basin_threshold(
        pair_counts,
        required_frames=len(frames),
        row=row,
        sequence=row.sequence,
        dca_scores=dca_scores,
        dca_threshold=dca_thr,
        chemical_threshold=chemical_threshold,
        topology_mode=topology_mode,
        domain_boundaries=domain_boundaries,
        control_count=control_count,
        max_control_overlap_fraction=max_control_overlap_fraction,
        max_degree=max_degree,
        profile=profile,
        enabled_gates=enabled_gates,
    )
    if metadata.get("budgeted_pair_count", 0) == 0 and selected_pairs:
        metadata["budgeted_pair_count"] = len(selected_pairs)

    pre_control_count = int(metadata.get("pre_control_pair_count", 0))
    dca_supported_count = int(metadata.get("dca_supported_pair_count", 0))
    budget = int(metadata.get("frequency_budget_target", 0))
    dca_ratio = (len(selected_pairs) / max(1, dca_supported_count)) if dca_supported_count > 0 else (1.0 if not dca_scores else 0.0)

    ablation_counts: dict[str, int] = {}
    ablation_profiles: dict[str, dict[str, object]] = {}
    if ablation_gates:
        full_gates = {"frequency", "chemical", "dca", "topology", "control", "survival", "degree", "budget"}
        base = {
            "frequency_only": {"frequency"},
            "budget_only": {"budget"},
            "dca_only": {"frequency", "dca"},
            "topology_only": {"frequency", "chemical", "topology"},
            "control_gap_only": {"frequency", "chemical", "dca", "topology", "control"},
            "collapse_only": set(),
            "frequency_plus_budget": {"frequency", "budget"},
            "frequency_plus_budget_plus_collapse": {"frequency", "budget"},
            "frequency_plus_budget_plus_dca": {"frequency", "budget", "dca"},
            "full": full_gates,
        }
        for name in sorted(ablation_gates):
            if name not in base:
                continue
            if name == "collapse_only":
                continue
            gates = base[name]
            if name == "full":
                gates = full_gates
            one_gate_summary = _pair_gate_eval_for_ablation(
                pair_counts=pair_counts,
                required_frames=len(frames),
                row=row,
                sequence=row.sequence,
                dca_scores=dca_scores,
                dca_threshold=dca_thr,
                chemical_threshold=chemical_threshold,
                topology_mode=topology_mode,
                domain_boundaries=domain_boundaries,
                control_count=control_count,
                max_control_overlap_fraction=max_control_overlap_fraction,
                max_degree=max_degree,
                frequency_threshold=threshold,
                budget=budget,
                profile=profile,
                enabled_gates=gates,
            )
            cumulative_gates = _resolve_cumulative_ablation_gates(name)
            cumulative_summary = _pair_gate_eval_for_ablation(
                pair_counts=pair_counts,
                required_frames=len(frames),
                row=row,
                sequence=row.sequence,
                dca_scores=dca_scores,
                dca_threshold=dca_thr,
                chemical_threshold=chemical_threshold,
                topology_mode=topology_mode,
                domain_boundaries=domain_boundaries,
                control_count=control_count,
                max_control_overlap_fraction=max_control_overlap_fraction,
                max_degree=max_degree,
                frequency_threshold=threshold,
                budget=budget,
                profile=profile,
                enabled_gates=cumulative_gates,
            )
            profile_payload: dict[str, object] = {
                "one_gate_only": {
                    "selected_pair_count": int(one_gate_summary.get("selected_pair_count", 0)),
                    "selected_by_budget": int(one_gate_summary.get("budgeted_pair_count", 0)),
                    "survival_scored_pair_count": int(one_gate_summary.get("survival_scored_pair_count", 0)),
                    "degree_scored_pair_count": int(one_gate_summary.get("degree_scored_pair_count", 0)),
                },
                "cumulative_gate": {
                    "selected_pair_count": int(cumulative_summary.get("selected_pair_count", 0)),
                    "selected_by_budget": int(cumulative_summary.get("budgeted_pair_count", 0)),
                    "survival_scored_pair_count": int(cumulative_summary.get("survival_scored_pair_count", 0)),
                    "degree_scored_pair_count": int(cumulative_summary.get("degree_scored_pair_count", 0)),
                },
            }
            if name == "dca_only":
                soft_summary = _pair_gate_eval_for_ablation(
                    pair_counts=pair_counts,
                    required_frames=len(frames),
                    row=row,
                    sequence=row.sequence,
                    dca_scores=dca_scores,
                    dca_threshold=dca_thr,
                    chemical_threshold=chemical_threshold,
                    topology_mode=topology_mode,
                    domain_boundaries=domain_boundaries,
                    control_count=control_count,
                    max_control_overlap_fraction=max_control_overlap_fraction,
                    max_degree=max_degree,
                    frequency_threshold=threshold,
                    budget=budget,
                    profile="balanced",
                    enabled_gates={"frequency", "dca"},
                    dca_mode="soft",
                )
                hard_summary = _pair_gate_eval_for_ablation(
                    pair_counts=pair_counts,
                    required_frames=len(frames),
                    row=row,
                    sequence=row.sequence,
                    dca_scores=dca_scores,
                    dca_threshold=dca_thr,
                    chemical_threshold=chemical_threshold,
                    topology_mode=topology_mode,
                    domain_boundaries=domain_boundaries,
                    control_count=control_count,
                    max_control_overlap_fraction=max_control_overlap_fraction,
                    max_degree=max_degree,
                    frequency_threshold=threshold,
                    budget=budget,
                    profile="strict",
                    enabled_gates={"frequency", "dca"},
                    dca_mode="strict",
                )
                profile_payload["soft_score_only"] = {
                    "selected_pair_count": int(soft_summary.get("selected_pair_count", 0)),
                    "selected_by_budget": int(soft_summary.get("budgeted_pair_count", 0)),
                }
                profile_payload["hard_gate_final"] = {
                    "selected_pair_count": int(hard_summary.get("selected_pair_count", 0)),
                    "selected_by_budget": int(hard_summary.get("budgeted_pair_count", 0)),
                }
                profile_payload["dca_as_soft_score_topK"] = {
                    "selected_pair_count": int(soft_summary.get("selected_pair_count", 0)),
                    "selected_by_budget": int(soft_summary.get("budgeted_pair_count", 0)),
                }
                profile_payload["dca_as_hard_gate_topK"] = {
                    "selected_pair_count": int(hard_summary.get("selected_pair_count", 0)),
                    "selected_by_budget": int(hard_summary.get("budgeted_pair_count", 0)),
                }

            ablation_profiles[name] = profile_payload
            ablation_counts[name] = int(cumulative_summary.get("budgeted_pair_count", 0))

    selected = [pair for pair, _ in selected_pairs]
    metric = evaluate_contact_prediction(
        native_pairs=row.native_contact_pairs(),
        predicted_pairs=set(selected),
        long_range_threshold=24,
        short_range_threshold=12,
    ).to_dict()

    return BasinPrediction(
        basin_id=-1,
        frame_count=len(frames),
        raw_pair_count=raw_pair_count,
        predicted_count=len(selected),
        metric=metric,
        contacts=selected,
        frequencies=dict(pair_counts),
        confidence=float(len(selected) / budget) if budget else float(len(selected) > 0),
        collapse_score=0.0,
        feature_centroid={},
        selected_frequency=threshold,
        survival_mean_frequency=float(metadata.get("best_mean_frequency", 0.0)),
        survival_mean_score=float(metadata.get("best_mean_survival_score", 0.0)),
        internal_quality=0.0,
        dca_supported_count=dca_supported_count,
        dca_met_ratio=dca_ratio,
        gate_counts={
            "raw_pair_count": raw_pair_count,
            "pre_frequency_pair_count": int(metadata.get("raw_pair_count", raw_pair_count)),
            "frequency_pass_pair_count": int(metadata.get("frequency_gate_pair_count", 0)),
            "chemical_pass_pair_count": int(metadata.get("chemical_gate_pair_count", 0)),
            "dca_pass_pair_count": int(metadata.get("dca_gate_pair_count", 0)),
            "topology_pass_pair_count": int(metadata.get("topology_gate_pair_count", 0)),
            "pre_control_pair_count": pre_control_count,
            "control_removed_pair_count": int(metadata.get("control_removed_pair_count", 0)),
            "survival_scored_pair_count": int(metadata.get("survival_scored_pair_count", 0)),
            "degree_culled_pair_count": int(metadata.get("degree_culled_pair_count", 0)),
            "budgeted_pair_count": int(metadata.get("budgeted_pair_count", 0)),
            "final_selected_pair_count": int(metadata.get("selected_count", len(selected_pairs))),
        },
        ablation_counts=ablation_counts,
        ablation_profiles=ablation_profiles,
    )


def _contact_set_overlap(left: set[tuple[int, int]], right: set[tuple[int, int]]) -> float:
    if not left and not right:
        return 1.0
    union = left | right
    if not union:
        return 0.0
    return len(left & right) / len(union)


def _to_pair_obj(raw: object) -> tuple[int, int] | None:
    if isinstance(raw, (tuple, list)) and len(raw) >= 2:
        try:
            left = int(raw[0])
            right = int(raw[1])
            if left < right:
                return (left, right)
            return (right, left)
        except Exception:
            return None
    return None


def _extract_ranked_pairs(payload: object) -> list[tuple[tuple[int, int], float]]:
    out: list[tuple[tuple[int, int], float]] = []
    if not isinstance(payload, list):
        return out
    for item in payload:
        if isinstance(item, dict):
            pair = _to_pair_obj(item.get("pair"))
            if pair is None:
                continue
            score = float(item.get("count", 0.0) or 0.0)
            out.append((pair, score))
            continue
        if not isinstance(item, (tuple, list)) or len(item) < 2:
            continue
        pair = _to_pair_obj(item[:2])
        if pair is None:
            continue
        score = float(item[2]) if len(item) > 2 else 1.0
        out.append((pair, score))
    if not out:
        return out
    out.sort(key=lambda item: (-item[1], item[0][0], item[0][1]))
    return out


def _profile_top_k(sequence_length: int, profile: str) -> int:
    if profile == "strict":
        return max(1, sequence_length // PLATEAU_PROFILE_BUDGET_DIVISORS["strict"])
    if profile == "balanced":
        return max(1, sequence_length // PLATEAU_PROFILE_BUDGET_DIVISORS["balanced"])
    if profile == "broad":
        return max(1, sequence_length // PLATEAU_PROFILE_BUDGET_DIVISORS["broad"])
    if profile == "diagnostic":
        return -1
    return sequence_length


def _normalize_rank_payload(
    ranked: list[tuple[tuple[int, int], float]],
    top_k: int,
) -> tuple[list[tuple[int, int]], dict[tuple[int, int], int]]:
    if top_k < 0:
        top_k = len(ranked)
    top = [pair for pair, _ in ranked[:top_k]]
    positions: dict[tuple[int, int], int] = {pair: idx + 1 for idx, pair in enumerate(top)}
    return top, positions


def _rank_similarity(
    left: dict[tuple[int, int], int],
    right: dict[tuple[int, int], int],
    top_k: int,
) -> float:
    if top_k < 0:
        top_k = max(
            max(left.values(), default=1),
            max(right.values(), default=1),
        )
    if not left and not right:
        return 1.0
    union = set(left.keys()) | set(right.keys())
    if not union:
        return 1.0
    scale = float(top_k + 1)
    score = 0.0
    for pair in union:
        lpos = float(left.get(pair, top_k + 1))
        rpos = float(right.get(pair, top_k + 1))
        score += max(0.0, (scale - abs(lpos - rpos)) / scale)
    return score / len(union)


def _mean_matrix_off_diagonal(values: list[list[float]]) -> float:
    if len(values) < 2:
        return 1.0
    total = 0.0
    pairs = 0
    for i in range(len(values)):
        for j in range(i + 1, len(values[i])):
            total += values[i][j]
            pairs += 1
    return total / pairs if pairs else 1.0


def _neighbor_overlap_stability(
    entries: list[dict],
    dimension: str,
    profile_sets: list[set[tuple[int, int]]],
) -> float:
    if not entries:
        return 0.0
    groups: dict[int, list[int]] = defaultdict(list)
    for idx, entry in enumerate(entries):
        key = int(entry.get("tail_index" if dimension == "cutoff" else "cutoff_index", 0)) if dimension in {"tail", "cutoff"} else 0
        groups[key].append(idx)
    if not groups:
        return 1.0
    overlaps: list[float] = []
    for _, idxs in groups.items():
        idxs = sorted(set(idxs))
        if len(idxs) < 2:
            continue
        for a, b in zip(idxs, idxs[1:]):
            overlaps.append(_contact_set_overlap(profile_sets[a], profile_sets[b]))
    if not overlaps:
        return 1.0
    return sum(overlaps) / len(overlaps)


def _extract_gate_support_metrics(entries: list[dict]) -> tuple[float, float]:
    topo_support_values: list[float] = []
    control_gap_support_values: list[float] = []
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        internal_summary = entry.get("internal_summary")
        if not isinstance(internal_summary, dict):
            continue
        mean_gate = internal_summary.get("mean_gate_cull_profile")
        if not isinstance(mean_gate, dict):
            continue
        freq_pass = float(mean_gate.get("frequency_pass_pair_count", 0.0))
        topo_pass = float(mean_gate.get("topology_pass_pair_count", 0.0))
        pre_control = float(mean_gate.get("pre_control_pair_count", 0.0))
        control_removed = float(mean_gate.get("control_removed_pair_count", 0.0))

        topology = topo_pass / freq_pass if freq_pass > 0 else 0.0
        control_gap = 1.0 if pre_control <= 0.0 else max(0.0, 1.0 - (control_removed / pre_control))
        topo_support_values.append(min(1.0, max(0.0, topology)))
        control_gap_support_values.append(min(1.0, max(0.0, control_gap)))

    topology_support = stats.mean(topo_support_values) if topo_support_values else 1.0
    control_gap_support = stats.mean(control_gap_support_values) if control_gap_support_values else 1.0
    return topology_support, control_gap_support


def _evaluate_plateau_ranked_stability(
    configs: list[dict],
    *,
    sequence_length: int,
    row: RealCoordinateVisualRow,
    dca_scores: dict[tuple[int, int], float],
    dca_threshold: float,
    profile_names: Sequence[str] | None = None,
) -> dict[str, object]:
    active_profiles = tuple(profile_names) if profile_names else PLATEAU_PROFILE_NAMES
    profile_set = list(active_profiles)
    if "diagnostic" not in profile_set:
        profile_set.append("diagnostic")
    profile_set = [name for name in profile_set if name in PLATEAU_PROFILE_NAMES]
    if not profile_set:
        profile_set = list(PLATEAU_PROFILE_NAMES)

    topology_support, control_gap_support = _extract_gate_support_metrics(configs)
    if not configs:
        return {
            "selected_profile": "diagnostic",
            "profiles": {},
            "topK_overlap_matrix": [],
            "rank_stability": 0.0,
            "score_stability": 0.0,
            "cutoff_sensitivity": 0.0,
            "tail_sensitivity": 0.0,
            "basin_membership_stability": 0.0,
            "control_gap_support": float(control_gap_support),
            "topology_support": float(topology_support),
            "stable_core_count": 0,
            "stable_core_DCA_fraction": 0.0,
            "stable_core_collapse_score": 0.0,
            "stable_native_dashboard": {
                "native_contact_precision": 0.0,
                "native_contact_recall": 0.0,
                "long_range_contact_recall": 0.0,
                "contact_map_f1": 0.0,
                "false_contact_rate": 0.0,
            },
            "failure_hint": "empty_sweep",
            "profile_core_sizes": {},
        }

    profile_payloads: dict[str, object] = {}

    for profile_name in profile_set:
        top_k = _profile_top_k(sequence_length, profile_name)
        profile_sets: list[set[tuple[int, int]]] = []
        rank_positions: list[dict[tuple[int, int], int]] = []
        basin_sets: list[set[int]] = []
        for entry in configs:
            ranked_pairs = _extract_ranked_pairs(entry.get("ensemble_contact_ranked_pairs"))
            pairs, positions = _normalize_rank_payload(ranked_pairs, top_k)
            if not pairs:
                contact_pairs = [tuple(pair) for pair in entry.get("ensemble_contacts", []) if isinstance(pair, (list, tuple)) and len(pair) >= 2]
                pairs = [(_to_pair_obj(pair) or (0, 0)) for pair in contact_pairs]
                pairs = [pair for pair in pairs if pair != (0, 0)]
                if len(pairs) > top_k:
                    pairs = pairs[:top_k]
                positions = {pair: idx + 1 for idx, pair in enumerate(pairs)}
            profile_sets.append(set(pairs))
            rank_positions.append(positions)
            compatible_ids = entry.get("compatible_basin_ids", [])
            basin_sets.append(set(int(v) for v in compatible_ids if isinstance(v, int)) if isinstance(compatible_ids, list) else set())

        overlap_matrix: list[list[float]] = [
            [_contact_set_overlap(profile_sets[i], profile_sets[j]) for j in range(len(profile_sets))]
            for i in range(len(profile_sets))
        ]

        rank_sims: list[float] = []
        for i in range(len(profile_sets)):
            for j in range(i + 1, len(profile_sets)):
                rank_sims.append(_rank_similarity(rank_positions[i], rank_positions[j], top_k))

        basin_overlaps: list[float] = []
        for i in range(len(profile_sets)):
            for j in range(i + 1, len(profile_sets)):
                basin_overlaps.append(_contact_set_overlap(basin_sets[i], basin_sets[j]))

        top_k_pairs = [pair for pair in profile_sets for pair in pair]
        stable_counter: Counter[tuple[int, int]] = Counter(top_k_pairs)
        min_support = max(1, math.ceil(PLATEAU_STABLE_CORE_SUPPORT * len(configs)))
        stable_core = {pair for pair, freq in stable_counter.items() if freq >= min_support}
        dca_supported = sum(
            1
            for pair in stable_core
            if (dca_threshold <= 0 and dca_scores.get(pair, 0.0) > 0.0)
            or (dca_threshold > 0.0 and dca_scores.get(pair, 0.0) >= dca_threshold)
        )
        dca_fraction = float(dca_supported / len(stable_core)) if stable_core else 0.0
        collapse_scores = [
            float(entry.get("internal_summary", {}).get("mean_collapse_score", 0.0))
            for entry in configs
            if isinstance(entry.get("internal_summary"), dict)
        ]
        collapse_score = stats.mean(collapse_scores) if collapse_scores else 0.0

        profile_payloads[profile_name] = {
            "top_k": int(top_k),
            "overlap_matrix": overlap_matrix,
            "rank_stability": _mean_matrix_off_diagonal(overlap_matrix),
            "score_stability": sum(rank_sims) / len(rank_sims) if rank_sims else 1.0,
            "cutoff_sensitivity": _neighbor_overlap_stability(configs, "cutoff", profile_sets),
            "tail_sensitivity": _neighbor_overlap_stability(configs, "tail", profile_sets),
            "basin_membership_stability": sum(basin_overlaps) / len(basin_overlaps) if basin_overlaps else 1.0,
            "stable_core_count": int(len(stable_core)),
            "stable_core_DCA_fraction": float(dca_fraction),
            "stable_core_collapse_score": float(collapse_score),
            "stable_core_pairs": sorted(stable_core),
            "collapse_association": float(collapse_score),
            "control_gap_support": float(control_gap_support),
            "topology_support": float(topology_support),
        }

    selected_profile = "diagnostic"
    ranked_scores: list[tuple[str, float]] = []
    for profile in profile_set:
        if profile == "diagnostic" and "diagnostic" not in active_profiles:
            continue
        payload = profile_payloads.get(profile)
        if not isinstance(payload, dict):
            continue
        score = 0.55 * float(payload.get("rank_stability", 0.0)) + 0.45 * float(payload.get("score_stability", 0.0))
        ranked_scores.append((profile, score))
    if ranked_scores:
        selected_profile = max(ranked_scores, key=lambda item: item[1])[0]

    selected_payload = profile_payloads.get(selected_profile, {})
    if not isinstance(selected_payload, dict):
        selected_payload = {}

    stable_pairs_raw = selected_payload.get("stable_core_pairs", [])
    stable_pairs: set[tuple[int, int]] = set()
    if isinstance(stable_pairs_raw, list):
        for item in stable_pairs_raw:
            pair = _to_pair_obj(item)
            if pair is not None:
                stable_pairs.add(pair)

    stable_metric = {
        "native_contact_precision": 0.0,
        "native_contact_recall": 0.0,
        "long_range_contact_recall": 0.0,
        "contact_map_f1": 0.0,
        "false_contact_rate": 0.0,
    }
    if stable_pairs:
        stable_metric = evaluate_contact_prediction(
            native_pairs=row.native_contact_pairs(),
            predicted_pairs=stable_pairs,
            long_range_threshold=24,
            short_range_threshold=12,
        ).to_dict()

    selected_stable_dca_fraction = float(selected_payload.get("stable_core_DCA_fraction", 0.0))
    diagnostic_payload = profile_payloads.get("diagnostic", {})
    diagnostic_dca_fraction = float(diagnostic_payload.get("stable_core_DCA_fraction", 0.0) if isinstance(diagnostic_payload, dict) else 0.0)
    for profile_name in profile_set:
        candidate = profile_payloads.get(profile_name)
        if not isinstance(candidate, dict):
            continue
        candidate_dca_fraction = float(candidate.get("stable_core_DCA_fraction", 0.0))
        candidate_dca_enrichment = 0.0
        if diagnostic_dca_fraction > 0.0:
            candidate_dca_enrichment = candidate_dca_fraction / diagnostic_dca_fraction
        candidate["DCA_enrichment_vs_full"] = candidate_dca_enrichment
        dashboard = {
            "native_contact_precision": 0.0,
            "native_contact_recall": 0.0,
            "long_range_contact_recall": 0.0,
            "contact_map_f1": 0.0,
            "false_contact_rate": 0.0,
        }
        profile_pairs = set(candidate.get("stable_core_pairs", []))
        if profile_pairs:
            profile_metric = evaluate_contact_prediction(
                native_pairs=row.native_contact_pairs(),
                predicted_pairs=set(profile_pairs),
                long_range_threshold=24,
                short_range_threshold=12,
            ).to_dict()
            dashboard = {
                "native_contact_precision": float(profile_metric.get("native_contact_precision", 0.0)),
                "native_contact_recall": float(profile_metric.get("native_contact_recall", 0.0)),
                "long_range_contact_recall": float(profile_metric.get("long_range_contact_recall", 0.0)),
                "contact_map_f1": float(profile_metric.get("contact_map_f1", 0.0)),
                "false_contact_rate": float(profile_metric.get("false_contact_rate", 0.0)),
            }
        candidate["stable_native_dashboard"] = dashboard
        candidate["native_dashboard_precision"] = dashboard.get("native_contact_precision", 0.0)
        candidate["native_dashboard_recall"] = dashboard.get("native_contact_recall", 0.0)
        candidate["native_dashboard_long_range_recall"] = dashboard.get("long_range_contact_recall", 0.0)
        candidate["rank_stability"] = float(candidate.get("rank_stability", 0.0))
        candidate["score_stability"] = float(candidate.get("score_stability", 0.0))
        candidate["cutoff_sensitivity"] = float(candidate.get("cutoff_sensitivity", 0.0))
        candidate["tail_sensitivity"] = float(candidate.get("tail_sensitivity", 0.0))
        candidate["collapse_association"] = float(candidate.get("collapse_association", candidate.get("stable_core_collapse_score", 0.0)))
        candidate["control_gap_support"] = float(candidate.get("control_gap_support", control_gap_support))
        candidate["topology_support"] = float(candidate.get("topology_support", topology_support))

    selected_stable_payload = profile_payloads.get(selected_profile, {})

    return {
        "selected_profile": selected_profile,
        "profiles": profile_payloads,
        "topK_overlap_matrix": selected_payload.get("overlap_matrix", []),
        "rank_stability": float(selected_stable_payload.get("rank_stability", 0.0)),
        "score_stability": float(selected_stable_payload.get("score_stability", 0.0)),
        "cutoff_sensitivity": float(selected_stable_payload.get("cutoff_sensitivity", 0.0)),
        "tail_sensitivity": float(selected_stable_payload.get("tail_sensitivity", 0.0)),
        "basin_membership_stability": float(selected_payload.get("basin_membership_stability", 0.0)),
        "control_gap_support": float(selected_stable_payload.get("control_gap_support", control_gap_support)),
        "topology_support": float(selected_stable_payload.get("topology_support", topology_support)),
        "selected_dca_enrichment_vs_full": selected_stable_dca_fraction / diagnostic_dca_fraction if diagnostic_dca_fraction > 0.0 else 0.0,
        "stable_core_count": int(selected_payload.get("stable_core_count", 0)),
        "stable_core_DCA_fraction": float(selected_payload.get("stable_core_DCA_fraction", 0.0)),
        "stable_core_collapse_score": float(selected_payload.get("stable_core_collapse_score", 0.0)),
        "collapse_association": float(selected_payload.get("collapse_association", selected_payload.get("stable_core_collapse_score", 0.0))),
        "selected_native_dashboard": stable_metric,
        "stable_native_dashboard": stable_metric,
        "profile_core_sizes": {name: int(profile_payloads.get(name, {}).get("stable_core_count", 0)) for name in active_profiles},
        "failure_hint": "no_plateau_found" if len(configs) > 1 else "insufficient_configs",
    }


def _classify_no_plateau_by_readout(plateau_diag: dict[str, object]) -> tuple[str, str, str]:
    selected_profile = plateau_diag.get("selected_profile")
    profiles = plateau_diag.get("profiles")
    selected = {}
    if isinstance(profiles, dict) and isinstance(selected_profile, str):
        candidate = profiles.get(selected_profile)
        if isinstance(candidate, dict):
            selected = candidate

    stable_core_count = int(selected.get("stable_core_count", 0))
    if stable_core_count > 0:
        return "no_plateau_but_stable_core_found", "stable_core_detected", "review_stable_core_only"

    extinction = plateau_diag.get("extinction")
    if isinstance(extinction, dict) and extinction.get("is_extinct"):
        reason = str(extinction.get("reason") or "")
        subtype = str(extinction.get("subtype") or "")
        if reason in {"collapse_extinction", "topology_extinction", "dca_extinction", "control_gap_extinction", "voting_extinction"}:
            if reason == "collapse_extinction":
                return "gate_extinction", "collapse_extinction", f"inspect_extinction:{subtype}" if subtype else "inspect_gate_flow"
            if reason == "topology_extinction":
                return "gate_extinction", "topology_extinction", f"inspect_extinction:{subtype}" if subtype else "inspect_gate_flow"
            if reason == "dca_extinction":
                return "gate_extinction", "dca_extinction", f"inspect_extinction:{subtype}" if subtype else "inspect_gate_flow"
            if reason == "control_gap_extinction":
                return "gate_extinction", "control_gap_extinction", f"inspect_extinction:{subtype}" if subtype else "inspect_gate_flow"
            if reason == "voting_extinction":
                return "gate_extinction", "voting_extinction", f"inspect_extinction:{subtype}" if subtype else "inspect_gate_flow"

    stable_flag = float(selected.get("rank_stability", 0.0))
    tail_sens = float(selected.get("tail_sensitivity", 0.0))
    cutoff_sens = float(selected.get("cutoff_sensitivity", 0.0))
    basin_stability = float(selected.get("basin_membership_stability", 0.0))
    if tail_sens < PLATEAU_TAIL_SENSITIVITY_THRESHOLD and tail_sens <= cutoff_sens:
        return "no_plateau_due_to_broad_tail_noise", "weak_tail_stability", "sweep_finer_tail_grid"
    if cutoff_sens < PLATEAU_CUTOFF_SENSITIVITY_THRESHOLD and cutoff_sens <= tail_sens:
        return "no_plateau_due_to_cutoff_sensitivity", "weak_cutoff_stability", "sweep_finer_cutoff_grid"
    if stable_flag < PLATEAU_STABILITY_LOW_THRESHOLD and basin_stability < PLATEAU_BASIN_STABILITY_THRESHOLD:
        return "no_plateau_due_to_basin_instability", "bassin_membership_jitter", "repair_compatible_basins"
    if stable_flag < PLATEAU_STABILITY_LOW_THRESHOLD:
        return "no_plateau_due_to_basin_instability", "score_rank_jitter", "tighten_ensemble_or_gates"
    return "no_plateau_found", "insufficient_plateau_support", "expand_sweep_or_relax_gates"


def _parse_float_list(raw: str | None, *, default: Sequence[float], validate: str = "any") -> list[float]:
    if not raw:
        return sorted(set(float(value) for value in default))
    values: list[float] = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        value = float(item)
        if validate == "tail":
            if not (0 < value <= 1.0):
                raise ValueError(f"tail fractions must be in (0, 1], got {value}")
        elif validate == "cutoff":
            if value <= 0:
                raise ValueError(f"contact cutoffs must be > 0, got {value}")
        values.append(value)
    if not values:
        return sorted(set(float(value) for value in default))
    return sorted(set(values))


def _parse_gate_set(raw: str | None) -> set[str]:
    if not raw:
        return set()
    out = set(item.strip().lower() for item in raw.split(",") if item.strip())
    normalized: set[str] = set()
    for gate in out:
        if gate in {"all", "full"}:
            normalized.update({"frequency", "chemical", "dca", "topology", "control", "survival", "degree", "budget"})
            continue
        if gate in {
            "frequency_only",
            "budget_only",
            "dca_only",
            "topology_only",
            "control_gap_only",
            "collapse_only",
            "frequency_plus_budget",
            "frequency_plus_budget_plus_collapse",
            "frequency_plus_budget_plus_dca",
            "full",
        }:
            normalized.add(gate)
            continue
        if gate in {"frequency", "chem", "chemical", "dca", "topology", "control", "survival", "degree", "budget", "collapse"}:
            normalized.add("collapse" if gate == "collapse" else gate)
            if gate == "chem":
                normalized.add("chemical")
    return normalized


def _pair_gate_eval_for_ablation(
    pair_counts: Counter[tuple[int, int]],
    required_frames: int,
    row: RealCoordinateVisualRow,
    sequence: str,
    dca_scores: dict[tuple[int, int], float],
    dca_threshold: float,
    chemical_threshold: float,
    topology_mode: str,
    domain_boundaries: tuple[tuple[int, int], ...],
    control_count: int,
    max_control_overlap_fraction: float,
    max_degree: int,
    frequency_threshold: float,
    budget: int,
    profile: str,
    enabled_gates: set[str],
    dca_mode: str | None = None,
) -> dict[str, int]:
    dca_mode = dca_mode or _dca_mode_for_profile(profile)
    pairs, metadata = _build_stable_pairs_for_frequency_threshold(
        pair_counts,
        required_frames=required_frames,
        row=row,
        sequence=sequence,
        dca_scores=dca_scores,
        dca_threshold=dca_threshold,
        chemical_threshold=chemical_threshold,
        topology_mode=topology_mode,
        domain_boundaries=domain_boundaries,
        control_count=control_count,
        max_control_overlap_fraction=max_control_overlap_fraction,
        max_degree=max_degree,
        profile=profile,
        dca_mode=dca_mode,
        frequency_threshold=frequency_threshold,
        enabled_gates=enabled_gates,
    )
    budgeted = len(pairs[:budget]) if "budget" in enabled_gates else len(pairs)
    return {
        "frequency_gate_pair_count": int(metadata.get("frequency_gate_pair_count", 0)),
        "chemical_gate_pair_count": int(metadata.get("chemical_gate_pair_count", 0)),
        "dca_gate_pair_count": int(metadata.get("dca_gate_pair_count", 0)),
        "topology_gate_pair_count": int(metadata.get("topology_gate_pair_count", 0)),
        "pre_control_pair_count": int(metadata.get("pre_control_pair_count", 0)),
        "control_removed_pair_count": int(metadata.get("control_removed_pair_count", 0)),
        "degree_culled_pair_count": int(metadata.get("degree_culled_pair_count", 0)),
        "survival_scored_pair_count": int(metadata.get("survival_scored_pair_count", 0)),
        "degree_scored_pair_count": int(metadata.get("degree_scored_pair_count", 0)),
        "cumulative_after_frequency_pair_count": int(metadata.get("cumulative_after_frequency_pair_count", 0)),
        "cumulative_after_chemical_pair_count": int(metadata.get("cumulative_after_chemical_pair_count", 0)),
        "cumulative_after_dca_pair_count": int(metadata.get("cumulative_after_dca_pair_count", 0)),
        "cumulative_after_topology_pair_count": int(metadata.get("cumulative_after_topology_pair_count", 0)),
        "cumulative_after_control_pair_count": int(metadata.get("cumulative_after_control_pair_count", 0)),
        "cumulative_after_survival_pair_count": int(metadata.get("cumulative_after_survival_pair_count", 0)),
        "cumulative_after_degree_pair_count": int(metadata.get("cumulative_after_degree_pair_count", 0)),
        "selected_pair_count": int(len(pairs)),
        "budgeted_pair_count": int(budgeted),
    }


def _mean_gate_profile(predictions: list[BasinPrediction]) -> dict[str, float]:
    if not predictions:
        return {}
    keys = sorted({key for pred in predictions for key in pred.gate_counts})
    return {
        key: float(stats.mean([pred.gate_counts.get(key, 0) for pred in predictions]))
        for key in keys
    }


def _aggregate_replica_gate_profiles(replica_runs: list[dict]) -> dict[str, object]:
    profile_samples: dict[str, object] = {}

    def _collect(samples: dict[str, object], values: dict[str, object]) -> None:
        for key, value in values.items():
            if isinstance(value, dict):
                target = samples.setdefault(key, {})
                if isinstance(target, dict):
                    _collect(target, value)
                else:
                    continue
            elif isinstance(value, (int, float)):
                bucket = samples.setdefault(key, [])
                if isinstance(bucket, list):
                    bucket.append(float(value))
            elif isinstance(value, list):
                bucket = samples.setdefault(key, [])
                if isinstance(bucket, list):
                    total = sum(float(v) for v in value if isinstance(v, (int, float)))
                    if total:
                        bucket.append(total)

    for run in replica_runs:
        if "basins" not in run:
            continue
        for basin in run["basins"]:
            profiles = basin.get("ablation_profiles")
            if not isinstance(profiles, dict):
                continue
            _collect(profile_samples, profiles)

    output: dict[str, object] = {}

    def _to_means(samples: object) -> object:
        if isinstance(samples, list):
            return float(stats.mean(samples)) if samples else 0.0
        if not isinstance(samples, dict):
            return samples
        return {key: _to_means(value) for key, value in samples.items()}

    for key, value in profile_samples.items():
        output[key] = _to_means(value)
    return output


def _classify_extinction(selected_run: dict, *, collapse_reject_threshold: float = 0.75) -> dict[str, object]:
    gate_summary = selected_run.get("contact_gate_summary", {})
    if not isinstance(gate_summary, dict):
        return {"is_extinct": False, "reason": None, "subtype": None}

    final_count = int(gate_summary.get("final_predicted_count", 0))
    final_no_collapse = int(gate_summary.get("final_predicted_count_no_collapse", 0))
    if final_count > 0:
        return {"is_extinct": False, "reason": None, "subtype": None}

    basins = selected_run.get("basins", [])
    compatible_count = int(selected_run.get("compatible_basin_count", 0))
    internal_summary = selected_run.get("internal_summary", {})
    mean_profile = {}
    if isinstance(internal_summary, dict):
        mean_profile = internal_summary.get("mean_gate_cull_profile", {})
    if not isinstance(mean_profile, dict):
        mean_profile = {}

    raw = float(mean_profile.get("raw_pair_count", 0))
    dca_pass = float(mean_profile.get("dca_pass_pair_count", 0))
    topo_pass = float(mean_profile.get("topology_pass_pair_count", 0))
    freq_pass = float(mean_profile.get("frequency_pass_pair_count", 0))
    control_gap = float(mean_profile.get("control_removed_pair_count", 0))
    control_total = float(mean_profile.get("pre_control_pair_count", 0))
    survival_scored = float(mean_profile.get("survival_scored_pair_count", 0))
    degree_culled = float(mean_profile.get("degree_culled_pair_count", 0))

    basin_contact_total = sum(int(item.get("predicted_count", 0)) for item in basins if isinstance(item, dict))
    candidate_mean = float(stats.mean([float(item.get("collapse_score", 1.0)) for item in basins if isinstance(item, dict)]) if basins else 1.0)

    if dca_pass <= 0 and freq_pass > 0 and topo_pass > 0:
        return {"is_extinct": True, "reason": "dca_extinction", "subtype": "dca_gate_zero"}
    if topo_pass <= 0 and freq_pass > 0 and dca_pass > 0:
        return {"is_extinct": True, "reason": "topology_extinction", "subtype": "topology_cumulative_inconsistency"}
    if control_total > 0 and control_gap >= 0.75 * control_total:
        return {"is_extinct": True, "reason": "control_gap_extinction", "subtype": "control_overlap_filter"}
    if degree_culled > 0 and freq_pass > 0 and dca_pass > 0 and topo_pass > 0 and survival_scored > 0:
        return {"is_extinct": True, "reason": "collapse_extinction", "subtype": "degree_cap_extinction"}
    if basin_contact_total > 0 and (compatible_count >= 0) and final_no_collapse == 0:
        return {"is_extinct": True, "reason": "voting_extinction", "subtype": "basin_vote_all_or_nothing"}
    if raw <= 0:
        return {"is_extinct": True, "reason": "raw_signal_extinction", "subtype": "empty_input"}
    if candidate_mean >= collapse_reject_threshold and basin_contact_total > 0 and compatible_count > 0:
        return {"is_extinct": True, "reason": "collapse_extinction", "subtype": "collapse_threshold"}
    return {"is_extinct": False, "reason": None, "subtype": None}


def _classify_sweep_state(
    selected_config: dict,
    extinction: dict[str, object] | None = None,
    plateau_diagnostics: dict[str, object] | None = None,
) -> dict[str, object]:
    selected_contacts = int(selected_config.get("ensemble_contact_count", selected_config.get("ensemble_predicted_count", 0)))
    internal_score = float(selected_config.get("internal_stability_score", 0.0))
    plateau_size = int(selected_config.get("plateau_size", 0))
    extinction = extinction or {"is_extinct": False, "reason": None, "subtype": None}
    if selected_contacts <= 0 and not extinction.get("is_extinct"):
        return {
            "sweep_status": "degenerate_null_plateau",
            "signal_status": "no_surviving_contacts",
            "claim_allowed": False,
            "next_action": "stage_readout_ablation",
            "extinction": extinction,
        }
    if extinction.get("is_extinct"):
        reason = str(extinction.get("reason") or "gate_extinction")
        subtype = str(extinction.get("subtype") or "")
        next_action = "inspect_extinction"
        if subtype:
            next_action = f"inspect_extinction:{subtype}"
        return {
            "sweep_status": "gate_extinction",
            "signal_status": reason,
            "claim_allowed": False,
            "next_action": next_action,
            "extinction": extinction,
        }
    if plateau_size < 1:
        status, signal_status, next_action = _classify_no_plateau_by_readout(plateau_diagnostics or {})
        return {
            "sweep_status": status,
            "signal_status": signal_status,
            "claim_allowed": False,
            "next_action": next_action,
            "extinction": extinction,
        }
    if internal_score < 0.25:
        return {
            "sweep_status": "low_internal_stability",
            "signal_status": "weak_internal_signal",
            "claim_allowed": False,
            "next_action": "pipeline_tuning_readout",
            "extinction": extinction,
        }
    if selected_contacts > 900:
        return {
            "sweep_status": "degenerate_contact_soup",
            "signal_status": "contact_set_too_large",
            "claim_allowed": False,
            "next_action": "gate_tuning_ablation",
            "extinction": extinction,
        }
    return {
        "sweep_status": "stable_plateau",
        "signal_status": "surviving_contacts",
        "claim_allowed": True,
        "next_action": "proceed_to_follow_up",
        "extinction": extinction,
    }


def run_replica(
    replica_dir: Path,
    args,
    row: RealCoordinateVisualRow,
    dca_scores: dict[tuple[int, int], float],
    *,
    tail_fraction: float,
    contact_cutoff_ang: float,
) -> dict:
    trajectory_path = replica_dir / args.trajectory_name
    if not trajectory_path.exists():
        # fall back common names
        alt = replica_dir / "openmm_dca_restrained_trajectory.pdb"
        if alt.exists():
            trajectory_path = alt
        else:
            raise FileNotFoundError(f"Missing trajectory: {trajectory_path}")

    frames = _parse_ca_trajectory(trajectory_path)
    if not frames:
        raise FileNotFoundError(f"No frames in {trajectory_path}")

    tail_count = max(1, math.ceil(len(frames) * tail_fraction))
    frames = frames[-tail_count:]

    anchor_path = replica_dir / "target_anchors.csv"
    if not anchor_path.exists():
        anchor_path = replica_dir / "dca_anchors.csv"
    anchors = _load_anchor_csv(anchor_path)

    domain_boundaries = _parse_domain_boundaries(args.domain_boundaries)
    feature_vectors, frame_features, frame_contacts = _frame_features(
        frames,
        contact_cutoff=contact_cutoff_ang,
        domain_boundaries=domain_boundaries,
        anchor_pairs=anchors,
        min_separation=args.min_separation,
        seq_len=row.sequence_length,
    )

    labels = _cluster_basin_labels(feature_vectors, args.basin_count)
    basin_frames: dict[int, list[dict[int, tuple[float, float, float]]]] = defaultdict(list)
    ablation_gates = _parse_gate_set(args.gate_ablation_cases) if args.run_gate_readout and args.gate_ablation_cases else set()
    for frame, basin in zip(frames, labels):
        basin_frames[basin].append(frame)

    basin_results: list[BasinPrediction] = []
    for basin_id in sorted(basin_frames):
        pred = _select_contacts_for_frames(
            basin_frames[basin_id],
            row=row,
            dca_scores=dca_scores,
            dca_threshold=args.dca_threshold,
            contact_cutoff=contact_cutoff_ang,
            chemical_threshold=args.chemical_threshold,
            topology_mode=args.topology_mode,
            domain_boundaries=domain_boundaries,
            min_separation=args.min_separation,
            profile=args.frequency_auto_profile,
            control_count=args.control_count,
            max_control_overlap_fraction=args.max_control_overlap_fraction,
            max_degree=args.max_degree,
            enabled_gates={"frequency", "chemical", "dca", "topology", "control", "survival", "degree", "budget"},
            ablation_gates=ablation_gates if args.run_gate_readout else None,
        )
        if pred is None:
            continue
        pred_basin_features = [frame_features[i] for i, value in enumerate(labels) if value == basin_id]
        if not pred_basin_features:
            continue
        feature_centroid = {
            k: stats.mean([f[k] for f in pred_basin_features]) if pred_basin_features else 0.0
            for k in ["radius_gyration", "domain_gap", "contact_density", "anchor_satisfaction", "max_degree"]
        }
        pred = BasinPrediction(
            basin_id=basin_id,
            frame_count=len(basin_frames[basin_id]),
            raw_pair_count=pred.raw_pair_count,
            predicted_count=pred.predicted_count,
            metric=pred.metric,
            contacts=pred.contacts,
            frequencies=pred.frequencies,
            confidence=pred.confidence,
            collapse_score=0.0,
            feature_centroid=feature_centroid,
            selected_frequency=pred.selected_frequency,
            survival_mean_frequency=pred.survival_mean_frequency,
            survival_mean_score=pred.survival_mean_score,
            internal_quality=pred.internal_quality,
            dca_supported_count=pred.dca_supported_count,
            dca_met_ratio=pred.dca_met_ratio,
            gate_counts=pred.gate_counts,
            ablation_counts=pred.ablation_counts,
            ablation_profiles=pred.ablation_profiles,
        )
        basin_results.append(pred)

    # evaluate collapse scores for rejection
    all_feature_dicts = frame_features
    if not all_feature_dicts:
        all_feature_dicts = []
    for pred in basin_results:
        score, stats_map = _collapse_from_feature_stats(all_feature_dicts, list(range(len(labels))), labels, pred.basin_id)
        pred.feature_centroid.update(stats_map)
        pred.collapse_score = score
        if "collapse_only" in ablation_gates or "collapse_only" in _parse_gate_set(args.gate_ablation_cases):
            collapse_only_count = len(pred.contacts) if score < args.collapse_reject_threshold else 0
            pred.ablation_profiles["collapse_only"] = {
                "hard_gate_final": {
                    "selected_pair_count": int(collapse_only_count),
                    "selected_by_budget": int(collapse_only_count),
                },
                "one_gate_only": {
                    "selected_pair_count": int(collapse_only_count),
                    "selected_by_budget": int(collapse_only_count),
                },
            }
            pred.ablation_counts["collapse_only"] = int(collapse_only_count)

    baseline = sorted(
        [r for r in basin_results],
        key=lambda x: (
            x.dca_met_ratio,
            1.0 - x.collapse_score,
            x.survival_mean_score,
            x.survival_mean_frequency,
        ),
        reverse=True,
    )

    # internal per-basin score: avoid native labels at this point.
    for pred in basin_results:
        pred.internal_quality = 0.45 * (1.0 - pred.collapse_score)
        pred.internal_quality += 0.35 * min(1.0, pred.survival_mean_score)
        pred.internal_quality += 0.12 * min(1.0, pred.survival_mean_frequency)
        pred.internal_quality += 0.08 * min(1.0, pred.dca_met_ratio)
    baseline = sorted(basin_results, key=lambda x: x.internal_quality, reverse=True)

    # reject collapse-like basins and choose compatible ones for voting.
    if args.collapse_hard_reject:
        compatible = [r for r in basin_results if r.collapse_score < args.collapse_reject_threshold and r.internal_quality > 0.0]
        if not compatible:
            compatible = basin_results
    else:
        compatible = [r for r in basin_results if r.internal_quality > 0.0]

    compatible = sorted(compatible, key=lambda x: x.internal_quality, reverse=True)
    for pred in compatible:
        pred.confidence = pred.internal_quality

    # per-replica final voted contacts from compatible basins
    all_vote: Counter[tuple[int, int]] = Counter()
    weighted_vote: dict[tuple[int, int], float] = defaultdict(float)
    all_vote_no_collapse: Counter[tuple[int, int]] = Counter()
    for pred in basin_results:
        for p in pred.contacts:
            all_vote_no_collapse[p] += 1
    for pred in compatible[: max(1, args.compatible_basin_cap)]:
        for p in pred.contacts:
            all_vote[p] += 1
            weighted_vote[p] += pred.internal_quality

    final_pairs = set(pair for pair, count in all_vote.items() if count >= max(1, min(args.basin_vote_min, len(compatible))))
    final_pairs_no_collapse = set(pair for pair, count in all_vote_no_collapse.items() if count >= max(1, min(args.basin_vote_min, len(basin_results))))
    weighted_final_pairs = set(
        pair
        for pair, score in weighted_vote.items()
        if score >= max(0.5, min(args.basin_vote_min, max(1, len(compatible))) / 2.0)
    )
    final_metric = evaluate_contact_prediction(
        native_pairs=row.native_contact_pairs(),
        predicted_pairs=final_pairs,
        long_range_threshold=24,
        short_range_threshold=12,
    ).to_dict()
    weighted_metric = evaluate_contact_prediction(
        native_pairs=row.native_contact_pairs(),
        predicted_pairs=weighted_final_pairs,
        long_range_threshold=24,
        short_range_threshold=12,
    ).to_dict()

    mean_gate_cull_profile = _mean_gate_profile(basin_results)

    # basin-level and final writes
    output = {
        "source_accession": row.source_accession,
        "row_id": row.row_id,
        "replica": replica_dir.name,
        "trajectory": str(trajectory_path),
        "frame_count": len(frames),
        "tail_fraction": tail_fraction,
        "tail_frames": tail_count,
        "contact_cutoff_angstrom": contact_cutoff_ang,
        "basin_count": len(basin_frames),
        "compatible_basin_count": len(compatible),
        "basin_vote_threshold": max(1, min(args.basin_vote_min, len(compatible))),
        "frequency_profile": args.frequency_auto_profile,
        "internal_summary": {
            "compatible_basin_count": len(compatible),
            "mean_internal_quality": stats.mean([r.internal_quality for r in compatible]) if compatible else 0.0,
            "mean_collapse_score": stats.mean([r.collapse_score for r in compatible]) if compatible else 0.0,
            "mean_selected_frequency": stats.mean([r.selected_frequency for r in compatible]) if compatible else 0.0,
            "mean_survival_mean_frequency": stats.mean([r.survival_mean_frequency for r in compatible]) if compatible else 0.0,
            "mean_survival_mean_score": stats.mean([r.survival_mean_score for r in compatible]) if compatible else 0.0,
            "dca_supported_ratio": stats.mean([r.dca_met_ratio for r in compatible]) if compatible else 0.0,
            "dca_supported_count_mean": stats.mean([r.dca_supported_count for r in compatible]) if compatible else 0.0,
            "mean_gate_cull_profile": mean_gate_cull_profile,
        },
        "basins": [
            {
                "basin_id": pred.basin_id,
                "frame_count": pred.frame_count,
                "raw_pair_count": pred.raw_pair_count,
                "predicted_count": pred.predicted_count,
                "selected_frequency": pred.selected_frequency,
                "survival_mean_frequency": pred.survival_mean_frequency,
                "survival_mean_score": pred.survival_mean_score,
                "internal_quality": pred.internal_quality,
                "dca_supported_count": pred.dca_supported_count,
                "dca_met_ratio": pred.dca_met_ratio,
                "collapse_score": pred.collapse_score,
                "feature_centroid": pred.feature_centroid,
                **{k: pred.metric.get(k, 0.0) for k in [
                    "native_contact_precision",
                    "native_contact_recall",
                    "long_range_contact_recall",
                    "contact_map_f1",
                    "false_contact_rate",
                ]},
                "top_contacts": [list(pair) for pair in pred.contacts[: min(25, len(pred.contacts))]],
                "gate_counts": pred.gate_counts,
                "ablation_counts": pred.ablation_counts,
                "ablation_profiles": pred.ablation_profiles,
            }
            for pred in sorted(
                basin_results,
                key=lambda x: (x.internal_quality, 1.0 - x.collapse_score),
                reverse=True,
            )
        ],
        "contact_gate_summary": {
            "raw_pairs_total": int(sum(pred.raw_pair_count for pred in basin_results)),
            "compatible_basins": int(len(compatible)),
            "incompatible_basins": int(len(basin_results) - len(compatible)),
            "internal_summary_mean": float(stats.mean([pred.internal_quality for pred in compatible])) if compatible else 0.0,
            "mean_gate_cull_profile": mean_gate_cull_profile,
            "final_predicted_count": int(len(final_pairs)),
            "final_predicted_count_no_collapse": int(len(final_pairs_no_collapse)),
            "weighted_predicted_count": int(len(weighted_final_pairs)),
        },
        "ablation_profiles": _aggregate_replica_gate_profiles([{
            "basins": [
                {
                    "ablation_profiles": pred.ablation_profiles
                }
                for pred in basin_results
            ]
        }]),
        "compatible_basin_ids": [pred.basin_id for pred in compatible[: max(1, args.compatible_basin_cap)]],
        "final_vote_counts": {
            "with_collapse_filter": int(len(final_pairs)),
            "without_collapse_filter": int(len(final_pairs_no_collapse)),
            "weighted_vote_count": int(len(weighted_final_pairs)),
            "weighted_vote_threshold": float(max(0.5, min(args.basin_vote_min, max(1, len(compatible))) / 2.0)),
        },
        "predicted_contacts": sorted([list(pair) for pair in final_pairs]),
        "weighted_vote_score": weighted_metric,
        "contact_map_f1": final_metric["contact_map_f1"],
        "native_contact_precision": final_metric["native_contact_precision"],
        "native_contact_recall": final_metric["native_contact_recall"],
        "long_range_contact_recall": final_metric["long_range_contact_recall"],
        "false_contact_rate": final_metric["false_contact_rate"],
    }
    output["extinction"] = _classify_extinction(
        output,
        collapse_reject_threshold=args.collapse_reject_threshold if args.collapse_hard_reject else 1.0,
    )
    return output


def _run_ensemble(
    replica_paths: list[Path],
    args,
    row: RealCoordinateVisualRow,
    dca_scores: dict[tuple[int, int], float],
    *,
    tail_fraction: float,
    contact_cutoff_ang: float,
) -> tuple[dict, list[dict], set[tuple[int, int]]]:
    all_runs = []
    for rep in replica_paths:
        if not rep.is_dir():
            continue
        try:
            payload = run_replica(
                rep,
                args,
                row,
                dca_scores,
                tail_fraction=tail_fraction,
                contact_cutoff_ang=contact_cutoff_ang,
            )
        except Exception as exc:
            payload = {"replica": rep.name, "error": str(exc), "tail_fraction": tail_fraction, "contact_cutoff_angstrom": contact_cutoff_ang}
        all_runs.append(payload)

    # ensemble across replicas from compatible basins
    all_votes: Counter[tuple[int, int]] = Counter()
    internal_quality = []
    for item in all_runs:
        if "predicted_contacts" not in item:
            continue
        for pair in item["predicted_contacts"]:
            all_votes[tuple(pair)] += 1
        if "internal_summary" in item:
            internal_quality.append(float(item["internal_summary"]["mean_internal_quality"]))
    final_size = len(replica_paths) or 1
    ensemble_pairs = [pair for pair, count in all_votes.items() if count >= max(1, args.basin_vote_min)]
    ranked_pairs = [[pair[0], pair[1], count] for pair, count in all_votes.items()]
    if ranked_pairs:
        ranked_pairs.sort(key=lambda item: (-item[2], item[0], item[1]))
    ensemble_metric = evaluate_contact_prediction(
        native_pairs=row.native_contact_pairs(),
        predicted_pairs=set(tuple(pair) for pair in ensemble_pairs),
        long_range_threshold=24,
        short_range_threshold=12,
    ).to_dict()

    report = {
        "source_accession": args.source_accession,
        "row_id": row.row_id,
        "tail_fraction": tail_fraction,
        "contact_cutoff_angstrom": contact_cutoff_ang,
        "replica_count": len(replica_paths),
        "run_summaries": all_runs,
        "ensemble_predicted_count": len(ensemble_pairs),
        "replica_compatibility_basin_vote_min": args.basin_vote_min,
        "frequency_profile": args.frequency_auto_profile,
        "native_contact_precision": ensemble_metric["native_contact_precision"],
        "native_contact_recall": ensemble_metric["native_contact_recall"],
        "long_range_contact_recall": ensemble_metric["long_range_contact_recall"],
        "contact_map_f1": ensemble_metric["contact_map_f1"],
        "false_contact_rate": ensemble_metric["false_contact_rate"],
        "internal_stability_score": stats.mean(internal_quality) if internal_quality else 0.0,
        "ensemble_contacts": sorted([list(pair) for pair in ensemble_pairs]),
        "ensemble_contact_count": len(ensemble_pairs),
        "ensemble_contact_ranked_pairs": ranked_pairs,
        "ensemble_contact_vote_total": final_size,
        "ensemble_vote_threshold": max(1, args.basin_vote_min),
        "ablation_profiles": _aggregate_replica_gate_profiles(all_runs),
    }
    if args.run_gate_readout:
        run_extinctions = [r.get("extinction") for r in all_runs if isinstance(r, dict) and "extinction" in r]
        if run_extinctions and report["ensemble_contact_count"] <= 0:
            report["extinction"] = run_extinctions[0]
        else:
            report["extinction"] = {"is_extinct": False, "reason": None, "subtype": None}
    else:
        report["extinction"] = {"is_extinct": False, "reason": None, "subtype": None}
    return report, all_runs, set(ensemble_pairs)


def _find_plateau_by_overlap(
    configs: list[dict],
    *,
    overlap_threshold: float,
    min_plateau: int,
) -> tuple[dict | None, list[dict]]:
    if not configs:
        return None, []

    # adjacency by nearest-neighbor grid walk with overlap support
    index = {(item["tail_index"], item["cutoff_index"]): item for item in configs}
    neighbors = lambda i, j: ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1))

    visited: set[tuple[int, int]] = set()
    plateaus: list[dict] = []
    for key, item in index.items():
        if key in visited:
            continue
        queue = deque([key])
        component: list[dict] = []
        visited.add(key)
        while queue:
            node = queue.popleft()
            node_item = index[node]
            component.append(node_item)
            for neighbor in neighbors(*node):
                if neighbor not in index or neighbor in visited:
                    continue
                candidate = index[neighbor]
                overlap = _contact_set_overlap(
                    set(map(tuple, node_item["ensemble_contacts"])),
                    set(map(tuple, candidate["ensemble_contacts"])),
                )
                if overlap >= overlap_threshold:
                    visited.add(neighbor)
                    queue.append(neighbor)

        if len(component) >= max(1, min_plateau):
            mean_internal = stats.mean([float(entry["internal_stability_score"]) for entry in component])
            component_has_nonzero_contacts = any(len(entry.get("ensemble_contact_set", [])) > 0 for entry in component)
            if not component_has_nonzero_contacts:
                continue
            tail_values = sorted({entry["tail_fraction"] for entry in component})
            cutoff_values = sorted({entry["contact_cutoff_angstrom"] for entry in component})
            rep_pairs = [set(map(tuple, entry["ensemble_contacts"])) for entry in component]
            overlaps: list[float] = []
            for i, left in enumerate(rep_pairs):
                for right in rep_pairs[i + 1 :]:
                    overlaps.append(_contact_set_overlap(left, right))
            mean_overlap = stats.mean(overlaps) if overlaps else 1.0
            score = mean_internal * (0.75 + 0.25 * mean_overlap)
            tail_span = tail_values[-1] - tail_values[0] if len(tail_values) > 1 else 0.0
            cutoff_span = cutoff_values[-1] - cutoff_values[0] if len(cutoff_values) > 1 else 0.0

            center_tail = sum(tail_values) / len(tail_values)
            center_cutoff = sum(cutoff_values) / len(cutoff_values)
            representative = min(
                component,
                key=lambda entry: abs(entry["tail_fraction"] - center_tail) + abs(entry["contact_cutoff_angstrom"] - center_cutoff),
            )

            plateaus.append(
                {
                    "size": len(component),
                    "mean_internal_stability": mean_internal,
                    "mean_overlap": mean_overlap,
                    "tail_span": tail_span,
                    "cutoff_span": cutoff_span,
                    "score": score,
                    "tail_values": tail_values,
                    "contact_cutoff_values": cutoff_values,
                    "representative_tail": representative["tail_fraction"],
                    "representative_contact_cutoff": representative["contact_cutoff_angstrom"],
                    "config_indexes": [[entry["tail_index"], entry["cutoff_index"]] for entry in component],
                }
            )

    if not plateaus:
        return None, []
    plateaus.sort(key=lambda item: item["score"], reverse=True)
    winner = plateaus[0]
    selected_key = (winner["representative_tail"], winner["representative_contact_cutoff"])
    selected = next(item for item in configs if (item["tail_fraction"], item["contact_cutoff_angstrom"]) == selected_key)
    return selected, plateaus


def _load_row(rows_file: Path, source_accession: str) -> RealCoordinateVisualRow:
    for row in load_real_coordinate_visual_rows(rows_file):
        if row.source_accession == source_accession:
            return row
    raise SystemExit(f"No row for source_accession={source_accession}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Basin-aware openmm trajectory postprocessing")
    parser.add_argument("--source-accession", default="4AKE:A")
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK))
    parser.add_argument("--external-coupling-file", default=str(DEFAULT_COUPLING_FILE))
    parser.add_argument("--replica-root", required=True)
    parser.add_argument("--trajectory-name", default="openmm_dca_restrained_trajectory.pdb")
    parser.add_argument("--out-dir", default="first_contact_clean_pharmacotopology_layer_run/openmm_basin_survival_v0")
    parser.add_argument("--tail-fraction", type=float, default=0.20)
    parser.add_argument("--contact-cutoff-ang", type=float, default=7.0)
    parser.add_argument("--min-separation", type=int, default=2)
    parser.add_argument("--basin-count", type=int, default=5)
    parser.add_argument("--basin-vote-min", type=int, default=1)
    parser.add_argument("--compatible-basin-cap", type=int, default=2)
    parser.add_argument("--collapse-reject-threshold", type=float, default=0.75)
    parser.add_argument("--frequency-auto-profile", default="balanced", choices=("strict", "balanced", "broad"))
    parser.add_argument("--chemical-threshold", type=float, default=0.50)
    parser.add_argument("--dca-threshold", type=float, default=0.0)
    parser.add_argument("--topology-mode", choices=("none", "interdomain", "intradomain"), default="none")
    parser.add_argument("--domain-boundaries", default="")
    parser.add_argument("--max-degree", type=int, default=10)
    parser.add_argument("--control-count", type=int, default=0)
    parser.add_argument("--max-control-overlap-fraction", type=float, default=0.20)
    parser.add_argument("--run-one-replica", default="")
    parser.add_argument("--sweep", action="store_true", help="Run tail/cutoff sweep as collapse-aware calibration")
    parser.add_argument("--sweep-tail-fractions", default="0.10,0.15,0.20,0.30")
    parser.add_argument("--sweep-contact-cutoffs", default="6.5,7.0,7.5")
    parser.add_argument("--sweep-overlap-threshold", type=float, default=0.70)
    parser.add_argument("--sweep-min-plateau-size", type=int, default=2)
    parser.add_argument("--sweep-top-k", type=int, default=3)
    parser.add_argument("--run-gate-readout", action="store_true", help="Keep per-phase contact counts and optional one-gate ablation diagnostics.")
    parser.add_argument("--gate-ablation-cases", default="frequency_only,budget_only,dca_only,topology_only,control_gap_only,collapse_only,frequency_plus_budget,frequency_plus_budget_plus_collapse,frequency_plus_budget_plus_dca,full", help="Comma-separated ablation cases to evaluate when run-gate-readout is enabled.")
    parser.add_argument("--collapse-hard-reject", action="store_true", help="Apply hard collapse threshold rejection (legacy default behavior).")
    parser.add_argument("--stable-core-only", action="store_true", help="Evaluate only strict/balanced/broad ranked core envelopes first.")
    args = parser.parse_args()

    row = _load_row(Path(args.benchmark_file), args.source_accession)
    coupling_file = Path(args.external_coupling_file) if args.external_coupling_file else None
    dca_scores = _load_dca_scores(coupling_file=coupling_file, row=row, sequence_length=row.sequence_length)

    replica_root = Path(args.replica_root)
    if args.run_one_replica:
        replica_paths = [replica_root / args.run_one_replica]
    else:
        replica_paths = sorted(replica_root.glob("replica_*"))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.sweep:
        tail_values = _parse_float_list(args.sweep_tail_fractions, default=(0.10, 0.20, 0.30), validate="tail")
        cutoff_values = _parse_float_list(args.sweep_contact_cutoffs, default=(6.5, 7.0, 7.5), validate="cutoff")
        sweep_runs: list[dict] = []
        for tail_i, tail_fraction in enumerate(tail_values):
            for cutoff_j, contact_cutoff in enumerate(cutoff_values):
                report, runs, contact_set = _run_ensemble(
                    replica_paths=replica_paths,
                    args=args,
                    row=row,
                    dca_scores=dca_scores,
                    tail_fraction=tail_fraction,
                    contact_cutoff_ang=contact_cutoff,
                )
                report["tail_index"] = tail_i
                report["cutoff_index"] = cutoff_j
                report["ensemble_contact_count"] = len(contact_set)
                report["ensemble_contact_set"] = sorted([list(pair) for pair in contact_set])
                report["config_size"] = len(sweep_runs) + 1
                report["plateau_size"] = 0
                report["run_summaries"] = runs
                sweep_runs.append(report)

        selected_by_plateau, plateau_summaries = _find_plateau_by_overlap(
            sweep_runs,
            overlap_threshold=args.sweep_overlap_threshold,
            min_plateau=max(1, args.sweep_min_plateau_size),
        )
        plateau_found = selected_by_plateau is not None
        if selected_by_plateau is None:
            selected_by_plateau = max(sweep_runs, key=lambda item: item["internal_stability_score"] if "internal_stability_score" in item else 0.0)
            plateau_summaries = []
        selected_tail = selected_by_plateau["tail_fraction"]
        selected_cutoff = selected_by_plateau["contact_cutoff_angstrom"]
        selected_report, selected_runs, _ = _run_ensemble(
            replica_paths=replica_paths,
            args=args,
            row=row,
            dca_scores=dca_scores,
            tail_fraction=selected_tail,
            contact_cutoff_ang=selected_cutoff,
        )
        selected_report["plateau_size"] = 0
        for entry in plateau_summaries:
            if entry["representative_tail"] == selected_tail and entry["representative_contact_cutoff"] == selected_cutoff:
                selected_report["plateau_size"] = int(entry["size"])
        if not selected_report["plateau_size"] and plateau_found:
            selected_report["plateau_size"] = int(len(plateau_summaries[0].get("config_indexes", [])))
        profile_probe = PLATEAU_CORE_PROFILE_NAMES if args.stable_core_only else PLATEAU_PROFILE_NAMES
        plateau_diagnostics = _evaluate_plateau_ranked_stability(
            sweep_runs,
            sequence_length=len(row.sequence),
            row=row,
            dca_scores=dca_scores,
            dca_threshold=args.dca_threshold,
            profile_names=profile_probe,
        )
        plateau_diagnostics["extinction"] = selected_report.get("extinction", {"is_extinct": False, "reason": None, "subtype": None})
        selected_report["plateau_ranked_readout"] = plateau_diagnostics
        stable_profiles = plateau_diagnostics.get("profiles", {})
        selected_report["stable_core_profiles"] = {
            name: {
                "stable_core_count": int(data.get("stable_core_count", 0)),
                "stable_core_DCA_fraction": float(data.get("stable_core_DCA_fraction", 0.0)),
                "DCA_enrichment_vs_full": float(data.get("DCA_enrichment_vs_full", 0.0)),
                "rank_stability": float(data.get("rank_stability", 0.0)),
                "score_stability": float(data.get("score_stability", 0.0)),
                "cutoff_sensitivity": float(data.get("cutoff_sensitivity", 0.0)),
                "tail_sensitivity": float(data.get("tail_sensitivity", 0.0)),
                "collapse_association": float(data.get("collapse_association", data.get("stable_core_collapse_score", 0.0))),
                "control_gap_support": float(data.get("control_gap_support", 0.0)),
                "topology_support": float(data.get("topology_support", 0.0)),
                "native_dashboard_precision": float(data.get("native_dashboard_precision", data.get("stable_native_dashboard", {}).get("native_contact_precision", 0.0))),
                "native_dashboard_recall": float(data.get("native_dashboard_recall", data.get("stable_native_dashboard", {}).get("native_contact_recall", 0.0))),
                "native_dashboard_long_range_recall": float(data.get("native_dashboard_long_range_recall", data.get("stable_native_dashboard", {}).get("long_range_contact_recall", 0.0))),
            }
            for name, data in stable_profiles.items()
            if isinstance(data, dict) and name in profile_probe
        }
        classification = _classify_sweep_state(
            selected_report,
            extinction=selected_report.get("extinction"),
            plateau_diagnostics=plateau_diagnostics,
        )
        if args.stable_core_only:
            classification["claim_allowed"] = False
        selected_report["calibration_note"] = (
            "This is a garage calibration run. Native metrics are used only as "
            "diagnostic dashboard, not as selection authority."
        )
        selected_report.update(classification)
        summary_payload = {
            "calibration_note": selected_report["calibration_note"],
            **selected_report,
            "selected_tail_fraction": selected_tail,
            "selected_contact_cutoff_angstrom": selected_cutoff,
            "plateau_summaries": plateau_summaries,
            "plateau_found": bool(plateau_found),
            **classification,
            "sweep_best_config_index": selected_by_plateau.get("config_size", 1),
            "sweep_overlap_threshold": args.sweep_overlap_threshold,
            "sweep_min_plateau_size": args.sweep_min_plateau_size,
        }
        summary_path = out_dir / "openmm_basin_survival_summary.json"
        summary_path.write_text(json.dumps(selected_runs, indent=2, sort_keys=True), encoding="utf-8")
        sweep_payload = {
            "sweep_tail_fractions": tail_values,
            "sweep_contact_cutoffs": cutoff_values,
            "sweep_overlap_threshold": args.sweep_overlap_threshold,
            "sweep_min_plateau_size": args.sweep_min_plateau_size,
            "sweep_top_k": args.sweep_top_k,
            "plateau_rank_diagnostics": plateau_diagnostics,
            "plateau_summaries": plateau_summaries,
            "configs": [
                {
                    k: v
                    for k, v in entry.items()
                    if k
                    in {
                        "tail_fraction",
                        "contact_cutoff_angstrom",
                        "tail_index",
                        "cutoff_index",
                        "internal_stability_score",
                        "ensemble_predicted_count",
                        "ensemble_contact_count",
                        "plateau_size",
                        "contact_map_f1",
                        "sweep_status",
                        "signal_status",
                        "claim_allowed",
                        "next_action",
                        "native_contact_recall",
                        "native_contact_precision",
                        "weighted_vote_count",
                        "weighted_vote_score",
                    }
                }
                for entry in sweep_runs
            ],
            "selected": {
                "tail_fraction": selected_tail,
                "contact_cutoff_angstrom": selected_cutoff,
                "internal_stability_score": selected_report.get("internal_stability_score", 0.0),
                "ensemble_predicted_count": int(selected_report.get("ensemble_predicted_count", 0)),
                "ensemble_contact_count": int(selected_report.get("ensemble_contact_count", 0)),
                "contact_map_f1": selected_report.get("contact_map_f1", 0.0),
                "native_contact_recall": selected_report.get("native_contact_recall", 0.0),
                "native_contact_precision": selected_report.get("native_contact_precision", 0.0),
                **classification,
                "weighted_vote_count": int(selected_report.get("final_vote_counts", {}).get("weighted_vote_count", 0)),
                "weighted_vote_precision": selected_report.get("weighted_vote_score", {}).get("native_contact_precision", 0.0),
                "weighted_vote_recall": selected_report.get("weighted_vote_score", {}).get("native_contact_recall", 0.0),
                "weighted_vote_f1": selected_report.get("weighted_vote_score", {}).get("contact_map_f1", 0.0),
                "topK_overlap_matrix": plateau_diagnostics.get("topK_overlap_matrix", []),
                "rank_stability": plateau_diagnostics.get("rank_stability", 0.0),
                "score_stability": plateau_diagnostics.get("score_stability", 0.0),
                "cutoff_sensitivity": plateau_diagnostics.get("cutoff_sensitivity", 0.0),
                "tail_sensitivity": plateau_diagnostics.get("tail_sensitivity", 0.0),
                "basin_membership_stability": plateau_diagnostics.get("basin_membership_stability", 0.0),
                "stable_core_count": plateau_diagnostics.get("stable_core_count", 0),
                "stable_core_DCA_fraction": plateau_diagnostics.get("stable_core_DCA_fraction", 0.0),
                "DCA_enrichment_vs_full": plateau_diagnostics.get("selected_dca_enrichment_vs_full", 0.0),
                "stable_core_collapse_score": plateau_diagnostics.get("stable_core_collapse_score", 0.0),
                "collapse_association": plateau_diagnostics.get("collapse_association", plateau_diagnostics.get("stable_core_collapse_score", 0.0)),
                "control_gap_support": plateau_diagnostics.get("control_gap_support", 0.0),
                "topology_support": plateau_diagnostics.get("topology_support", 0.0),
                "stable_native_dashboard": plateau_diagnostics.get("stable_native_dashboard", {}),
                "stable_core_profiles": selected_report.get("stable_core_profiles", {}),
                "stable_core_only_mode": bool(args.stable_core_only),
                "extinction": selected_report.get("extinction", {"is_extinct": False, "reason": None, "subtype": None}),
                "ablation_profiles": selected_report.get("ablation_profiles", {}),
                "ranked_plateau_readout": plateau_diagnostics,
            },
        }
        sweep_path = out_dir / "openmm_basin_survival_sweep_summary.json"
        sweep_path.write_text(json.dumps(sweep_payload, indent=2, sort_keys=True), encoding="utf-8")
        report_path = out_dir / "openmm_basin_survival_report.json"
        report_path.write_text(json.dumps(summary_payload, indent=2, sort_keys=True), encoding="utf-8")
        print(json.dumps(summary_payload, indent=2, sort_keys=True))
    else:
        selected_report, selected_runs, _ = _run_ensemble(
            replica_paths=replica_paths,
            args=args,
            row=row,
            dca_scores=dca_scores,
            tail_fraction=args.tail_fraction,
            contact_cutoff_ang=args.contact_cutoff_ang,
        )
        selected_report["calibration_note"] = (
            "This is a garage calibration run. Native metrics are used only as "
            "diagnostic dashboard, not as selection authority."
        )
        summary_path = out_dir / "openmm_basin_survival_summary.json"
        summary_path.write_text(json.dumps(selected_runs, indent=2, sort_keys=True), encoding="utf-8")
        report_path = out_dir / "openmm_basin_survival_report.json"
        report_path.write_text(json.dumps(selected_report, indent=2, sort_keys=True), encoding="utf-8")
        print(json.dumps(selected_report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
