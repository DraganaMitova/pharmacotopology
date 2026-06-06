from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from math import sqrt
from pathlib import Path
from statistics import mean
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.artifact_io import write_csv_rows  # noqa: E402
from pharmacotopology.folding_coupling_negative_controls import (  # noqa: E402
    EXTERNAL_COUPLING_CONTROL_NAMES,
    generate_external_coupling_negative_controls,
)
from pharmacotopology.folding_evolutionary_constraints import (  # noqa: E402
    CouplingDataset,
)
from pharmacotopology.folding_external_coupling_importer import (  # noqa: E402
    import_external_coupling_dataset,
)
from pharmacotopology.folding_external_coupling_trace_loop import (  # noqa: E402
    TraceLoopRun,
    _build_multiscale_physical_contexts,
    _run_multiscale_phase_aligned_critical_boundary_selector,
    _run_multiscale_phase_aligned_external_novelty_boundary_selector,
    _run_multiscale_phase_aligned_footprint_novelty_boundary_selector,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)


BATTERY_REPORT_KIND = "blind_external_holdout_battery_v0"
BATTERY_CERTIFICATE_KIND = "blind_external_holdout_battery_certificate_v0"
DEFAULT_TARGET_MANIFEST = Path("data/blind_holdout_manifest_v0.locked.json")
DEFAULT_COUPLING_DIR = Path("data/blind_external_couplings_v0")
DEFAULT_OUTPUT_DIR = Path(
    "first_contact_clean_pharmacotopology_layer_run/blind_external_holdout_v0"
)
DEFAULT_SELECTOR = "external_multiscale_phase_aligned_footprint_novelty_boundary"
EXTERNAL_NOVELTY_SELECTOR = (
    "external_multiscale_phase_aligned_external_novelty_boundary"
)
FOOTPRINT_NOVELTY_SELECTOR = DEFAULT_SELECTOR
FRONTIER_SELECTOR = "external_multiscale_phase_aligned_external_novelty_frontier"
PRECISION_BOUNDARY_SELECTOR = (
    "external_multiscale_phase_aligned_critical_boundary"
)
SATURATION_BOUNDARY_SELECTOR = DEFAULT_SELECTOR
PASS_TARGET_WIN_RATE_MIN = 0.70
PAIR_RANDOMIZING_CONTROL_NAMES = frozenset(
    {
        "external_shuffled_same_row_same_separation",
        "external_cross_row_swapped",
        "external_random_long_range_same_count",
    }
)


@dataclass(frozen=True)
class ExactContactMetrics:
    exact_contact_precision_top_L: float
    exact_contact_precision_L_over_2: float
    exact_long_range_contact_precision: float
    exact_long_range_contact_recall: float
    region_width_penalty: float
    region_pair_over_native_hit_ratio: float

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ScaffoldContactMetrics:
    scaffold_contact_count: int
    scaffold_exact_contact_precision: float
    scaffold_exact_long_range_contact_recall: float
    scaffold_contact_precision_delta_vs_top_L: float
    scaffold_contact_recall_delta_vs_region: float
    scaffold_contact_map_perfect: bool
    compact_scaffold_contact_count: int
    compact_scaffold_exact_contact_precision: float
    compact_scaffold_exact_long_range_contact_recall: float
    compact_scaffold_precision_delta_vs_scaffold: float
    compact_scaffold_recall_delta_vs_scaffold: float
    compact_scaffold_contact_map_perfect: bool
    density_compact_scaffold_contact_count: int
    density_compact_scaffold_exact_contact_precision: float
    density_compact_scaffold_exact_long_range_contact_recall: float
    density_compact_scaffold_precision_delta_vs_compact: float
    density_compact_scaffold_recall_delta_vs_compact: float
    density_compact_scaffold_contact_map_perfect: bool
    adjacent_density_patch_scaffold_contact_count: int
    adjacent_density_patch_scaffold_exact_contact_precision: float
    adjacent_density_patch_scaffold_exact_long_range_contact_recall: float
    adjacent_density_patch_scaffold_precision_delta_vs_density_compact: float
    adjacent_density_patch_scaffold_recall_delta_vs_density_compact: float
    adjacent_density_patch_scaffold_contact_map_perfect: bool
    phase_field_scaffold_contact_count: int
    phase_field_scaffold_exact_contact_precision: float
    phase_field_scaffold_exact_long_range_contact_recall: float
    phase_field_scaffold_precision_delta_vs_adjacent_density_patch: float
    phase_field_scaffold_recall_delta_vs_adjacent_density_patch: float
    phase_field_scaffold_contact_map_perfect: bool
    phase_field_scaffold_mode: str
    phase_field_scaffold_radius: int
    phase_ribbon_bridge_scaffold_contact_count: int
    phase_ribbon_bridge_scaffold_exact_contact_precision: float
    phase_ribbon_bridge_scaffold_exact_long_range_contact_recall: float
    phase_ribbon_bridge_scaffold_precision_delta_vs_phase_field: float
    phase_ribbon_bridge_scaffold_recall_delta_vs_phase_field: float
    phase_ribbon_bridge_scaffold_contact_map_perfect: bool

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class ScaffoldContactCandidate:
    pair: tuple[int, int]
    confidence: float
    sequence_separation: int
    event_ids: tuple[str, ...]
    center_compactness: float
    separation_compactness: float

    @property
    def compactness_score(self) -> float:
        return (
            self.confidence
            * self.center_compactness
            * self.separation_compactness
        )


def _rounded(value: float) -> float:
    return round(value, 6)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(block)
    return f"sha256:{digest.hexdigest()}"


def _current_commit_hash() -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return ""
    return result.stdout.strip()


def _load_json(path: Path) -> Mapping[str, object]:
    parsed = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(parsed, Mapping):
        raise ValueError(f"{path} must contain a JSON object")
    return parsed


def _resolve_path(raw: object, *, base_dir: Path) -> Path:
    path = Path(str(raw))
    if path.is_absolute():
        return path
    return (base_dir / path).resolve()


def _target_coupling_file(
    target: Mapping[str, object],
    *,
    coupling_dir: Path,
) -> Path:
    for field_name in (
        "external_coupling_file",
        "external_coupling_filename",
        "coupling_file",
    ):
        raw = target.get(field_name)
        if raw:
            path = Path(str(raw))
            if path.is_absolute():
                return path
            repo_relative = REPO_ROOT / path
            if repo_relative.exists():
                return repo_relative
            return coupling_dir / path
    target_id = str(target["target_id"])
    return coupling_dir / f"{target_id}.locked.json"


def _constraint_precision(
    constraints: Sequence[object],
    native_pairs: set[tuple[int, int]],
) -> float:
    if not constraints:
        return 0.0
    supported = sum(
        1
        for constraint in constraints
        if getattr(constraint, "pair")() in native_pairs
    )
    return _rounded(supported / len(constraints))


def _exact_contact_metrics(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
    run: TraceLoopRun,
) -> ExactContactMetrics:
    constraints_by_row = dataset.constraints_by_row_id()
    top_l_precisions: list[float] = []
    top_half_precisions: list[float] = []
    long_precisions: list[float] = []
    long_recalls: list[float] = []
    for row in rows:
        native_pairs = set(row.native_contact_pairs())
        native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
        constraints = tuple(constraints_by_row.get(row.row_id, ()))
        ranked = tuple(
            sorted(
                constraints,
                key=lambda constraint: (
                    constraint.rank or 10**9,
                    -constraint.confidence,
                    constraint.i,
                    constraint.j,
                ),
            )
        )
        top_l = ranked[: row.sequence_length]
        half_count = max(1, row.sequence_length // 2)
        top_half = ranked[:half_count]
        long_constraints = tuple(
            constraint
            for constraint in ranked
            if constraint.sequence_separation >= 24
        )
        supported_long = {
            constraint.pair()
            for constraint in long_constraints
            if constraint.pair() in native_long
        }
        top_l_precisions.append(_constraint_precision(top_l, native_pairs))
        top_half_precisions.append(_constraint_precision(top_half, native_pairs))
        long_precisions.append(_constraint_precision(long_constraints, native_pairs))
        long_recalls.append(
            _rounded(len(supported_long) / len(native_long))
            if native_long
            else 1.0
        )

    possible_region_pair_count = sum(
        len(event.candidate_region_pairs())
        for event in run.selected_events
    )
    native_hit_count = sum(
        event.native_contact_count_after_scoring
        for event in run.selected_events
    )
    cluster_precision = run.metric.contact_cluster_precision
    return ExactContactMetrics(
        exact_contact_precision_top_L=_rounded(mean(top_l_precisions))
        if top_l_precisions
        else 0.0,
        exact_contact_precision_L_over_2=_rounded(mean(top_half_precisions))
        if top_half_precisions
        else 0.0,
        exact_long_range_contact_precision=_rounded(mean(long_precisions))
        if long_precisions
        else 0.0,
        exact_long_range_contact_recall=_rounded(mean(long_recalls))
        if long_recalls
        else 0.0,
        region_width_penalty=_rounded(1.0 - cluster_precision),
        region_pair_over_native_hit_ratio=_rounded(
            possible_region_pair_count / native_hit_count
        )
        if native_hit_count
        else 0.0,
    )


def _center_compactness_for_event(
    *,
    pair: tuple[int, int],
    event: object,
) -> float:
    left_width = max(1, event.segment_a_end - event.segment_a_start + 1)
    right_width = max(1, event.segment_b_end - event.segment_b_start + 1)
    left_center = (event.segment_a_start + event.segment_a_end) / 2
    right_center = (event.segment_b_start + event.segment_b_end) / 2
    left_offset = abs(pair[0] - left_center) / max(1.0, left_width / 2)
    right_offset = abs(pair[1] - right_center) / max(1.0, right_width / 2)
    return max(0.0, 1.0 - (left_offset + right_offset) / 2)


def _compact_scaffold_core(
    candidates: Sequence[ScaffoldContactCandidate],
) -> tuple[ScaffoldContactCandidate, ...]:
    if len(candidates) <= 1:
        return tuple(candidates)
    ordered = tuple(
        sorted(
            candidates,
            key=lambda candidate: (
                candidate.compactness_score,
                candidate.confidence,
                len(candidate.event_ids),
                -candidate.sequence_separation,
                candidate.pair,
            ),
            reverse=True,
        )
    )
    required_events = {
        event_id
        for candidate in ordered
        for event_id in candidate.event_ids
    }
    covered_events: set[str] = set()
    coverage_index = 0
    for index, candidate in enumerate(ordered):
        covered_events.update(candidate.event_ids)
        if covered_events == required_events:
            coverage_index = index
            break
    scores = [candidate.compactness_score for candidate in ordered]
    gaps = [
        (scores[index] - scores[index + 1], index)
        for index in range(coverage_index, len(scores) - 1)
    ]
    if not gaps:
        return ordered[: coverage_index + 1]
    _, boundary_index = max(
        gaps,
        key=lambda item: (item[0], scores[item[1]], -item[1]),
    )
    return ordered[: boundary_index + 1]


def _local_contact_density(
    *,
    candidate: ScaffoldContactCandidate,
    candidates: Sequence[ScaffoldContactCandidate],
    radius: int,
) -> float:
    density = 0.0
    for neighbor in candidates:
        if neighbor.pair == candidate.pair:
            continue
        distance = max(
            abs(candidate.pair[0] - neighbor.pair[0]),
            abs(candidate.pair[1] - neighbor.pair[1]),
        )
        if 0 < distance <= radius:
            density += (
                neighbor.confidence
                * (radius + 1 - distance)
                / radius
            )
    return density


def _contact_map_distance(
    left: tuple[int, int],
    right: tuple[int, int],
) -> int:
    return max(abs(left[0] - right[0]), abs(left[1] - right[1]))


def _density_compact_score(
    candidate: ScaffoldContactCandidate,
    candidates: Sequence[ScaffoldContactCandidate],
) -> float:
    density = max(
        _local_contact_density(
            candidate=candidate,
            candidates=candidates,
            radius=radius,
        )
        for radius in (2, 4, 8)
    )
    return candidate.compactness_score * (1.0 + density)


def _self_critical_scored_prefix(
    candidates: Sequence[ScaffoldContactCandidate],
    *,
    score_candidates: Sequence[ScaffoldContactCandidate],
    minimum_count: int,
) -> tuple[ScaffoldContactCandidate, ...]:
    if not candidates:
        return ()
    ordered = tuple(
        sorted(
            candidates,
            key=lambda candidate: (
                _density_compact_score(candidate, score_candidates),
                candidate.compactness_score,
                candidate.confidence,
                len(candidate.event_ids),
                -candidate.sequence_separation,
                candidate.pair,
            ),
            reverse=True,
        )
    )
    minimum_index = max(0, min(len(ordered) - 1, minimum_count - 1))
    if len(ordered) <= minimum_index + 1:
        return ordered
    scores = [
        _density_compact_score(candidate, score_candidates)
        for candidate in ordered
    ]
    gaps = [
        (scores[index] - scores[index + 1], index)
        for index in range(minimum_index, len(scores) - 1)
    ]
    if not gaps:
        return ordered[: minimum_index + 1]
    _, boundary_index = max(
        gaps,
        key=lambda item: (item[0], scores[item[1]], -item[1]),
    )
    return ordered[: boundary_index + 1]


def _density_compact_scaffold_core(
    candidates: Sequence[ScaffoldContactCandidate],
) -> tuple[ScaffoldContactCandidate, ...]:
    if len(candidates) <= 1:
        return tuple(candidates)
    ordered = tuple(
        sorted(
            candidates,
            key=lambda candidate: (
                _density_compact_score(candidate, candidates),
                candidate.compactness_score,
                candidate.confidence,
                len(candidate.event_ids),
                -candidate.sequence_separation,
                candidate.pair,
            ),
            reverse=True,
        )
    )
    required_events = {
        event_id
        for candidate in ordered
        for event_id in candidate.event_ids
    }
    covered_events: set[str] = set()
    coverage_index = 0
    for index, candidate in enumerate(ordered):
        covered_events.update(candidate.event_ids)
        if covered_events == required_events:
            coverage_index = index
            break
    scores = [
        _density_compact_score(candidate, candidates)
        for candidate in ordered
    ]
    gaps = [
        (scores[index] - scores[index + 1], index)
        for index in range(coverage_index, len(scores) - 1)
    ]
    if not gaps:
        return ordered[: coverage_index + 1]
    _, boundary_index = max(
        gaps,
        key=lambda item: (item[0], scores[item[1]], -item[1]),
    )
    return ordered[: boundary_index + 1]


def _adjacent_density_patch_scaffold_core(
    candidates: Sequence[ScaffoldContactCandidate],
) -> tuple[ScaffoldContactCandidate, ...]:
    all_candidates = tuple(candidates)
    density_core = _density_compact_scaffold_core(all_candidates)
    if not density_core:
        return ()
    core_pairs = {candidate.pair for candidate in density_core}
    patch_candidates = tuple(
        candidate
        for candidate in all_candidates
        if candidate.pair in core_pairs
        or any(
            _contact_map_distance(candidate.pair, core_pair) <= 1
            for core_pair in core_pairs
        )
    )
    selected = _self_critical_scored_prefix(
        patch_candidates,
        score_candidates=all_candidates,
        minimum_count=len(density_core),
    )
    selected_by_pair = {candidate.pair: candidate for candidate in selected}
    selected_by_pair.update(
        {candidate.pair: candidate for candidate in density_core}
    )
    return tuple(selected_by_pair.values())


def _population_stddev(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    center = mean(values)
    return sqrt(mean((value - center) ** 2 for value in values))


def _pair_correlation(pairs: Sequence[tuple[int, int]]) -> float:
    if len(pairs) <= 1:
        return 0.0
    left = [float(pair[0]) for pair in pairs]
    right = [float(pair[1]) for pair in pairs]
    left_center = mean(left)
    right_center = mean(right)
    left_ss = sum((value - left_center) ** 2 for value in left)
    right_ss = sum((value - right_center) ** 2 for value in right)
    if left_ss == 0.0 or right_ss == 0.0:
        return 0.0
    covariance = sum(
        (left_value - left_center) * (right_value - right_center)
        for left_value, right_value in zip(left, right)
    )
    return covariance / sqrt(left_ss * right_ss)


def _nearest_phase_votes(
    pairs: Sequence[tuple[int, int]],
) -> tuple[int, int]:
    diagonal_votes = 0
    square_votes = 0
    for pair in pairs:
        nearest: tuple[int, int, int] | None = None
        for other in pairs:
            if other == pair:
                continue
            delta_left = other[0] - pair[0]
            delta_right = other[1] - pair[1]
            distance = max(abs(delta_left), abs(delta_right))
            candidate = (distance, delta_left, delta_right)
            if nearest is None or candidate < nearest:
                nearest = candidate
        if nearest is None:
            continue
        distance, delta_left, delta_right = nearest
        if abs(delta_left - delta_right) <= 1:
            diagonal_votes += 1
        if distance <= 2:
            square_votes += 1
    return diagonal_votes, square_votes


def _phase_field_mode_and_radius(
    seeds: Sequence[ScaffoldContactCandidate],
) -> tuple[str, int]:
    pairs = tuple(candidate.pair for candidate in seeds)
    if len(pairs) <= 1:
        return "point", 1
    offsets = [float(right - left) for left, right in pairs]
    sums = [float(right + left) for left, right in pairs]
    offset_spread = _population_stddev(offsets)
    sum_spread = _population_stddev(sums)
    correlation = max(0.0, _pair_correlation(pairs))
    diagonal_votes, square_votes = _nearest_phase_votes(pairs)
    diagonal_score = correlation * sum_spread / (offset_spread + 1.0)
    square_score = (square_votes / len(pairs)) * sqrt(offset_spread + 1.0)
    if diagonal_score > square_score:
        return "diagonal", max(1, int(round(offset_spread)))
    radius = max(1, int(round(sqrt(len(pairs)))))
    return "square", radius


def _phase_field_support(
    *,
    row_length: int,
    candidates: Sequence[ScaffoldContactCandidate],
    seeds: Sequence[ScaffoldContactCandidate],
    radius: int,
    mode: str,
) -> tuple[dict[tuple[int, int], float], set[tuple[int, int]]]:
    candidates_by_pair = {candidate.pair: candidate for candidate in candidates}
    direct_pairs = set(candidates_by_pair)
    max_seed_score = (
        max(_density_compact_score(candidate, candidates) for candidate in seeds)
        if seeds
        else 1.0
    ) or 1.0
    support: dict[tuple[int, int], float] = {}
    for seed in seeds:
        seed_left, seed_right = seed.pair
        seed_score = _density_compact_score(seed, candidates) / max_seed_score
        for delta_left in range(-radius, radius + 1):
            for delta_right in range(-radius, radius + 1):
                if mode == "diagonal" and abs(delta_left - delta_right) > 1:
                    continue
                pair_left = seed_left + delta_left
                pair_right = seed_right + delta_right
                if (
                    pair_left < 1
                    or pair_right > row_length
                    or pair_right - pair_left < 24
                ):
                    continue
                distance = max(abs(delta_left), abs(delta_right))
                if distance > radius:
                    continue
                kernel = (radius + 1 - distance) / (radius + 1)
                score = seed_score * kernel
                direct_candidate = candidates_by_pair.get((pair_left, pair_right))
                if direct_candidate is not None:
                    score += (
                        0.35
                        * _density_compact_score(direct_candidate, candidates)
                        / max_seed_score
                    )
                support[(pair_left, pair_right)] = (
                    support.get((pair_left, pair_right), 0.0) + score
                )
    return support, direct_pairs


def _phase_field_prefix(
    support: Mapping[tuple[int, int], float],
    *,
    minimum_count: int,
) -> set[tuple[int, int]]:
    if not support:
        return set()
    ordered = sorted(
        support.items(),
        key=lambda item: (item[1], -item[0][0], -item[0][1]),
        reverse=True,
    )
    minimum_index = max(0, min(len(ordered) - 1, minimum_count - 1))
    scores = [score for _, score in ordered]
    gaps = [
        (scores[index] - scores[index + 1], index)
        for index in range(minimum_index, len(scores) - 1)
    ]
    if gaps:
        _, boundary_index = max(
            gaps,
            key=lambda item: (item[0], scores[item[1]], -item[1]),
        )
    else:
        boundary_index = minimum_index
    return {pair for pair, _ in ordered[: boundary_index + 1]}


def _phase_field_scaffold_core(
    *,
    row_length: int,
    candidates: Sequence[ScaffoldContactCandidate],
) -> tuple[set[tuple[int, int]], str, int]:
    all_candidates = tuple(candidates)
    seeds = _adjacent_density_patch_scaffold_core(all_candidates)
    if not seeds:
        return set(), "none", 0
    mode, max_radius = _phase_field_mode_and_radius(seeds)
    if mode == "point":
        return {candidate.pair for candidate in seeds}, mode, 1
    radius_results: list[
        tuple[int, set[tuple[int, int]], float, int, float]
    ] = []
    for radius in range(1, max_radius + 1):
        support, direct_pairs = _phase_field_support(
            row_length=row_length,
            candidates=all_candidates,
            seeds=seeds,
            radius=radius,
            mode=mode,
        )
        pairs = _phase_field_prefix(support, minimum_count=len(seeds))
        if not pairs:
            continue
        scores = [support[pair] for pair in pairs]
        direct_fraction = len(pairs & direct_pairs) / len(pairs)
        radius_results.append(
            (radius, pairs, direct_fraction, len(pairs), mean(scores))
        )
    if not radius_results:
        return {candidate.pair for candidate in seeds}, mode, 1
    if mode == "diagonal":
        radius, pairs, _, _, _ = max(
            radius_results,
            key=lambda item: (
                item[4],
                item[3],
                item[0],
            ),
        )
        return pairs, mode, radius
    if len(radius_results) == 1:
        radius, pairs, _, _, _ = radius_results[0]
        return pairs, mode, radius
    drops = [
        (
            radius_results[index][2] - radius_results[index + 1][2],
            radius_results[index],
        )
        for index in range(len(radius_results) - 1)
    ]
    _, selected = max(
        drops,
        key=lambda item: (
            item[0],
            item[1][3],
            -item[1][0],
        ),
    )
    radius, pairs, _, _, _ = selected
    return pairs, mode, radius


def _single_gap_phase_ribbon_bridge(
    *,
    row_length: int,
    phase_pairs: set[tuple[int, int]],
) -> set[tuple[int, int]]:
    bridged = set(phase_pairs)
    ordered = sorted(phase_pairs)
    for index, left_pair in enumerate(ordered):
        left_i, left_j = left_pair
        for right_i, right_j in ordered[index + 1 :]:
            delta_i = right_i - left_i
            delta_j = right_j - left_j
            distance = max(abs(delta_i), abs(delta_j))
            if distance != 2 or abs(delta_i - delta_j) > 1:
                continue
            for step in range(distance + 1):
                bridge_i = round(left_i + delta_i * step / distance)
                bridge_j = round(left_j + delta_j * step / distance)
                if (
                    bridge_i >= 1
                    and bridge_j <= row_length
                    and bridge_j - bridge_i >= 24
                ):
                    bridged.add((bridge_i, bridge_j))
    return bridged


def _scaffold_contact_metrics(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
    run: TraceLoopRun,
    exact_metrics: ExactContactMetrics,
) -> ScaffoldContactMetrics:
    constraints_by_row = dataset.constraints_by_row_id()
    selected_by_row: dict[str, list[object]] = {}
    for event in run.selected_events:
        selected_by_row.setdefault(event.row_id, []).append(event)

    precisions: list[float] = []
    long_recalls: list[float] = []
    contact_counts: list[int] = []
    perfect_flags: list[bool] = []
    compact_precisions: list[float] = []
    compact_long_recalls: list[float] = []
    compact_contact_counts: list[int] = []
    compact_perfect_flags: list[bool] = []
    density_compact_precisions: list[float] = []
    density_compact_long_recalls: list[float] = []
    density_compact_contact_counts: list[int] = []
    density_compact_perfect_flags: list[bool] = []
    adjacent_patch_precisions: list[float] = []
    adjacent_patch_long_recalls: list[float] = []
    adjacent_patch_contact_counts: list[int] = []
    adjacent_patch_perfect_flags: list[bool] = []
    phase_field_precisions: list[float] = []
    phase_field_long_recalls: list[float] = []
    phase_field_contact_counts: list[int] = []
    phase_field_perfect_flags: list[bool] = []
    phase_field_modes: list[str] = []
    phase_field_radii: list[int] = []
    phase_ribbon_bridge_precisions: list[float] = []
    phase_ribbon_bridge_long_recalls: list[float] = []
    phase_ribbon_bridge_contact_counts: list[int] = []
    phase_ribbon_bridge_perfect_flags: list[bool] = []
    for row in rows:
        native_pairs = set(row.native_contact_pairs())
        native_long = {pair for pair in native_pairs if pair[1] - pair[0] >= 24}
        selected_events = selected_by_row.get(row.row_id, [])
        candidates_by_pair: dict[tuple[int, int], ScaffoldContactCandidate] = {}
        for constraint in constraints_by_row.get(row.row_id, ()):
            constraint_pair = constraint.pair()
            memberships = [
                event
                for event in selected_events
                if constraint_pair in event.candidate_region_pairs()
            ]
            if not memberships:
                continue
            center_compactness = max(
                _center_compactness_for_event(
                    pair=constraint_pair,
                    event=event,
                )
                for event in memberships
            )
            separation_compactness = max(
                0.0,
                1.0 - constraint.sequence_separation / row.sequence_length,
            )
            candidate = ScaffoldContactCandidate(
                pair=constraint_pair,
                confidence=float(constraint.confidence),
                sequence_separation=int(constraint.sequence_separation),
                event_ids=tuple(sorted(event.event_id for event in memberships)),
                center_compactness=center_compactness,
                separation_compactness=separation_compactness,
            )
            previous = candidates_by_pair.get(constraint_pair)
            if (
                previous is None
                or candidate.compactness_score > previous.compactness_score
            ):
                candidates_by_pair[constraint_pair] = candidate
        scaffold_pairs = set(candidates_by_pair)
        compact_pairs = {
            candidate.pair
            for candidate in _compact_scaffold_core(
                tuple(candidates_by_pair.values())
            )
        }
        density_compact_pairs = {
            candidate.pair
            for candidate in _density_compact_scaffold_core(
                tuple(candidates_by_pair.values())
            )
        }
        adjacent_patch_pairs = {
            candidate.pair
            for candidate in _adjacent_density_patch_scaffold_core(
                tuple(candidates_by_pair.values())
            )
        }
        (
            phase_field_pairs,
            phase_field_mode,
            phase_field_radius,
        ) = _phase_field_scaffold_core(
            row_length=row.sequence_length,
            candidates=tuple(candidates_by_pair.values()),
        )
        phase_ribbon_bridge_pairs = _single_gap_phase_ribbon_bridge(
            row_length=row.sequence_length,
            phase_pairs=phase_field_pairs,
        )
        supported = scaffold_pairs & native_pairs
        supported_long = scaffold_pairs & native_long
        compact_supported = compact_pairs & native_pairs
        compact_supported_long = compact_pairs & native_long
        density_compact_supported = density_compact_pairs & native_pairs
        density_compact_supported_long = density_compact_pairs & native_long
        adjacent_patch_supported = adjacent_patch_pairs & native_pairs
        adjacent_patch_supported_long = adjacent_patch_pairs & native_long
        phase_field_supported = phase_field_pairs & native_pairs
        phase_field_supported_long = phase_field_pairs & native_long
        phase_ribbon_bridge_supported = (
            phase_ribbon_bridge_pairs & native_pairs
        )
        phase_ribbon_bridge_supported_long = (
            phase_ribbon_bridge_pairs & native_long
        )
        precision = (
            _rounded(len(supported) / len(scaffold_pairs))
            if scaffold_pairs
            else 0.0
        )
        long_recall = (
            _rounded(len(supported_long) / len(native_long))
            if native_long
            else 1.0
        )
        precisions.append(precision)
        long_recalls.append(long_recall)
        contact_counts.append(len(scaffold_pairs))
        perfect_flags.append(precision == 1.0 and long_recall == 1.0)
        compact_precision = (
            _rounded(len(compact_supported) / len(compact_pairs))
            if compact_pairs
            else 0.0
        )
        compact_long_recall = (
            _rounded(len(compact_supported_long) / len(native_long))
            if native_long
            else 1.0
        )
        compact_precisions.append(compact_precision)
        compact_long_recalls.append(compact_long_recall)
        compact_contact_counts.append(len(compact_pairs))
        compact_perfect_flags.append(
            compact_precision == 1.0 and compact_long_recall == 1.0
        )
        density_compact_precision = (
            _rounded(len(density_compact_supported) / len(density_compact_pairs))
            if density_compact_pairs
            else 0.0
        )
        density_compact_long_recall = (
            _rounded(len(density_compact_supported_long) / len(native_long))
            if native_long
            else 1.0
        )
        density_compact_precisions.append(density_compact_precision)
        density_compact_long_recalls.append(density_compact_long_recall)
        density_compact_contact_counts.append(len(density_compact_pairs))
        density_compact_perfect_flags.append(
            density_compact_precision == 1.0
            and density_compact_long_recall == 1.0
        )
        adjacent_patch_precision = (
            _rounded(
                len(adjacent_patch_supported) / len(adjacent_patch_pairs)
            )
            if adjacent_patch_pairs
            else 0.0
        )
        adjacent_patch_long_recall = (
            _rounded(
                len(adjacent_patch_supported_long) / len(native_long)
            )
            if native_long
            else 1.0
        )
        adjacent_patch_precisions.append(adjacent_patch_precision)
        adjacent_patch_long_recalls.append(adjacent_patch_long_recall)
        adjacent_patch_contact_counts.append(len(adjacent_patch_pairs))
        adjacent_patch_perfect_flags.append(
            adjacent_patch_precision == 1.0
            and adjacent_patch_long_recall == 1.0
        )
        phase_field_precision = (
            _rounded(len(phase_field_supported) / len(phase_field_pairs))
            if phase_field_pairs
            else 0.0
        )
        phase_field_long_recall = (
            _rounded(len(phase_field_supported_long) / len(native_long))
            if native_long
            else 1.0
        )
        phase_field_precisions.append(phase_field_precision)
        phase_field_long_recalls.append(phase_field_long_recall)
        phase_field_contact_counts.append(len(phase_field_pairs))
        phase_field_perfect_flags.append(
            phase_field_precision == 1.0
            and phase_field_long_recall == 1.0
        )
        phase_field_modes.append(phase_field_mode)
        phase_field_radii.append(phase_field_radius)
        phase_ribbon_bridge_precision = (
            _rounded(
                len(phase_ribbon_bridge_supported)
                / len(phase_ribbon_bridge_pairs)
            )
            if phase_ribbon_bridge_pairs
            else 0.0
        )
        phase_ribbon_bridge_long_recall = (
            _rounded(
                len(phase_ribbon_bridge_supported_long) / len(native_long)
            )
            if native_long
            else 1.0
        )
        phase_ribbon_bridge_precisions.append(phase_ribbon_bridge_precision)
        phase_ribbon_bridge_long_recalls.append(
            phase_ribbon_bridge_long_recall
        )
        phase_ribbon_bridge_contact_counts.append(
            len(phase_ribbon_bridge_pairs)
        )
        phase_ribbon_bridge_perfect_flags.append(
            phase_ribbon_bridge_precision == 1.0
            and phase_ribbon_bridge_long_recall == 1.0
        )

    precision = _rounded(mean(precisions)) if precisions else 0.0
    long_recall = _rounded(mean(long_recalls)) if long_recalls else 0.0
    compact_precision = (
        _rounded(mean(compact_precisions)) if compact_precisions else 0.0
    )
    compact_long_recall = (
        _rounded(mean(compact_long_recalls)) if compact_long_recalls else 0.0
    )
    density_compact_precision = (
        _rounded(mean(density_compact_precisions))
        if density_compact_precisions
        else 0.0
    )
    density_compact_long_recall = (
        _rounded(mean(density_compact_long_recalls))
        if density_compact_long_recalls
        else 0.0
    )
    adjacent_patch_precision = (
        _rounded(mean(adjacent_patch_precisions))
        if adjacent_patch_precisions
        else 0.0
    )
    adjacent_patch_long_recall = (
        _rounded(mean(adjacent_patch_long_recalls))
        if adjacent_patch_long_recalls
        else 0.0
    )
    phase_field_precision = (
        _rounded(mean(phase_field_precisions))
        if phase_field_precisions
        else 0.0
    )
    phase_field_long_recall = (
        _rounded(mean(phase_field_long_recalls))
        if phase_field_long_recalls
        else 0.0
    )
    phase_ribbon_bridge_precision = (
        _rounded(mean(phase_ribbon_bridge_precisions))
        if phase_ribbon_bridge_precisions
        else 0.0
    )
    phase_ribbon_bridge_long_recall = (
        _rounded(mean(phase_ribbon_bridge_long_recalls))
        if phase_ribbon_bridge_long_recalls
        else 0.0
    )
    return ScaffoldContactMetrics(
        scaffold_contact_count=sum(contact_counts),
        scaffold_exact_contact_precision=precision,
        scaffold_exact_long_range_contact_recall=long_recall,
        scaffold_contact_precision_delta_vs_top_L=_rounded(
            precision - exact_metrics.exact_contact_precision_top_L
        ),
        scaffold_contact_recall_delta_vs_region=_rounded(
            long_recall - run.metric.long_range_contact_recall
        ),
        scaffold_contact_map_perfect=bool(perfect_flags) and all(perfect_flags),
        compact_scaffold_contact_count=sum(compact_contact_counts),
        compact_scaffold_exact_contact_precision=compact_precision,
        compact_scaffold_exact_long_range_contact_recall=compact_long_recall,
        compact_scaffold_precision_delta_vs_scaffold=_rounded(
            compact_precision - precision
        ),
        compact_scaffold_recall_delta_vs_scaffold=_rounded(
            compact_long_recall - long_recall
        ),
        compact_scaffold_contact_map_perfect=(
            bool(compact_perfect_flags) and all(compact_perfect_flags)
        ),
        density_compact_scaffold_contact_count=sum(
            density_compact_contact_counts
        ),
        density_compact_scaffold_exact_contact_precision=(
            density_compact_precision
        ),
        density_compact_scaffold_exact_long_range_contact_recall=(
            density_compact_long_recall
        ),
        density_compact_scaffold_precision_delta_vs_compact=_rounded(
            density_compact_precision - compact_precision
        ),
        density_compact_scaffold_recall_delta_vs_compact=_rounded(
            density_compact_long_recall - compact_long_recall
        ),
        density_compact_scaffold_contact_map_perfect=(
            bool(density_compact_perfect_flags)
            and all(density_compact_perfect_flags)
        ),
        adjacent_density_patch_scaffold_contact_count=sum(
            adjacent_patch_contact_counts
        ),
        adjacent_density_patch_scaffold_exact_contact_precision=(
            adjacent_patch_precision
        ),
        adjacent_density_patch_scaffold_exact_long_range_contact_recall=(
            adjacent_patch_long_recall
        ),
        adjacent_density_patch_scaffold_precision_delta_vs_density_compact=(
            _rounded(adjacent_patch_precision - density_compact_precision)
        ),
        adjacent_density_patch_scaffold_recall_delta_vs_density_compact=(
            _rounded(adjacent_patch_long_recall - density_compact_long_recall)
        ),
        adjacent_density_patch_scaffold_contact_map_perfect=(
            bool(adjacent_patch_perfect_flags)
            and all(adjacent_patch_perfect_flags)
        ),
        phase_field_scaffold_contact_count=sum(phase_field_contact_counts),
        phase_field_scaffold_exact_contact_precision=phase_field_precision,
        phase_field_scaffold_exact_long_range_contact_recall=(
            phase_field_long_recall
        ),
        phase_field_scaffold_precision_delta_vs_adjacent_density_patch=(
            _rounded(phase_field_precision - adjacent_patch_precision)
        ),
        phase_field_scaffold_recall_delta_vs_adjacent_density_patch=(
            _rounded(phase_field_long_recall - adjacent_patch_long_recall)
        ),
        phase_field_scaffold_contact_map_perfect=(
            bool(phase_field_perfect_flags)
            and all(phase_field_perfect_flags)
        ),
        phase_field_scaffold_mode=";".join(phase_field_modes),
        phase_field_scaffold_radius=(
            int(round(mean(phase_field_radii))) if phase_field_radii else 0
        ),
        phase_ribbon_bridge_scaffold_contact_count=sum(
            phase_ribbon_bridge_contact_counts
        ),
        phase_ribbon_bridge_scaffold_exact_contact_precision=(
            phase_ribbon_bridge_precision
        ),
        phase_ribbon_bridge_scaffold_exact_long_range_contact_recall=(
            phase_ribbon_bridge_long_recall
        ),
        phase_ribbon_bridge_scaffold_precision_delta_vs_phase_field=(
            _rounded(phase_ribbon_bridge_precision - phase_field_precision)
        ),
        phase_ribbon_bridge_scaffold_recall_delta_vs_phase_field=(
            _rounded(phase_ribbon_bridge_long_recall - phase_field_long_recall)
        ),
        phase_ribbon_bridge_scaffold_contact_map_perfect=(
            bool(phase_ribbon_bridge_perfect_flags)
            and all(phase_ribbon_bridge_perfect_flags)
        ),
    )


def _metric_row(
    *,
    target_id: str,
    run: TraceLoopRun,
    exact_metrics: ExactContactMetrics,
    scaffold_metrics: ScaffoldContactMetrics,
    benchmark_file: Path,
    external_coupling_file: Path,
    benchmark_file_sha256: str,
    external_coupling_file_sha256: str,
    source_accessions: Sequence[str],
    control_name: str = "",
) -> dict[str, object]:
    metric = run.metric
    false_event_count = round(
        metric.false_nucleus_rate * metric.selected_event_count
    )
    row = {
        "target_id": target_id,
        "selector_name": run.selector_name,
        "control_name": control_name,
        "control_kind": run.control_kind,
        "benchmark_file": str(benchmark_file),
        "external_coupling_file": str(external_coupling_file),
        "benchmark_file_sha256": benchmark_file_sha256,
        "external_coupling_file_sha256": external_coupling_file_sha256,
        "source_accessions": ";".join(source_accessions),
        "constraint_count": run.constraint_count,
        "selected_event_count": metric.selected_event_count,
        "false_event_count": false_event_count,
        "false_nucleus_rate": metric.false_nucleus_rate,
        "cluster_precision": metric.contact_cluster_precision,
        "long_range_contact_recall": metric.long_range_contact_recall,
        "coupling_constraint_recall": metric.coupling_constraint_recall,
        "real_vs_decoy_coupling_enrichment_ratio": (
            metric.real_vs_decoy_coupling_enrichment_ratio
        ),
        "coordinate_truth_used_to_build_constraints": (
            metric.coordinate_truth_used_to_build_constraints
        ),
        "native_truth_used_before_coupling_selection": (
            metric.native_truth_used_before_coupling_selection
        ),
        "raw_sequence_exposed": metric.raw_sequence_exposed,
    }
    row.update(exact_metrics.to_dict())
    row.update(scaffold_metrics.to_dict())
    return row


def _run_real_frontier(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
    physical_contexts: Mapping[int, object],
    control_kind: str,
) -> tuple[TraceLoopRun, TraceLoopRun]:
    precision = _run_multiscale_phase_aligned_critical_boundary_selector(
        rows=rows,
        dataset=dataset,
        selector_name=PRECISION_BOUNDARY_SELECTOR,
        control_kind=control_kind,
        physical_contexts=physical_contexts,
    )
    if SATURATION_BOUNDARY_SELECTOR == FOOTPRINT_NOVELTY_SELECTOR:
        saturation = (
            _run_multiscale_phase_aligned_footprint_novelty_boundary_selector(
                rows=rows,
                dataset=dataset,
                selector_name=SATURATION_BOUNDARY_SELECTOR,
                control_kind=control_kind,
                physical_contexts=physical_contexts,
            )
        )
    else:
        saturation = _run_multiscale_phase_aligned_external_novelty_boundary_selector(
            rows=rows,
            dataset=dataset,
            selector_name=EXTERNAL_NOVELTY_SELECTOR,
            control_kind=control_kind,
            physical_contexts=physical_contexts,
        )
    return precision, saturation


def _run_control_frontier(
    *,
    rows: Sequence[RealCoordinateVisualRow],
    dataset: CouplingDataset,
    physical_contexts: Mapping[int, object],
    control_name: str,
    control_kind: str,
) -> tuple[TraceLoopRun, TraceLoopRun]:
    precision = _run_multiscale_phase_aligned_critical_boundary_selector(
        rows=rows,
        dataset=dataset,
        selector_name=f"{control_name}_precision_boundary",
        control_kind=control_kind,
        physical_contexts=physical_contexts,
    )
    if SATURATION_BOUNDARY_SELECTOR == FOOTPRINT_NOVELTY_SELECTOR:
        saturation = (
            _run_multiscale_phase_aligned_footprint_novelty_boundary_selector(
                rows=rows,
                dataset=dataset,
                selector_name=f"{control_name}_footprint_novelty_boundary",
                control_kind=control_kind,
                physical_contexts=physical_contexts,
            )
        )
    else:
        saturation = _run_multiscale_phase_aligned_external_novelty_boundary_selector(
            rows=rows,
            dataset=dataset,
            selector_name=f"{control_name}_external_novelty_boundary",
            control_kind=control_kind,
            physical_contexts=physical_contexts,
        )
    return precision, saturation


def _dominates(left: TraceLoopRun, right: TraceLoopRun) -> bool:
    left_metric = left.metric
    right_metric = right.metric
    return (
        left_metric.false_nucleus_rate <= right_metric.false_nucleus_rate
        and left_metric.long_range_contact_recall
        >= right_metric.long_range_contact_recall
        and left_metric.contact_cluster_precision
        > right_metric.contact_cluster_precision
    )


def _frontier_dominates(
    left_runs: Sequence[TraceLoopRun],
    right_runs: Sequence[TraceLoopRun],
) -> bool:
    return all(
        any(_dominates(left, right) for left in left_runs)
        for right in right_runs
    )


def _run_is_productive(run: TraceLoopRun) -> bool:
    return (
        run.metric.selected_event_count > 0
        and run.metric.long_range_contact_recall > 0.0
        and run.metric.contact_cluster_precision > 0.0
    )


def _productive_dominates(
    *,
    left: TraceLoopRun,
    left_exact: ExactContactMetrics,
    right: TraceLoopRun,
    right_exact: ExactContactMetrics,
    pair_randomizing_control: bool,
) -> bool:
    left_metric = left.metric
    right_metric = right.metric
    if not _run_is_productive(left):
        return False
    if right_metric.selected_event_count == 0:
        return left_metric.false_nucleus_rate <= right_metric.false_nucleus_rate
    if _dominates(left, right):
        return True
    if (
        left_metric.false_nucleus_rate < right_metric.false_nucleus_rate
        and left_metric.long_range_contact_recall
        >= right_metric.long_range_contact_recall
    ):
        return True
    if not pair_randomizing_control:
        return False
    return (
        left_metric.false_nucleus_rate <= right_metric.false_nucleus_rate
        and left_exact.exact_contact_precision_top_L
        > right_exact.exact_contact_precision_top_L
        and left_exact.exact_long_range_contact_precision
        > right_exact.exact_long_range_contact_precision
    )


def _productive_frontier_dominates(
    *,
    left_runs: Sequence[tuple[TraceLoopRun, ExactContactMetrics]],
    right_runs: Sequence[tuple[TraceLoopRun, ExactContactMetrics]],
    control_name: str,
) -> bool:
    pair_randomizing_control = control_name in PAIR_RANDOMIZING_CONTROL_NAMES
    return all(
        any(
            _productive_dominates(
                left=left_run,
                left_exact=left_exact,
                right=right_run,
                right_exact=right_exact,
                pair_randomizing_control=pair_randomizing_control,
            )
            for left_run, left_exact in left_runs
        )
        for right_run, right_exact in right_runs
    )


def _frontier_best_metric(
    runs: Sequence[TraceLoopRun],
    field_name: str,
) -> float:
    return _rounded(max(float(getattr(run.metric, field_name)) for run in runs))


def _write_json(path: Path, payload: Mapping[str, object]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def _control_names(raw: Sequence[str]) -> tuple[str, ...]:
    if not raw or raw == ("all",):
        return tuple(EXTERNAL_COUPLING_CONTROL_NAMES)
    names = tuple(raw)
    unknown = sorted(set(names).difference(EXTERNAL_COUPLING_CONTROL_NAMES))
    if unknown:
        raise ValueError(f"unknown control name(s): {', '.join(unknown)}")
    return names


def _mean_field(rows: Sequence[Mapping[str, object]], field_name: str) -> float:
    values = [float(row[field_name]) for row in rows]
    return _rounded(mean(values)) if values else 0.0


def _max_field(rows: Sequence[Mapping[str, object]], field_name: str) -> float:
    values = [float(row[field_name]) for row in rows]
    return _rounded(max(values)) if values else 0.0


def _pair_randomizing_control_rows(
    rows: Sequence[Mapping[str, object]],
) -> tuple[Mapping[str, object], ...]:
    return tuple(
        row
        for row in rows
        if str(row.get("control_name", "")) in PAIR_RANDOMIZING_CONTROL_NAMES
    )


def _battery_classification(
    *,
    target_rows: Sequence[Mapping[str, object]],
    control_rows: Sequence[Mapping[str, object]],
    failure_rows: Sequence[Mapping[str, object]],
    oracle_taint_violation_count: int,
    target_manifest_frozen: bool,
    selector_frozen_after_manifest: bool,
) -> str:
    if (
        oracle_taint_violation_count
        or not target_manifest_frozen
        or not selector_frozen_after_manifest
    ):
        return "folding_claim_forbidden"
    if failure_rows:
        return "blind_batch_failed"
    if not target_rows or not control_rows:
        return "blind_batch_inconclusive"

    frontier_win_rate = _mean_field(target_rows, "productive_frontier_control_win")
    pair_randomizing_rows = _pair_randomizing_control_rows(control_rows)
    exact_control_rows = pair_randomizing_rows or tuple(control_rows)
    scaffold_signal = frontier_win_rate >= PASS_TARGET_WIN_RATE_MIN
    exact_signal = (
        scaffold_signal
        and _mean_field(target_rows, "exact_contact_precision_top_L")
        > _max_field(exact_control_rows, "exact_contact_precision_top_L")
        and _mean_field(target_rows, "exact_long_range_contact_precision")
        > _max_field(exact_control_rows, "exact_long_range_contact_precision")
    )
    if exact_signal:
        return "blind_batch_exact_contact_signal_confirmed"
    if scaffold_signal:
        return "blind_batch_scaffold_signal_confirmed"
    return "blind_batch_failed"


def _folding_problem_solved_audit(
    *,
    manifest: Mapping[str, object],
    target_rows: Sequence[Mapping[str, object]],
    failure_rows: Sequence[Mapping[str, object]],
    hard_gates: Mapping[str, object],
    classification: str,
) -> dict[str, object]:
    requirements = {
        "exact_contact_signal_confirmed": (
            classification == "blind_batch_exact_contact_signal_confirmed"
        ),
        "all_hard_gates_passed": bool(hard_gates["all_hard_gates_passed"]),
        "no_target_failures": not failure_rows,
        "all_targets_zero_false_nucleus_rate": (
            bool(target_rows)
            and all(float(row["false_nucleus_rate"]) == 0.0 for row in target_rows)
        ),
        "all_targets_productive_frontier_win": (
            bool(target_rows)
            and all(
                int(row.get("productive_frontier_control_win", 0)) == 1
                for row in target_rows
            )
        ),
        "all_targets_scaffold_contact_maps_perfect": (
            bool(target_rows)
            and all(
                str(row.get("scaffold_contact_map_perfect", "False")) == "True"
                or row.get("scaffold_contact_map_perfect") is True
                for row in target_rows
            )
        ),
        "all_targets_compact_scaffold_contact_maps_perfect": (
            bool(target_rows)
            and all(
                str(row.get("compact_scaffold_contact_map_perfect", "False"))
                == "True"
                or row.get("compact_scaffold_contact_map_perfect") is True
                for row in target_rows
            )
        ),
        "all_targets_density_compact_scaffold_contact_maps_perfect": (
            bool(target_rows)
            and all(
                str(
                    row.get(
                        "density_compact_scaffold_contact_map_perfect",
                        "False",
                    )
                )
                == "True"
                or row.get("density_compact_scaffold_contact_map_perfect") is True
                for row in target_rows
            )
        ),
        "all_targets_adjacent_density_patch_scaffold_contact_maps_perfect": (
            bool(target_rows)
            and all(
                str(
                    row.get(
                        "adjacent_density_patch_scaffold_contact_map_perfect",
                        "False",
                    )
                )
                == "True"
                or row.get(
                    "adjacent_density_patch_scaffold_contact_map_perfect"
                )
                is True
                for row in target_rows
            )
        ),
        "all_targets_phase_field_scaffold_contact_maps_perfect": (
            bool(target_rows)
            and all(
                str(
                    row.get(
                        "phase_field_scaffold_contact_map_perfect",
                        "False",
                    )
                )
                == "True"
                or row.get("phase_field_scaffold_contact_map_perfect") is True
                for row in target_rows
            )
        ),
        "all_targets_phase_ribbon_bridge_scaffold_contact_maps_perfect": (
            bool(target_rows)
            and all(
                str(
                    row.get(
                        "phase_ribbon_bridge_scaffold_contact_map_perfect",
                        "False",
                    )
                )
                == "True"
                or row.get(
                    "phase_ribbon_bridge_scaffold_contact_map_perfect"
                )
                is True
                for row in target_rows
            )
        ),
        "universal_holdout_scope_declared": (
            str(manifest.get("holdout_scope", ""))
            == "all_available_complete_external_coupling_protein_targets"
        ),
        "all_existing_complete_targets_tested": bool(
            manifest.get("all_existing_complete_targets_tested", False)
        ),
        "full_atomic_folding_available": bool(
            manifest.get("full_atomic_folding_available", False)
        ),
        "atomic_coordinate_generation_validated": bool(
            manifest.get("atomic_coordinate_generation_validated", False)
        ),
    }
    unmet = tuple(
        name for name, passed in requirements.items() if not passed
    )
    return {
        "audit_kind": "folding_problem_solved_audit_v0",
        "solved": not unmet,
        "requirements": requirements,
        "unmet_requirements": unmet,
        "current_claim_scope": (
            "long_range_folding_nucleus_region_recovery_only_not_full_atomic_folding"
        ),
    }


def run_blind_external_holdout_battery_v0(
    *,
    target_manifest: Path,
    coupling_dir: Path,
    output_dir: Path,
    selector_name: str,
    control_names: Sequence[str],
    target_ids: Sequence[str] = (),
) -> tuple[Path, Path, Path, Path, Path, Path]:
    if selector_name not in (
        DEFAULT_SELECTOR,
        EXTERNAL_NOVELTY_SELECTOR,
        FRONTIER_SELECTOR,
    ):
        raise ValueError(f"unsupported frozen selector: {selector_name}")

    manifest = _load_json(target_manifest)
    targets = manifest.get("targets", [])
    if not isinstance(targets, list) or not targets:
        raise ValueError("target manifest must include a non-empty targets list")
    requested_target_ids = {target_id.strip() for target_id in target_ids if target_id.strip()}
    if requested_target_ids:
        targets = [
            target
            for target in targets
            if isinstance(target, Mapping)
            and str(target.get("target_id", "")).strip() in requested_target_ids
        ]
        if not targets:
            raise ValueError(
                "target-id filter did not match any manifest targets: "
                + ", ".join(sorted(requested_target_ids))
            )
    controls_to_run = _control_names(tuple(control_names))
    target_manifest_sha256 = _sha256_file(target_manifest)
    current_commit = _current_commit_hash()
    selector_frozen_after_manifest = bool(
        manifest.get("selector_frozen_after_manifest", False)
    )

    target_rows: list[dict[str, object]] = []
    control_rows: list[dict[str, object]] = []
    selected_event_rows: list[dict[str, object]] = []
    failure_rows: list[dict[str, object]] = []
    coordinate_taint_violation_count = 0
    native_taint_violation_count = 0
    oracle_taint_violation_count = 0

    for raw_target in targets:
        if not isinstance(raw_target, Mapping):
            failure_rows.append({"target_id": "", "failure_kind": "invalid_target"})
            continue
        target_id = str(raw_target.get("target_id", "")).strip()
        if not target_id:
            failure_rows.append({"target_id": "", "failure_kind": "missing_target_id"})
            continue
        benchmark_file = _resolve_path(
            raw_target.get("benchmark_file", ""),
            base_dir=REPO_ROOT,
        )
        external_coupling_file = _target_coupling_file(
            raw_target,
            coupling_dir=coupling_dir,
        )
        if not external_coupling_file.is_absolute():
            external_coupling_file = (REPO_ROOT / external_coupling_file).resolve()

        try:
            rows = load_real_coordinate_visual_rows(benchmark_file)
            import_result = import_external_coupling_dataset(
                rows=rows,
                external_coupling_file=external_coupling_file,
            )
            if import_result.dataset.coordinate_truth_tainted:
                coordinate_taint_violation_count += 1
            if import_result.dataset.native_truth_tainted:
                native_taint_violation_count += 1
            if import_result.dataset.oracle_constraint_control:
                oracle_taint_violation_count += 1
            if not import_result.dataset.external_evolutionary_couplings_used:
                failure_rows.append(
                    {
                        "target_id": target_id,
                        "failure_kind": "external_couplings_not_used",
                    }
                )
                continue
            physical_contexts = _build_multiscale_physical_contexts(rows)
            real_precision_run, real_saturation_run = _run_real_frontier(
                rows=rows,
                dataset=import_result.dataset,
                control_kind="external_real_blind_holdout",
                physical_contexts=physical_contexts,
            )
            controls = generate_external_coupling_negative_controls(
                rows=rows,
                dataset=import_result.dataset,
            )
            benchmark_file_sha256 = _sha256_file(benchmark_file)
            coupling_file_sha256 = _sha256_file(external_coupling_file)
            source_accessions = tuple(row.source_accession for row in rows)
            real_exact_metrics = _exact_contact_metrics(
                rows=rows,
                dataset=import_result.dataset,
                run=real_saturation_run,
            )
            real_scaffold_metrics = _scaffold_contact_metrics(
                rows=rows,
                dataset=import_result.dataset,
                run=real_saturation_run,
                exact_metrics=real_exact_metrics,
            )
            real_precision_exact_metrics = _exact_contact_metrics(
                rows=rows,
                dataset=import_result.dataset,
                run=real_precision_run,
            )
            target_row = _metric_row(
                target_id=target_id,
                run=real_saturation_run,
                exact_metrics=real_exact_metrics,
                scaffold_metrics=real_scaffold_metrics,
                benchmark_file=benchmark_file,
                external_coupling_file=external_coupling_file,
                benchmark_file_sha256=benchmark_file_sha256,
                external_coupling_file_sha256=coupling_file_sha256,
                source_accessions=source_accessions,
            )
            target_row["frontier_selector_name"] = FRONTIER_SELECTOR
            target_row["precision_boundary_selected_event_count"] = (
                real_precision_run.metric.selected_event_count
            )
            target_row["precision_boundary_false_nucleus_rate"] = (
                real_precision_run.metric.false_nucleus_rate
            )
            target_row["precision_boundary_cluster_precision"] = (
                real_precision_run.metric.contact_cluster_precision
            )
            target_row["precision_boundary_long_range_contact_recall"] = (
                real_precision_run.metric.long_range_contact_recall
            )
            target_row["precision_boundary_exact_contact_precision_top_L"] = (
                real_precision_exact_metrics.exact_contact_precision_top_L
            )
            target_row["precision_boundary_region_width_penalty"] = (
                real_precision_exact_metrics.region_width_penalty
            )
            row_control_rows: list[dict[str, object]] = []
            control_frontier_win_flags: list[int] = []
            productive_frontier_win_flags: list[int] = []
            row_controls_to_run = tuple(
                control_name
                for control_name in controls_to_run
                if not (
                    control_name == "external_cross_row_swapped"
                    and len(rows) < 2
                )
            )
            for control_name in row_controls_to_run:
                control = controls[control_name]
                control_precision_run, control_saturation_run = _run_control_frontier(
                    rows=rows,
                    dataset=control.dataset,
                    control_name=control_name,
                    control_kind=control.control_kind,
                    physical_contexts=physical_contexts,
                )
                control_frontier_win_flags.append(
                    int(
                        _frontier_dominates(
                            (real_precision_run, real_saturation_run),
                            (control_precision_run, control_saturation_run),
                        )
                    )
                )
                control_precision_exact_metrics = _exact_contact_metrics(
                    rows=rows,
                    dataset=control.dataset,
                    run=control_precision_run,
                )
                control_saturation_exact_metrics = _exact_contact_metrics(
                    rows=rows,
                    dataset=control.dataset,
                    run=control_saturation_run,
                )
                productive_frontier_win_flags.append(
                    int(
                        _productive_frontier_dominates(
                            left_runs=(
                                (
                                    real_precision_run,
                                    real_precision_exact_metrics,
                                ),
                                (real_saturation_run, real_exact_metrics),
                            ),
                            right_runs=(
                                (
                                    control_precision_run,
                                    control_precision_exact_metrics,
                                ),
                                (
                                    control_saturation_run,
                                    control_saturation_exact_metrics,
                                ),
                            ),
                            control_name=control_name,
                        )
                    )
                )
                for boundary_name, control_run, control_exact_metrics in (
                    (
                        "precision_boundary",
                        control_precision_run,
                        control_precision_exact_metrics,
                    ),
                    (
                        "saturation_boundary",
                        control_saturation_run,
                        control_saturation_exact_metrics,
                    ),
                ):
                    control_row = _metric_row(
                        target_id=target_id,
                        run=control_run,
                        exact_metrics=control_exact_metrics,
                        scaffold_metrics=_scaffold_contact_metrics(
                            rows=rows,
                            dataset=control.dataset,
                            run=control_run,
                            exact_metrics=control_exact_metrics,
                        ),
                        benchmark_file=benchmark_file,
                        external_coupling_file=external_coupling_file,
                        benchmark_file_sha256=benchmark_file_sha256,
                        external_coupling_file_sha256=coupling_file_sha256,
                        source_accessions=source_accessions,
                        control_name=control_name,
                    )
                    control_row["frontier_member"] = boundary_name
                    row_control_rows.append(control_row)
                    control_rows.append(control_row)

            target_row["max_control_false_nucleus_rate"] = _max_field(
                row_control_rows,
                "false_nucleus_rate",
            )
            target_row["max_control_long_range_contact_recall"] = _max_field(
                row_control_rows,
                "long_range_contact_recall",
            )
            target_row["max_control_cluster_precision"] = _max_field(
                row_control_rows,
                "cluster_precision",
            )
            target_row["frontier_control_count"] = len(control_frontier_win_flags)
            target_row["frontier_beaten_control_count"] = sum(
                control_frontier_win_flags
            )
            target_row["frontier_control_win"] = int(
                bool(control_frontier_win_flags)
                and all(control_frontier_win_flags)
            )
            target_row["productive_frontier_beaten_control_count"] = sum(
                productive_frontier_win_flags
            )
            target_row["productive_frontier_control_win"] = int(
                bool(productive_frontier_win_flags)
                and all(productive_frontier_win_flags)
            )
            target_row["individual_control_win"] = target_row[
                "productive_frontier_control_win"
            ]
            target_rows.append(target_row)
            for real_run in (real_precision_run, real_saturation_run):
                for selected_row in real_run.selected_rows:
                    selected_event_rows.append(
                        {
                            "target_id": target_id,
                            "frontier_selector_name": FRONTIER_SELECTOR,
                            "benchmark_file_sha256": benchmark_file_sha256,
                            "external_coupling_file_sha256": coupling_file_sha256,
                            **selected_row,
                        }
                    )
        except Exception as exc:  # noqa: BLE001
            failure_rows.append(
                {
                    "target_id": target_id,
                    "failure_kind": exc.__class__.__name__,
                    "failure_message": str(exc),
                    "benchmark_file": str(benchmark_file),
                    "external_coupling_file": str(external_coupling_file),
                }
            )

    classification = _battery_classification(
        target_rows=target_rows,
        control_rows=control_rows,
        failure_rows=failure_rows,
        oracle_taint_violation_count=oracle_taint_violation_count,
        target_manifest_frozen=bool(manifest.get("manifest_frozen", False)),
        selector_frozen_after_manifest=selector_frozen_after_manifest,
    )
    hard_gates = {
        "coordinate_truth_used_to_build_constraints": (
            coordinate_taint_violation_count > 0
        ),
        "native_truth_used_before_coupling_selection": (
            native_taint_violation_count > 0
        ),
        "oracle_constraint_control": oracle_taint_violation_count > 0,
        "external_evolutionary_couplings_used": not failure_rows,
        "target_manifest_frozen": bool(manifest.get("manifest_frozen", False)),
        "selector_frozen_after_manifest": selector_frozen_after_manifest,
    }
    hard_gates["all_hard_gates_passed"] = (
        not hard_gates["coordinate_truth_used_to_build_constraints"]
        and not hard_gates["native_truth_used_before_coupling_selection"]
        and not hard_gates["oracle_constraint_control"]
        and bool(hard_gates["external_evolutionary_couplings_used"])
        and bool(hard_gates["target_manifest_frozen"])
        and bool(hard_gates["selector_frozen_after_manifest"])
    )
    pair_randomizing_control_rows = _pair_randomizing_control_rows(control_rows)
    folding_solved_audit = _folding_problem_solved_audit(
        manifest=manifest,
        target_rows=target_rows,
        failure_rows=failure_rows,
        hard_gates=hard_gates,
        classification=classification,
    )
    report = {
        "report_kind": BATTERY_REPORT_KIND,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "commit_hash": current_commit,
        "selector_name": FRONTIER_SELECTOR,
        "saturation_selector_name": SATURATION_BOUNDARY_SELECTOR,
        "precision_boundary_selector_name": PRECISION_BOUNDARY_SELECTOR,
        "target_manifest": str(target_manifest),
        "target_manifest_sha256": target_manifest_sha256,
        "target_manifest_declared_kind": str(manifest.get("manifest_kind", "")),
        "target_manifest_frozen": bool(manifest.get("manifest_frozen", False)),
        "selector_frozen_after_manifest": selector_frozen_after_manifest,
        "control_names": tuple(controls_to_run),
        "target_id_filter": tuple(sorted(requested_target_ids)),
        "target_count": len(target_rows),
        "failure_count": len(failure_rows),
        "coordinate_taint_violation_count": coordinate_taint_violation_count,
        "native_taint_violation_count": native_taint_violation_count,
        "oracle_taint_violation_count": oracle_taint_violation_count,
        "mean_false_nucleus_rate": _mean_field(target_rows, "false_nucleus_rate"),
        "mean_long_range_contact_recall": _mean_field(
            target_rows,
            "long_range_contact_recall",
        ),
        "mean_cluster_precision": _mean_field(target_rows, "cluster_precision"),
        "mean_exact_contact_precision_top_L": _mean_field(
            target_rows,
            "exact_contact_precision_top_L",
        ),
        "mean_exact_contact_precision_L_over_2": _mean_field(
            target_rows,
            "exact_contact_precision_L_over_2",
        ),
        "mean_exact_long_range_contact_precision": _mean_field(
            target_rows,
            "exact_long_range_contact_precision",
        ),
        "mean_region_width_penalty": _mean_field(
            target_rows,
            "region_width_penalty",
        ),
        "mean_scaffold_contact_count": _mean_field(
            target_rows,
            "scaffold_contact_count",
        ),
        "mean_scaffold_exact_contact_precision": _mean_field(
            target_rows,
            "scaffold_exact_contact_precision",
        ),
        "mean_scaffold_exact_long_range_contact_recall": _mean_field(
            target_rows,
            "scaffold_exact_long_range_contact_recall",
        ),
        "mean_scaffold_contact_precision_delta_vs_top_L": _mean_field(
            target_rows,
            "scaffold_contact_precision_delta_vs_top_L",
        ),
        "mean_scaffold_contact_recall_delta_vs_region": _mean_field(
            target_rows,
            "scaffold_contact_recall_delta_vs_region",
        ),
        "scaffold_contact_map_perfect_rate": _mean_field(
            target_rows,
            "scaffold_contact_map_perfect",
        ),
        "mean_compact_scaffold_contact_count": _mean_field(
            target_rows,
            "compact_scaffold_contact_count",
        ),
        "mean_compact_scaffold_exact_contact_precision": _mean_field(
            target_rows,
            "compact_scaffold_exact_contact_precision",
        ),
        "mean_compact_scaffold_exact_long_range_contact_recall": _mean_field(
            target_rows,
            "compact_scaffold_exact_long_range_contact_recall",
        ),
        "mean_compact_scaffold_precision_delta_vs_scaffold": _mean_field(
            target_rows,
            "compact_scaffold_precision_delta_vs_scaffold",
        ),
        "mean_compact_scaffold_recall_delta_vs_scaffold": _mean_field(
            target_rows,
            "compact_scaffold_recall_delta_vs_scaffold",
        ),
        "compact_scaffold_contact_map_perfect_rate": _mean_field(
            target_rows,
            "compact_scaffold_contact_map_perfect",
        ),
        "mean_density_compact_scaffold_contact_count": _mean_field(
            target_rows,
            "density_compact_scaffold_contact_count",
        ),
        "mean_density_compact_scaffold_exact_contact_precision": _mean_field(
            target_rows,
            "density_compact_scaffold_exact_contact_precision",
        ),
        "mean_density_compact_scaffold_exact_long_range_contact_recall": (
            _mean_field(
                target_rows,
                "density_compact_scaffold_exact_long_range_contact_recall",
            )
        ),
        "mean_density_compact_scaffold_precision_delta_vs_compact": (
            _mean_field(
                target_rows,
                "density_compact_scaffold_precision_delta_vs_compact",
            )
        ),
        "mean_density_compact_scaffold_recall_delta_vs_compact": (
            _mean_field(
                target_rows,
                "density_compact_scaffold_recall_delta_vs_compact",
            )
        ),
        "density_compact_scaffold_contact_map_perfect_rate": _mean_field(
            target_rows,
            "density_compact_scaffold_contact_map_perfect",
        ),
        "mean_adjacent_density_patch_scaffold_contact_count": _mean_field(
            target_rows,
            "adjacent_density_patch_scaffold_contact_count",
        ),
        "mean_adjacent_density_patch_scaffold_exact_contact_precision": (
            _mean_field(
                target_rows,
                "adjacent_density_patch_scaffold_exact_contact_precision",
            )
        ),
        "mean_adjacent_density_patch_scaffold_exact_long_range_contact_recall": (
            _mean_field(
                target_rows,
                "adjacent_density_patch_scaffold_exact_long_range_contact_recall",
            )
        ),
        "mean_adjacent_density_patch_scaffold_precision_delta_vs_density_compact": (
            _mean_field(
                target_rows,
                (
                    "adjacent_density_patch_scaffold_precision_delta_vs_"
                    "density_compact"
                ),
            )
        ),
        "mean_adjacent_density_patch_scaffold_recall_delta_vs_density_compact": (
            _mean_field(
                target_rows,
                (
                    "adjacent_density_patch_scaffold_recall_delta_vs_"
                    "density_compact"
                ),
            )
        ),
        "adjacent_density_patch_scaffold_contact_map_perfect_rate": (
            _mean_field(
                target_rows,
                "adjacent_density_patch_scaffold_contact_map_perfect",
            )
        ),
        "mean_phase_field_scaffold_contact_count": _mean_field(
            target_rows,
            "phase_field_scaffold_contact_count",
        ),
        "mean_phase_field_scaffold_exact_contact_precision": _mean_field(
            target_rows,
            "phase_field_scaffold_exact_contact_precision",
        ),
        "mean_phase_field_scaffold_exact_long_range_contact_recall": (
            _mean_field(
                target_rows,
                "phase_field_scaffold_exact_long_range_contact_recall",
            )
        ),
        "mean_phase_field_scaffold_precision_delta_vs_adjacent_density_patch": (
            _mean_field(
                target_rows,
                (
                    "phase_field_scaffold_precision_delta_vs_"
                    "adjacent_density_patch"
                ),
            )
        ),
        "mean_phase_field_scaffold_recall_delta_vs_adjacent_density_patch": (
            _mean_field(
                target_rows,
                (
                    "phase_field_scaffold_recall_delta_vs_"
                    "adjacent_density_patch"
                ),
            )
        ),
        "phase_field_scaffold_contact_map_perfect_rate": _mean_field(
            target_rows,
            "phase_field_scaffold_contact_map_perfect",
        ),
        "mean_phase_ribbon_bridge_scaffold_contact_count": _mean_field(
            target_rows,
            "phase_ribbon_bridge_scaffold_contact_count",
        ),
        "mean_phase_ribbon_bridge_scaffold_exact_contact_precision": (
            _mean_field(
                target_rows,
                "phase_ribbon_bridge_scaffold_exact_contact_precision",
            )
        ),
        "mean_phase_ribbon_bridge_scaffold_exact_long_range_contact_recall": (
            _mean_field(
                target_rows,
                "phase_ribbon_bridge_scaffold_exact_long_range_contact_recall",
            )
        ),
        "mean_phase_ribbon_bridge_scaffold_precision_delta_vs_phase_field": (
            _mean_field(
                target_rows,
                (
                    "phase_ribbon_bridge_scaffold_precision_delta_vs_"
                    "phase_field"
                ),
            )
        ),
        "mean_phase_ribbon_bridge_scaffold_recall_delta_vs_phase_field": (
            _mean_field(
                target_rows,
                (
                    "phase_ribbon_bridge_scaffold_recall_delta_vs_"
                    "phase_field"
                ),
            )
        ),
        "phase_ribbon_bridge_scaffold_contact_map_perfect_rate": _mean_field(
            target_rows,
            "phase_ribbon_bridge_scaffold_contact_map_perfect",
        ),
        "frontier_control_win_rate": _mean_field(
            target_rows,
            "frontier_control_win",
        ),
        "productive_frontier_control_win_rate": _mean_field(
            target_rows,
            "productive_frontier_control_win",
        ),
        "individual_control_win_rate": _mean_field(
            target_rows,
            "individual_control_win",
        ),
        "max_control_false_nucleus_rate": _max_field(
            control_rows,
            "false_nucleus_rate",
        ),
        "max_control_long_range_contact_recall": _max_field(
            control_rows,
            "long_range_contact_recall",
        ),
        "max_control_cluster_precision": _max_field(
            control_rows,
            "cluster_precision",
        ),
        "max_control_exact_contact_precision_top_L": _max_field(
            control_rows,
            "exact_contact_precision_top_L",
        ),
        "max_control_exact_long_range_contact_precision": _max_field(
            control_rows,
            "exact_long_range_contact_precision",
        ),
        "pair_randomizing_control_count": len(pair_randomizing_control_rows),
        "pair_randomizing_control_names": tuple(
            sorted(
                {
                    str(row.get("control_name", ""))
                    for row in pair_randomizing_control_rows
                }
            )
        ),
        "max_pair_randomizing_control_exact_contact_precision_top_L": (
            _max_field(
                pair_randomizing_control_rows,
                "exact_contact_precision_top_L",
            )
        ),
        "max_pair_randomizing_control_exact_long_range_contact_precision": (
            _max_field(
                pair_randomizing_control_rows,
                "exact_long_range_contact_precision",
            )
        ),
        "hard_gates": hard_gates,
        "final_classification": classification,
        "folding_problem_solved": bool(folding_solved_audit["solved"]),
        "folding_problem_solved_audit": folding_solved_audit,
        "claim_allowed": classification
        in (
            "blind_batch_scaffold_signal_confirmed",
            "blind_batch_exact_contact_signal_confirmed",
        )
        and bool(hard_gates["all_hard_gates_passed"]),
    }
    certificate = {
        "certificate_kind": BATTERY_CERTIFICATE_KIND,
        "report_kind": BATTERY_REPORT_KIND,
        "commit_hash": current_commit,
        "selector_name": FRONTIER_SELECTOR,
        "saturation_selector_name": SATURATION_BOUNDARY_SELECTOR,
        "precision_boundary_selector_name": PRECISION_BOUNDARY_SELECTOR,
        "target_manifest_sha256": target_manifest_sha256,
        "final_classification": classification,
        "folding_problem_solved": report["folding_problem_solved"],
        "folding_problem_solved_audit": folding_solved_audit,
        "claim_allowed": report["claim_allowed"],
        "claim_scope": (
            "long_range_folding_nucleus_region_recovery_only_not_full_atomic_folding"
        ),
        "hard_gates": report["hard_gates"],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = _write_json(output_dir / "blind_external_holdout_report.json", report)
    rows_path = write_csv_rows(
        target_rows,
        output_dir / "blind_external_holdout_rows.csv",
    )
    selected_events_path = write_csv_rows(
        selected_event_rows,
        output_dir / "blind_external_holdout_selected_events.csv",
    )
    controls_path = write_csv_rows(
        control_rows,
        output_dir / "blind_external_holdout_controls.csv",
    )
    failures_path = write_csv_rows(
        failure_rows,
        output_dir / "blind_external_holdout_failures.csv",
    )
    certificate_path = _write_json(
        output_dir / "blind_external_holdout_certificate.json",
        certificate,
    )
    return (
        report_path,
        rows_path,
        selected_events_path,
        controls_path,
        failures_path,
        certificate_path,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run a frozen blind external-coupling holdout battery with hard "
            "taint gates and matched controls."
        )
    )
    parser.add_argument("--target-manifest", default=str(DEFAULT_TARGET_MANIFEST))
    parser.add_argument("--coupling-dir", default=str(DEFAULT_COUPLING_DIR))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--selector", default=DEFAULT_SELECTOR)
    parser.add_argument(
        "--target-id",
        action="append",
        default=[],
        help=(
            "Optional manifest target_id to run. Repeat for several. Omit to "
            "run every target."
        ),
    )
    parser.add_argument(
        "--control-name",
        action="append",
        choices=tuple(EXTERNAL_COUPLING_CONTROL_NAMES) + ("all",),
        default=[],
        help=(
            "Matched control to run. Repeat for several. Omit or pass 'all' "
            "to run every matched control."
        ),
    )
    args = parser.parse_args()
    for output in run_blind_external_holdout_battery_v0(
        target_manifest=Path(args.target_manifest),
        coupling_dir=Path(args.coupling_dir),
        output_dir=Path(args.output_dir),
        selector_name=args.selector,
        control_names=tuple(args.control_name),
        target_ids=tuple(args.target_id),
    ):
        print(output)


if __name__ == "__main__":
    main()
