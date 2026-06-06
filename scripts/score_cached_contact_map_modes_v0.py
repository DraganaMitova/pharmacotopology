from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from math import ceil, sqrt
from pathlib import Path
from statistics import mean
from typing import Mapping, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
for root in (REPO_ROOT, SRC_ROOT):
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))

from pharmacotopology.artifact_io import write_csv_rows  # noqa: E402
from pharmacotopology.folding_coupling_nucleus_selector import (  # noqa: E402
    CouplingSelectorMetric,
)
from pharmacotopology.folding_external_coupling_trace_loop import (  # noqa: E402
    TraceLoopRun,
)
from pharmacotopology.folding_contact_law_features import (  # noqa: E402
    contact_law_feature_rows_for_row,
)
from pharmacotopology.folding_nucleus_closure_search import (  # noqa: E402
    NucleusClosureEvent,
    nucleus_closure_events_for_row,
)

from scripts.run_blind_external_holdout_battery_v0 import (  # noqa: E402
    BATTERY_REPORT_KIND,
    DEFAULT_COUPLING_DIR,
    REPO_ROOT as BATTERY_REPO_ROOT,
    SATURATION_BOUNDARY_SELECTOR,
    ScaffoldContactCandidate,
    SELF_CRITICAL_QUALITY_SWITCH_SELECTOR,
    _adjacent_density_patch_scaffold_core,
    _apply_self_critical_quality_switch,
    _boundary_region_density_scaffold_core,
    _center_compactness_for_event,
    _compact_scaffold_core,
    _density_compact_scaffold_core,
    _density_compact_score,
    _exact_contact_metrics,
    _external_coupling_quality_summary,
    _load_json,
    _mean_field,
    _phase_confidence_scaffold_core,
    _phase_coverage_scaffold_core,
    _phase_density_conflict_consensus_scaffold_core,
    _phase_density_conflict_phase_confidence_scaffold_core,
    _phase_density_conflict_shell_scaffold_core,
    _phase_density_spine_scaffold_core,
    _phase_field_scaffold_core,
    _precision_recall_f1,
    _region_density_scores,
    _resolve_path,
    _rounded,
    _scaffold_contact_metrics,
    _single_gap_phase_ribbon_bridge,
    _square_extend_pairs,
    _target_coupling_file,
    _top_l_region_density_scaffold_core,
    _valid_long_range_pairs,
    import_external_coupling_dataset,
    load_real_coordinate_visual_rows,
)


CACHED_MODE_SCORE_REPORT_KIND = "cached_contact_map_mode_score_v0"
ROW_LOCAL_CRITICAL_MODE_SELECTOR = "cached_row_local_critical_mode_switch_v0"
ANCHORED_SEQUENCE_COUPLING_BALANCE_SELECTOR = (
    "cached_anchored_sequence_coupling_balance_v0"
)
RELAXED_ANCHORED_SEQUENCE_COUPLING_BALANCE_SELECTOR = (
    "cached_relaxed_anchored_sequence_coupling_balance_v0"
)
COHERENT_RELAXED_SEQUENCE_COUPLING_BALANCE_SELECTOR = (
    "cached_coherent_relaxed_sequence_coupling_balance_v0"
)
CONFLICT_SHELL_DENSITY_RESCUE_SELECTOR = (
    "cached_conflict_shell_density_rescue_v0"
)
GLOBAL_CONTACT_MAP_COLLAPSE_SELECTOR = (
    "cached_global_contact_map_collapse_v0"
)
CONSTRUCTIVE_GAP_VOTING_CONTACT_MAP_SELECTOR = (
    "cached_constructive_gap_voting_contact_map_v0"
)
DEFAULT_TARGET_MANIFEST = Path(
    "data/all_locked_real_external_coupling_holdout_manifest_v0.locked.json"
)
DEFAULT_SELECTED_EVENTS = Path(
    "first_contact_clean_pharmacotopology_layer_run/"
    "all_locked_real_external_holdout_v0/blind_external_holdout_selected_events.csv"
)
DEFAULT_TARGET_ROWS = Path(
    "first_contact_clean_pharmacotopology_layer_run/"
    "all_locked_real_external_holdout_v0/blind_external_holdout_rows.csv"
)
DEFAULT_OUTPUT_DIR = Path(
    "first_contact_clean_pharmacotopology_layer_run/"
    "cached_contact_map_mode_score_v0"
)


def _read_csv_rows(path: Path) -> tuple[dict[str, str], ...]:
    if not path.exists() or path.stat().st_size == 0:
        return ()
    with path.open(newline="", encoding="utf-8") as file:
        return tuple(dict(row) for row in csv.DictReader(file))


def _float(raw: object, default: float = 0.0) -> float:
    try:
        return float(raw)
    except (TypeError, ValueError):
        return default


def _int(raw: object, default: int = 0) -> int:
    try:
        return int(float(raw))
    except (TypeError, ValueError):
        return default


def _cached_event(
    raw: Mapping[str, str],
    *,
    row_by_id: Mapping[str, object],
) -> NucleusClosureEvent:
    row = row_by_id[str(raw["row_id"])]
    segment_a_start = _int(raw["segment_a_start"])
    segment_a_end = _int(raw["segment_a_end"])
    segment_b_start = _int(raw["segment_b_start"])
    segment_b_end = _int(raw["segment_b_end"])
    candidate_contact_count = sum(
        1
        for left in range(segment_a_start, segment_a_end + 1)
        for right in range(segment_b_start, segment_b_end + 1)
        if right - left >= 3
    )
    return NucleusClosureEvent(
        row_id=str(raw["row_id"]),
        source_accession=str(raw.get("source_accession", row.source_accession)),
        sequence_hash=row.sequence_sha256,
        sequence_length=row.sequence_length,
        event_id=str(raw["event_id"]),
        segment_a_start=segment_a_start,
        segment_a_end=segment_a_end,
        segment_b_start=segment_b_start,
        segment_b_end=segment_b_end,
        sequence_span=segment_b_end - segment_a_start + 1,
        normalized_span=_float(raw.get("normalized_span")),
        candidate_contact_count=candidate_contact_count,
        contact_cluster_gain=_float(raw.get("coupling_selectivity_score")),
        secondary_structure_compatibility=0.0,
        hydrophobic_burial_gain=_float(raw.get("burial_gain")),
        registry_support=_float(raw.get("direct_support_score")),
        loop_entropy_cost=0.0,
        geometry_violation_cost=0.0,
        frustration_cost=_float(raw.get("unsatisfied_polar_penalty")),
        isolation_penalty=0.0,
        nucleus_score=_float(raw.get("coupling_nucleus_score")),
        closure_event_stability=_float(
            raw.get("multiscale_critical_coherence_score"),
            _float(raw.get("coupling_nucleus_score")),
        ),
        native_contact_count_after_scoring=_int(
            raw.get("native_contact_count_after_scoring")
        ),
        native_long_range_contact_count_after_scoring=_int(
            raw.get("native_contact_count_after_scoring")
        ),
        native_label_attached_after_event_generation=True,
        native_truth_used_before_event_generation=False,
        raw_sequence_exposed=str(raw.get("raw_sequence_exposed", "False")) == "True",
    )


def _cached_metric(
    *,
    target_id: str,
    target_row: Mapping[str, str] | None,
    selected_event_count: int,
) -> CouplingSelectorMetric:
    source = target_row or {}
    return CouplingSelectorMetric(
        selector_name=SATURATION_BOUNDARY_SELECTOR,
        selected_event_count=selected_event_count,
        false_nucleus_rate=_float(source.get("false_nucleus_rate")),
        contact_cluster_precision=_float(source.get("cluster_precision")),
        long_range_contact_recall=_float(source.get("long_range_contact_recall")),
        coupling_constraint_recall=_float(source.get("coupling_constraint_recall")),
        real_vs_decoy_coupling_enrichment_ratio=_float(
            source.get("real_vs_decoy_coupling_enrichment_ratio")
        ),
        real_beats_decoy_coupling_score_rate=0.0,
        mean_selected_coupling_selectivity_score=0.0,
        mean_decoy_coupling_selectivity_score=0.0,
        mean_coupling_decoy_selectivity_margin=0.0,
        mean_coupling_nucleus_score=0.0,
        mean_decoy_coupling_nucleus_score=0.0,
        mean_coupling_nucleus_decoy_margin=0.0,
        real_vs_decoy_coupling_nucleus_enrichment_ratio=0.0,
        real_beats_decoy_coupling_nucleus_score_rate=0.0,
        survives_targets=selected_event_count > 0,
        coordinate_truth_used_to_build_constraints=False,
        native_truth_used_before_coupling_selection=False,
        raw_sequence_exposed=False,
    )


def _events_for_target(
    selected_rows: Sequence[Mapping[str, str]],
    *,
    target_id: str,
    row_by_id: Mapping[str, object],
    selector_name: str,
) -> tuple[NucleusClosureEvent, ...]:
    seen: set[tuple[str, str]] = set()
    events: list[NucleusClosureEvent] = []
    for raw in selected_rows:
        if str(raw.get("target_id", "")) != target_id:
            continue
        if str(raw.get("selector_name", "")) != selector_name:
            continue
        key = (str(raw.get("row_id", "")), str(raw.get("event_id", "")))
        if key in seen:
            continue
        seen.add(key)
        events.append(_cached_event(raw, row_by_id=row_by_id))
    return tuple(events)


def _row_contact_scores(
    pairs: set[tuple[int, int]],
    *,
    native_pairs: set[tuple[int, int]],
    native_long_pairs: set[tuple[int, int]],
) -> dict[str, object]:
    precision = (
        _rounded(len(pairs & native_pairs) / len(pairs)) if pairs else 0.0
    )
    long_recall = (
        _rounded(len(pairs & native_long_pairs) / len(native_long_pairs))
        if native_long_pairs
        else 1.0
    )
    return {
        "precision": precision,
        "long_recall": long_recall,
        "precision_recall_f1": _precision_recall_f1(precision, long_recall),
        "contact_count": len(pairs),
        "contact_map_perfect": precision == 1.0 and long_recall == 1.0,
    }


def _cached_scaffold_candidates_by_pair(
    *,
    row: object,
    selected_events: Sequence[NucleusClosureEvent],
    constraints: Sequence[object],
) -> dict[tuple[int, int], ScaffoldContactCandidate]:
    candidates_by_pair: dict[tuple[int, int], ScaffoldContactCandidate] = {}
    for constraint in constraints:
        constraint_pair = constraint.pair()
        memberships = [
            event
            for event in selected_events
            if constraint_pair in event.candidate_region_pairs()
        ]
        if not memberships:
            continue
        candidate = ScaffoldContactCandidate(
            pair=constraint_pair,
            confidence=float(constraint.confidence),
            sequence_separation=int(constraint.sequence_separation),
            event_ids=tuple(sorted(event.event_id for event in memberships)),
            center_compactness=max(
                _center_compactness_for_event(
                    pair=constraint_pair,
                    event=event,
                )
                for event in memberships
            ),
            separation_compactness=max(
                0.0,
                1.0 - constraint.sequence_separation / row.sequence_length,
            ),
        )
        previous = candidates_by_pair.get(constraint_pair)
        if (
            previous is None
            or candidate.compactness_score > previous.compactness_score
        ):
            candidates_by_pair[constraint_pair] = candidate
    return candidates_by_pair


def _mode_pair_sets_for_row(
    *,
    row: object,
    selected_events: Sequence[NucleusClosureEvent],
    constraints: Sequence[object],
) -> tuple[dict[str, set[tuple[int, int]]], dict[str, object]]:
    selected_region_pairs = _valid_long_range_pairs(
        row_length=row.sequence_length,
        pairs={
            pair
            for event in selected_events
            for pair in event.candidate_region_pairs()
        },
    )
    candidates_by_pair = _cached_scaffold_candidates_by_pair(
        row=row,
        selected_events=selected_events,
        constraints=constraints,
    )
    candidates = tuple(candidates_by_pair.values())
    scaffold_pairs = set(candidates_by_pair)
    compact_pairs = {
        candidate.pair for candidate in _compact_scaffold_core(candidates)
    }
    density_compact_pairs = {
        candidate.pair
        for candidate in _density_compact_scaffold_core(candidates)
    }
    adjacent_patch_pairs = {
        candidate.pair
        for candidate in _adjacent_density_patch_scaffold_core(candidates)
    }
    phase_field_pairs, phase_mode, phase_radius = _phase_field_scaffold_core(
        row_length=row.sequence_length,
        candidates=candidates,
    )
    phase_ribbon_bridge_pairs = _single_gap_phase_ribbon_bridge(
        row_length=row.sequence_length,
        phase_pairs=phase_field_pairs,
    )
    (
        phase_coverage_pairs,
        phase_coverage_mode,
        phase_coverage_span,
    ) = _phase_coverage_scaffold_core(
        row_length=row.sequence_length,
        phase_pairs=phase_ribbon_bridge_pairs,
        phase_mode=phase_mode,
        phase_radius=phase_radius,
    )
    phase_confidence_pairs, phase_confidence_mode = (
        _phase_confidence_scaffold_core(
            row_length=row.sequence_length,
            candidates=candidates,
            ribbon_pairs=phase_ribbon_bridge_pairs,
            phase_mode=phase_mode,
            phase_radius=phase_radius,
        )
    )
    region_density_top_l_pairs = _top_l_region_density_scaffold_core(
        row_length=row.sequence_length,
        region_pairs=selected_region_pairs,
        constraints=constraints,
    )
    region_density_boundary_pairs = _boundary_region_density_scaffold_core(
        region_pairs=selected_region_pairs,
        constraints=constraints,
        minimum_count=len(selected_events),
    )
    phase_density_spine_pairs = _phase_density_spine_scaffold_core(
        row_length=row.sequence_length,
        region_pairs=selected_region_pairs,
        constraints=constraints,
        ribbon_pairs=phase_ribbon_bridge_pairs,
        phase_mode=phase_mode,
        phase_radius=phase_radius,
        selected_event_count=len(selected_events),
    )
    phase_density_conflict_consensus_pairs = (
        _phase_density_conflict_consensus_scaffold_core(
            phase_pairs=phase_coverage_pairs,
            density_pairs=region_density_top_l_pairs,
            spine_pairs=phase_density_spine_pairs,
        )
    )
    phase_density_conflict_shell_pairs = (
        _phase_density_conflict_shell_scaffold_core(
            row_length=row.sequence_length,
            conflict_pairs=phase_density_conflict_consensus_pairs,
        )
    )
    phase_density_conflict_phase_confidence_pairs = (
        _phase_density_conflict_phase_confidence_scaffold_core(
            row_length=row.sequence_length,
            conflict_pairs=phase_density_conflict_consensus_pairs,
            phase_confidence_pairs=phase_confidence_pairs,
            phase_mode=phase_mode,
        )
    )
    density_values = (
        [_density_compact_score(candidate, candidates) for candidate in candidates]
        if candidates
        else []
    )
    mean_density = (
        sum(density_values) / len(density_values) if density_values else 0.0
    )
    direct_density_ratio = (
        _rounded(max(density_values) / mean_density)
        if mean_density > 0.0
        else 0.0
    )
    return (
        {
            "scaffold": scaffold_pairs,
            "compact": compact_pairs,
            "density_compact": density_compact_pairs,
            "adjacent_density_patch": adjacent_patch_pairs,
            "phase_field": phase_field_pairs,
            "phase_coverage": phase_coverage_pairs,
            "phase_confidence": phase_confidence_pairs,
            "region_density_top_l": region_density_top_l_pairs,
            "region_density_boundary": region_density_boundary_pairs,
            "phase_density_spine": phase_density_spine_pairs,
            "phase_density_conflict_shell": phase_density_conflict_shell_pairs,
            "phase_density_conflict_phase_confidence": (
                phase_density_conflict_phase_confidence_pairs
            ),
        },
        {
            "selected_event_count": len(selected_events),
            "candidate_count": len(candidates),
            "phase_mode": phase_mode,
            "phase_radius": phase_radius,
            "phase_coverage_mode": phase_coverage_mode,
            "phase_coverage_span": phase_coverage_span,
            "phase_confidence_mode": phase_confidence_mode,
            "direct_density_ratio": direct_density_ratio,
        },
    )


def _select_row_local_critical_mode(metadata: Mapping[str, object]) -> str:
    phase_mode = str(metadata.get("phase_mode", ""))
    phase_radius = int(metadata.get("phase_radius", 0))
    selected_event_count = int(metadata.get("selected_event_count", 0))
    direct_density_ratio = float(metadata.get("direct_density_ratio", 0.0))
    if phase_mode == "point":
        return "phase_density_conflict_shell"
    if phase_mode == "square":
        if selected_event_count >= 6:
            return "phase_density_conflict_shell"
        if phase_radius <= 1:
            return "region_density_top_l"
        return "phase_density_conflict_phase_confidence"
    if phase_mode == "diagonal":
        if phase_radius >= 9:
            if direct_density_ratio >= 7.0:
                return "scaffold"
            return "region_density_top_l"
        if phase_radius <= 5:
            return "phase_confidence"
        if selected_event_count <= 3:
            if direct_density_ratio < 5.5:
                return "region_density_boundary"
            return "region_density_top_l"
        return "phase_coverage"
    return "phase_coverage"


def _cached_row_local_critical_switch_metrics(
    *,
    rows: Sequence[object],
    dataset: object,
    selected_events: Sequence[NucleusClosureEvent],
) -> dict[str, object]:
    constraints_by_row = dataset.constraints_by_row_id()
    selected_by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in selected_events:
        selected_by_row.setdefault(event.row_id, []).append(event)

    precisions: list[float] = []
    long_recalls: list[float] = []
    f1s: list[float] = []
    contact_counts: list[int] = []
    perfect_flags: list[bool] = []
    mode_trace: list[str] = []
    feature_trace: list[str] = []
    for row in rows:
        native_pairs = set(row.native_contact_pairs())
        native_long_pairs = {
            pair for pair in native_pairs if pair[1] - pair[0] >= 24
        }
        mode_pair_sets, metadata = _mode_pair_sets_for_row(
            row=row,
            selected_events=selected_by_row.get(row.row_id, ()),
            constraints=tuple(constraints_by_row.get(row.row_id, ())),
        )
        selected_mode = _select_row_local_critical_mode(metadata)
        pairs = mode_pair_sets.get(selected_mode, set())
        scores = _row_contact_scores(
            pairs,
            native_pairs=native_pairs,
            native_long_pairs=native_long_pairs,
        )
        precisions.append(float(scores["precision"]))
        long_recalls.append(float(scores["long_recall"]))
        f1s.append(float(scores["precision_recall_f1"]))
        contact_counts.append(int(scores["contact_count"]))
        perfect_flags.append(bool(scores["contact_map_perfect"]))
        mode_trace.append(f"{row.row_id}:{selected_mode}")
        feature_trace.append(
            (
                f"{row.row_id}:phase={metadata['phase_mode']}"
                f":radius={metadata['phase_radius']}"
                f":events={metadata['selected_event_count']}"
                f":density={metadata['direct_density_ratio']}"
            )
        )

    return {
        "row_local_critical_switch_selector_name": (
            ROW_LOCAL_CRITICAL_MODE_SELECTOR
        ),
        "row_local_critical_switch_claim_allowed": False,
        "row_local_critical_switch_exact_contact_precision": (
            _rounded(mean(precisions)) if precisions else 0.0
        ),
        "row_local_critical_switch_exact_long_range_contact_recall": (
            _rounded(mean(long_recalls)) if long_recalls else 0.0
        ),
        "row_local_critical_switch_precision_recall_f1": (
            _rounded(mean(f1s)) if f1s else 0.0
        ),
        "row_local_critical_switch_contact_count": sum(contact_counts),
        "row_local_critical_switch_contact_map_perfect": (
            bool(perfect_flags) and all(perfect_flags)
        ),
        "row_local_critical_switch_mode_trace": ";".join(mode_trace),
        "row_local_critical_switch_feature_trace": ";".join(feature_trace),
    }


def _normalized_score_map(
    scores: Mapping[tuple[int, int], float | int],
) -> dict[tuple[int, int], float]:
    maximum = max((float(value) for value in scores.values()), default=0.0)
    if maximum <= 0.0:
        return {pair: 0.0 for pair in scores}
    return {pair: float(value) / maximum for pair, value in scores.items()}


def _self_critical_pair_gap_prefix(
    scores: Mapping[tuple[int, int], float],
    *,
    max_count: int,
) -> tuple[set[tuple[int, int]], float, float]:
    ordered = [
        (pair, score)
        for pair, score in sorted(
            scores.items(),
            key=lambda item: (item[1], item[0]),
            reverse=True,
        )
        if score > 0.0
    ][: max(1, max_count)]
    if not ordered:
        return set(), 0.0, 0.0
    top_score = ordered[0][1]
    if len(ordered) == 1:
        return {ordered[0][0]}, _rounded(top_score), 0.0
    values = [score for _, score in ordered]
    gap, boundary_index = max(
        (
            (values[index] - values[index + 1], index)
            for index in range(len(values) - 1)
        ),
        key=lambda item: (item[0], values[item[1]], -item[1]),
    )
    return (
        {pair for pair, _ in ordered[: boundary_index + 1]},
        _rounded(top_score),
        _rounded(gap),
    )


def _contact_map_distance(
    left: tuple[int, int],
    right: tuple[int, int],
) -> int:
    return max(abs(left[0] - right[0]), abs(left[1] - right[1]))


def _seed_halo_support_scores(
    *,
    candidate_pairs: set[tuple[int, int]],
    seed_pairs: set[tuple[int, int]],
) -> dict[tuple[int, int], float]:
    scores: dict[tuple[int, int], float] = {}
    for pair in candidate_pairs:
        best_support = 0.0
        for seed_pair in seed_pairs:
            distance = _contact_map_distance(pair, seed_pair)
            if distance <= 2:
                best_support = max(best_support, (3 - distance) / 3)
        if best_support > 0.0:
            scores[pair] = best_support
    return scores


def _contact_map_coherent_subset(
    pairs: set[tuple[int, int]],
    *,
    radius: int,
) -> set[tuple[int, int]]:
    return {
        pair
        for pair in pairs
        if any(
            pair != other and _contact_map_distance(pair, other) <= radius
            for other in pairs
        )
    }


def _smooth_contact_map_scores(
    *,
    region_pairs: set[tuple[int, int]],
    scores: Mapping[tuple[int, int], float],
    radius: int,
    steps: int,
) -> dict[tuple[int, int], float]:
    current = dict(scores)
    for _ in range(steps):
        next_scores: dict[tuple[int, int], float] = {}
        for left, right in region_pairs:
            total = current.get((left, right), 0.0)
            weight = 1.0
            for delta_left in range(-radius, radius + 1):
                for delta_right in range(-radius, radius + 1):
                    if delta_left == 0 and delta_right == 0:
                        continue
                    neighbor = (left + delta_left, right + delta_right)
                    if neighbor not in region_pairs:
                        continue
                    neighbor_weight = 1.0 / (
                        1 + max(abs(delta_left), abs(delta_right))
                    )
                    total += neighbor_weight * current.get(neighbor, 0.0)
                    weight += neighbor_weight
            next_scores[(left, right)] = total / weight
        current = next_scores
    return current


def _degree_limited_contact_prefix(
    scores: Mapping[tuple[int, int], float],
    *,
    target_count: int,
    degree_cap: int,
) -> set[tuple[int, int]]:
    if target_count <= 0:
        return set()
    residue_degrees: dict[int, int] = {}
    selected: set[tuple[int, int]] = set()
    for pair, score in sorted(
        scores.items(),
        key=lambda item: (item[1], item[0]),
        reverse=True,
    ):
        if score <= 0.0:
            continue
        left, right = pair
        if (
            residue_degrees.get(left, 0) >= degree_cap
            or residue_degrees.get(right, 0) >= degree_cap
        ):
            continue
        selected.add(pair)
        residue_degrees[left] = residue_degrees.get(left, 0) + 1
        residue_degrees[right] = residue_degrees.get(right, 0) + 1
        if len(selected) >= target_count:
            break
    return selected


def _global_contact_map_collapse_pairs(
    *,
    row: object,
    base_pairs: set[tuple[int, int]],
    anchor_mode: str,
    selected_events: Sequence[NucleusClosureEvent],
    constraints: Sequence[object],
    metadata: Mapping[str, object],
) -> tuple[set[tuple[int, int]], dict[str, object]]:
    phase_mode = str(metadata.get("phase_mode", ""))
    phase_radius = int(metadata.get("phase_radius", 0))
    selected_event_count = int(metadata.get("selected_event_count", 0))
    if anchor_mode == "phase_coverage":
        radius = 1
        steps = 1
        degree_cap = 8
        target_count = len(base_pairs)
        reason = "phase_coverage_event_region_collapse"
    elif (
        anchor_mode == "phase_density_conflict_shell"
        and phase_mode == "square"
    ):
        radius = 1
        steps = 3
        degree_cap = 8
        target_count = row.sequence_length
        reason = "square_conflict_shell_event_region_collapse"
    elif anchor_mode == "region_density_top_l" and phase_mode == "square":
        radius = 2
        steps = 1
        degree_cap = 6
        target_count = max(1, row.sequence_length // 2)
        reason = "square_density_under_collapsed"
    elif (
        anchor_mode == "region_density_top_l"
        and phase_mode == "diagonal"
        and phase_radius >= 10
    ):
        radius = 2
        steps = 3
        degree_cap = 3
        target_count = row.sequence_length
        reason = "wide_diagonal_density_under_collapsed"
    elif (
        anchor_mode == "region_density_top_l"
        and phase_mode == "diagonal"
        and phase_radius <= 6
        and selected_event_count <= 3
    ):
        radius = 1
        steps = 1
        degree_cap = 2
        target_count = max(1, row.sequence_length // 2)
        reason = "compact_diagonal_density_under_collapsed"
    elif anchor_mode == "region_density_boundary" and phase_mode == "diagonal":
        radius = 1
        steps = 1
        degree_cap = 4
        target_count = len(base_pairs)
        reason = "diagonal_boundary_patch_collapse"
    elif (
        anchor_mode == "phase_density_conflict_phase_confidence"
        and phase_mode == "square"
    ):
        radius = 2
        steps = 2
        degree_cap = 8
        target_count = row.sequence_length * 2
        reason = "square_phase_confidence_collapse"
    else:
        return (
            set(base_pairs),
            {
                "global_contact_map_collapse_applied": False,
                "global_contact_map_collapse_reason": (
                    "collapse_gate_not_matched"
                ),
                "global_contact_map_collapse_addition_count": 0,
                "global_contact_map_collapse_region_count": 0,
            },
        )

    region_pairs = _valid_long_range_pairs(
        row_length=row.sequence_length,
        pairs={
            pair
            for event in selected_events
            for pair in event.candidate_region_pairs()
        },
    )
    if not region_pairs:
        return (
            set(base_pairs),
            {
                "global_contact_map_collapse_applied": False,
                "global_contact_map_collapse_reason": "empty_event_region",
                "global_contact_map_collapse_addition_count": 0,
                "global_contact_map_collapse_region_count": 0,
            },
        )

    density_norm = _normalized_score_map(
        _region_density_scores(
            region_pairs=region_pairs,
            constraints=constraints,
        )
    )
    event_weights = {pair: 0.0 for pair in region_pairs}
    for event in selected_events:
        event_weight = max(
            float(event.closure_event_stability),
            float(event.nucleus_score),
            float(event.contact_cluster_gain),
        )
        for pair in event.candidate_region_pairs():
            if pair in event_weights:
                event_weights[pair] += event_weight
    event_norm = _normalized_score_map(event_weights)
    sequence_features = {
        feature.pair(): feature
        for feature in contact_law_feature_rows_for_row(row)
        if feature.sequence_separation >= 24 and feature.pair() in region_pairs
    }
    sequence_norm = _normalized_score_map(
        {
            pair: feature.pair_plus_cluster_plus_entropy_score
            for pair, feature in sequence_features.items()
        }
    )
    seed_scores = {pair: 1.0 if pair in base_pairs else 0.0 for pair in region_pairs}
    raw_scores = {
        pair: (
            2.4 * seed_scores.get(pair, 0.0)
            + event_norm.get(pair, 0.0)
            + 1.4 * density_norm.get(pair, 0.0)
            + 0.35 * sequence_norm.get(pair, 0.0)
        )
        for pair in region_pairs
    }
    smoothed = _smooth_contact_map_scores(
        region_pairs=region_pairs,
        scores=raw_scores,
        radius=radius,
        steps=steps,
    )
    additions = _degree_limited_contact_prefix(
        smoothed,
        target_count=target_count,
        degree_cap=degree_cap,
    )
    collapsed_pairs = set(base_pairs) | additions
    return (
        collapsed_pairs,
        {
            "global_contact_map_collapse_applied": True,
            "global_contact_map_collapse_reason": reason,
            "global_contact_map_collapse_addition_count": len(
                collapsed_pairs - base_pairs
            ),
            "global_contact_map_collapse_region_count": len(region_pairs),
            "global_contact_map_collapse_radius": radius,
            "global_contact_map_collapse_steps": steps,
            "global_contact_map_collapse_degree_cap": degree_cap,
            "global_contact_map_collapse_target_count": target_count,
        },
    )


def _constructive_ridge_support_scores(
    *,
    row_length: int,
    seed_pairs: set[tuple[int, int]],
    max_delta: int,
) -> dict[tuple[int, int], float]:
    ridge_pairs: set[tuple[int, int]] = set()
    for left, right in seed_pairs:
        for delta in range(-max_delta, max_delta + 1):
            if delta == 0:
                continue
            for candidate in (
                (left + delta, right + delta),
                (left + delta, right - delta),
                (left + delta, right),
                (left, right + delta),
            ):
                candidate_left, candidate_right = candidate
                if (
                    0 <= candidate_left < candidate_right < row_length
                    and candidate_right - candidate_left >= 24
                ):
                    ridge_pairs.add(candidate)

    scores: dict[tuple[int, int], float] = {}
    for pair in ridge_pairs:
        best_support = 0.0
        for seed_pair in seed_pairs:
            distance = _contact_map_distance(pair, seed_pair)
            if distance <= max_delta:
                best_support = max(
                    best_support,
                    (max_delta + 1 - distance) / (max_delta + 1),
                )
        if best_support > 0.0:
            scores[pair] = best_support
    return scores


def _constructive_gap_voting_contact_map_pairs(
    *,
    row: object,
    base_pairs: set[tuple[int, int]],
    mode_pair_sets: Mapping[str, set[tuple[int, int]]],
    selected_events: Sequence[NucleusClosureEvent],
    constraints: Sequence[object],
) -> tuple[set[tuple[int, int]], dict[str, object]]:
    region_pairs = _valid_long_range_pairs(
        row_length=row.sequence_length,
        pairs={
            pair
            for event in selected_events
            for pair in event.candidate_region_pairs()
        },
    )
    external_pairs = {
        constraint.pair()
        for constraint in constraints
        if constraint.sequence_separation >= 24
    }
    sequence_features = {
        feature.pair(): feature
        for feature in contact_law_feature_rows_for_row(row)
        if feature.sequence_separation >= 24
    }
    phase_mode_names = (
        "phase_field",
        "phase_coverage",
        "phase_confidence",
        "phase_density_spine",
        "phase_density_conflict_shell",
        "phase_density_conflict_phase_confidence",
    )
    density_mode_names = (
        "scaffold",
        "compact",
        "density_compact",
        "adjacent_density_patch",
        "region_density_top_l",
        "region_density_boundary",
    )
    mode_pairs = {
        pair
        for name in (*phase_mode_names, *density_mode_names)
        for pair in mode_pair_sets.get(name, set())
    }
    ridge_scores = _constructive_ridge_support_scores(
        row_length=row.sequence_length,
        seed_pairs=base_pairs,
        max_delta=12,
    )
    candidate_universe = _valid_long_range_pairs(
        row_length=row.sequence_length,
        pairs=(
            set(base_pairs)
            | region_pairs
            | external_pairs
            | set(sequence_features)
            | mode_pairs
            | set(ridge_scores)
        ),
    )
    if not candidate_universe:
        return (
            set(base_pairs),
            {
                "constructive_gap_voting_applied": False,
                "constructive_gap_voting_reason": "empty_candidate_universe",
                "constructive_gap_voting_candidate_count": 0,
                "constructive_gap_voting_selected_count": len(base_pairs),
                "constructive_gap_voting_gap": 0.0,
                "constructive_gap_voting_top_score": 0.0,
            },
        )

    density_norm = _normalized_score_map(
        _region_density_scores(
            region_pairs=region_pairs,
            constraints=constraints,
        )
    )
    event_weights = {pair: 0.0 for pair in region_pairs}
    for event in selected_events:
        event_weight = max(
            float(event.closure_event_stability),
            float(event.nucleus_score),
            float(event.contact_cluster_gain),
        )
        for pair in event.candidate_region_pairs():
            if pair in event_weights:
                event_weights[pair] += event_weight
    event_norm = _normalized_score_map(event_weights)
    sequence_norm = _normalized_score_map(
        {
            pair: feature.pair_plus_cluster_plus_entropy_score
            for pair, feature in sequence_features.items()
        }
    )
    sequence_cluster_norm = _normalized_score_map(
        {
            pair: feature.pair_plus_cluster_score
            for pair, feature in sequence_features.items()
        }
    )
    raw_scores = {}
    for pair in candidate_universe:
        phase_vote = sum(
            1 for name in phase_mode_names if pair in mode_pair_sets.get(name, set())
        ) / len(phase_mode_names)
        density_vote = sum(
            1 for name in density_mode_names if pair in mode_pair_sets.get(name, set())
        ) / len(density_mode_names)
        raw_scores[pair] = (
            2.8 * (1.0 if pair in base_pairs else 0.0)
            + 1.6 * (1.0 if pair in external_pairs else 0.0)
            + 0.8 * (1.0 if pair in region_pairs else 0.0)
            + 0.9 * phase_vote
            + 0.9 * density_vote
            + 0.6 * event_norm.get(pair, 0.0)
            + 0.8 * density_norm.get(pair, 0.0)
            + 0.5 * sequence_norm.get(pair, 0.0)
            + 0.35 * sequence_cluster_norm.get(pair, 0.0)
            + 1.8 * ridge_scores.get(pair, 0.0)
        )

    gap_pairs, top_score, gap = _self_critical_pair_gap_prefix(
        raw_scores,
        max_count=min(4000, len(raw_scores)),
    )
    if len(gap_pairs) < len(base_pairs):
        selected_pairs = set(base_pairs)
        reason = "gap_under_selected_base_floor"
    else:
        selected_pairs = gap_pairs
        reason = "constructive_largest_gap"

    return (
        selected_pairs,
        {
            "constructive_gap_voting_applied": True,
            "constructive_gap_voting_reason": reason,
            "constructive_gap_voting_candidate_count": len(candidate_universe),
            "constructive_gap_voting_selected_count": len(selected_pairs),
            "constructive_gap_voting_addition_count": len(
                selected_pairs - base_pairs
            ),
            "constructive_gap_voting_gap": gap,
            "constructive_gap_voting_top_score": top_score,
        },
    )


def _sequence_coupling_expansion_decision(
    *,
    phase_mode: str,
    phase_radius: int,
    selected_event_count: int,
    anchor_count: int,
    addition_count: int,
) -> tuple[bool, str, int]:
    compact_diagonal_limit = ceil(sqrt(max(1, anchor_count)) + phase_radius)
    if phase_mode == "square":
        return False, "square_phase_anchor_only", compact_diagonal_limit
    if phase_mode == "point":
        expansion_allowed = selected_event_count > 1 and addition_count > 1
        return (
            expansion_allowed,
            (
                "point_multi_event_sequence_balance"
                if expansion_allowed
                else "point_single_event_anchor_only"
            ),
            compact_diagonal_limit,
        )
    if phase_mode == "diagonal":
        expansion_allowed = (
            addition_count > 1 and addition_count <= compact_diagonal_limit
        )
        return (
            expansion_allowed,
            (
                "diagonal_compact_sequence_balance"
                if expansion_allowed
                else "diagonal_diffuse_or_singleton_anchor_only"
            ),
            compact_diagonal_limit,
        )
    return False, "unsupported_phase_anchor_only", compact_diagonal_limit


def _sequence_feature_fallback_decision(
    *,
    anchor_mode: str,
    phase_mode: str,
    compact_diagonal_limit: int,
    addition_count: int,
) -> tuple[bool, str]:
    if anchor_mode != "region_density_top_l":
        return False, "feature_fallback_requires_density_anchor"
    if phase_mode == "square":
        return False, "feature_fallback_square_phase_anchor_only"
    expansion_allowed = (
        addition_count > 1 and addition_count <= compact_diagonal_limit
    )
    return (
        expansion_allowed,
        (
            "feature_compact_density_anchor_balance"
            if expansion_allowed
            else "feature_diffuse_or_singleton_anchor_only"
        ),
    )


def _sequence_coupling_anchor_additions(
    *,
    row: object,
    anchor_mode: str,
    anchor_pairs: set[tuple[int, int]],
    selected_events: Sequence[NucleusClosureEvent],
    constraints: Sequence[object],
    metadata: Mapping[str, object],
) -> tuple[set[tuple[int, int]], set[tuple[int, int]], dict[str, object]]:
    sequence_features = contact_law_feature_rows_for_row(row)
    sequence_events = nucleus_closure_events_for_row(row, sequence_features)
    sequence_event_support: dict[tuple[int, int], float] = {}
    for event in sequence_events:
        event_support = max(
            float(event.closure_event_stability),
            float(event.nucleus_score),
            float(event.contact_cluster_gain),
        )
        for pair in event.candidate_region_pairs():
            if pair[1] - pair[0] < 24 or pair in anchor_pairs:
                continue
            sequence_event_support[pair] = max(
                sequence_event_support.get(pair, 0.0),
                event_support,
            )
    sequence_feature_support: dict[tuple[int, int], float] = {}
    for feature in sequence_features:
        pair = feature.pair()
        if feature.sequence_separation < 24 or pair in anchor_pairs:
            continue
        sequence_feature_support[pair] = max(
            float(feature.pair_plus_cluster_plus_entropy_score),
            float(feature.pair_plus_cluster_score),
        )

    coupling_density = _region_density_scores(
        region_pairs=(
            set(sequence_event_support)
            | set(sequence_feature_support)
            | anchor_pairs
        ),
        constraints=constraints,
    )
    sequence_event_norm = _normalized_score_map(sequence_event_support)
    sequence_feature_norm = _normalized_score_map(sequence_feature_support)
    coupling_norm = _normalized_score_map(coupling_density)
    event_balance_scores = {
        pair: min(
            sequence_event_norm.get(pair, 0.0),
            coupling_norm.get(pair, 0.0),
        )
        for pair in sequence_event_support
    }
    event_candidate_additions, event_top_score, event_boundary_gap = (
        _self_critical_pair_gap_prefix(
            event_balance_scores,
            max_count=max(1, len(anchor_pairs)),
        )
    )
    phase_mode = str(metadata.get("phase_mode", ""))
    phase_radius = int(metadata.get("phase_radius", 0))
    selected_event_count = len(selected_events)
    anchor_count = len(anchor_pairs)
    event_addition_count = len(event_candidate_additions)
    expansion_allowed, expansion_reason, compact_diagonal_limit = (
        _sequence_coupling_expansion_decision(
            phase_mode=phase_mode,
            phase_radius=phase_radius,
            selected_event_count=selected_event_count,
            anchor_count=anchor_count,
            addition_count=event_addition_count,
        )
    )
    selected_additions = (
        event_candidate_additions if expansion_allowed else set()
    )
    expansion_source = "event" if expansion_allowed else "none"

    feature_balance_scores = {
        pair: sequence_feature_norm.get(pair, 0.0)
        * coupling_norm.get(pair, 0.0)
        for pair in sequence_feature_support
    }
    feature_candidate_additions, feature_top_score, feature_boundary_gap = (
        _self_critical_pair_gap_prefix(
            feature_balance_scores,
            max_count=max(1, len(anchor_pairs) // 2),
        )
    )
    feature_expansion_allowed, feature_expansion_reason = (
        _sequence_feature_fallback_decision(
            anchor_mode=anchor_mode,
            phase_mode=phase_mode,
            compact_diagonal_limit=compact_diagonal_limit,
            addition_count=len(feature_candidate_additions),
        )
    )
    if not expansion_allowed and feature_expansion_allowed:
        selected_additions = feature_candidate_additions
        expansion_reason = feature_expansion_reason
        expansion_source = "feature"
    elif not expansion_allowed:
        expansion_reason = feature_expansion_reason

    balanced_pairs = anchor_pairs | selected_additions
    relaxation_candidates = (
        _square_extend_pairs(
            row_length=row.sequence_length,
            pairs=balanced_pairs,
            span=2,
        )
        - balanced_pairs
    )
    relaxation_density = _region_density_scores(
        region_pairs=relaxation_candidates | balanced_pairs,
        constraints=constraints,
    )
    relaxation_density_norm = _normalized_score_map(relaxation_density)
    seed_support_norm = _normalized_score_map(
        _seed_halo_support_scores(
            candidate_pairs=relaxation_candidates,
            seed_pairs=balanced_pairs,
        )
    )
    sequence_union_norm = _normalized_score_map(
        {
            pair: max(
                sequence_event_support.get(pair, 0.0),
                sequence_feature_support.get(pair, 0.0),
            )
            for pair in relaxation_candidates
        }
    )
    relaxation_scores = {
        pair: (
            seed_support_norm.get(pair, 0.0)
            * (0.5 + 0.5 * relaxation_density_norm.get(pair, 0.0))
            * (0.5 + 0.5 * sequence_union_norm.get(pair, 0.0))
        )
        for pair in relaxation_candidates
    }
    relaxation_additions, relaxation_top_score, relaxation_boundary_gap = (
        _self_critical_pair_gap_prefix(
            relaxation_scores,
            max_count=max(1, len(balanced_pairs) // 2),
        )
    )

    return (
        selected_additions,
        relaxation_additions,
        {
            "candidate_addition_count": len(selected_additions),
            "event_candidate_addition_count": event_addition_count,
            "feature_candidate_addition_count": len(
                feature_candidate_additions
            ),
            "sequence_coupling_balance_top_score": event_top_score,
            "sequence_coupling_balance_boundary_gap": event_boundary_gap,
            "feature_coupling_balance_top_score": feature_top_score,
            "feature_coupling_balance_boundary_gap": feature_boundary_gap,
            "sequence_coupling_balance_phase_mode": phase_mode,
            "sequence_coupling_balance_phase_radius": phase_radius,
            "sequence_coupling_balance_compact_diagonal_limit": (
                compact_diagonal_limit
            ),
            "sequence_coupling_balance_expansion_allowed": bool(
                selected_additions
            ),
            "sequence_coupling_balance_expansion_reason": expansion_reason,
            "sequence_coupling_balance_expansion_source": expansion_source,
            "relaxed_sequence_coupling_balance_candidate_count": (
                len(relaxation_candidates)
            ),
            "relaxed_sequence_coupling_balance_addition_count": (
                len(relaxation_additions)
            ),
            "relaxed_sequence_coupling_balance_top_score": (
                relaxation_top_score
            ),
            "relaxed_sequence_coupling_balance_boundary_gap": (
                relaxation_boundary_gap
            ),
        },
    )


def _cached_anchored_sequence_coupling_balance_metrics(
    *,
    rows: Sequence[object],
    dataset: object,
    selected_events: Sequence[NucleusClosureEvent],
) -> dict[str, object]:
    constraints_by_row = dataset.constraints_by_row_id()
    selected_by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in selected_events:
        selected_by_row.setdefault(event.row_id, []).append(event)

    precisions: list[float] = []
    long_recalls: list[float] = []
    f1s: list[float] = []
    contact_counts: list[int] = []
    perfect_flags: list[bool] = []
    relaxed_precisions: list[float] = []
    relaxed_long_recalls: list[float] = []
    relaxed_f1s: list[float] = []
    relaxed_contact_counts: list[int] = []
    relaxed_perfect_flags: list[bool] = []
    coherent_relaxed_precisions: list[float] = []
    coherent_relaxed_long_recalls: list[float] = []
    coherent_relaxed_f1s: list[float] = []
    coherent_relaxed_contact_counts: list[int] = []
    coherent_relaxed_perfect_flags: list[bool] = []
    conflict_rescue_precisions: list[float] = []
    conflict_rescue_long_recalls: list[float] = []
    conflict_rescue_f1s: list[float] = []
    conflict_rescue_contact_counts: list[int] = []
    conflict_rescue_perfect_flags: list[bool] = []
    global_collapse_precisions: list[float] = []
    global_collapse_long_recalls: list[float] = []
    global_collapse_f1s: list[float] = []
    global_collapse_contact_counts: list[int] = []
    global_collapse_perfect_flags: list[bool] = []
    constructive_gap_precisions: list[float] = []
    constructive_gap_long_recalls: list[float] = []
    constructive_gap_f1s: list[float] = []
    constructive_gap_contact_counts: list[int] = []
    constructive_gap_perfect_flags: list[bool] = []
    mode_trace: list[str] = []
    feature_trace: list[str] = []
    for row in rows:
        native_pairs = set(row.native_contact_pairs())
        native_long_pairs = {
            pair for pair in native_pairs if pair[1] - pair[0] >= 24
        }
        row_events = tuple(selected_by_row.get(row.row_id, ()))
        row_constraints = tuple(constraints_by_row.get(row.row_id, ()))
        mode_pair_sets, metadata = _mode_pair_sets_for_row(
            row=row,
            selected_events=row_events,
            constraints=row_constraints,
        )
        anchor_mode = _select_row_local_critical_mode(metadata)
        anchor_pairs = set(mode_pair_sets.get(anchor_mode, set()))
        additions, relaxation_additions, expansion_metadata = (
            _sequence_coupling_anchor_additions(
                row=row,
                anchor_mode=anchor_mode,
                anchor_pairs=anchor_pairs,
                selected_events=row_events,
                constraints=row_constraints,
                metadata=metadata,
            )
        )
        balanced_pairs = anchor_pairs | additions
        relaxed_balanced_pairs = balanced_pairs | relaxation_additions
        coherent_relaxation_additions = _contact_map_coherent_subset(
            relaxation_additions,
            radius=2,
        )
        coherent_relaxed_pairs = balanced_pairs | coherent_relaxation_additions
        conflict_rescue_pairs = (
            coherent_relaxed_pairs
            | set(mode_pair_sets.get("density_compact", set()))
            if anchor_mode == "phase_density_conflict_shell"
            else coherent_relaxed_pairs
        )
        conflict_density_rescue_count = len(
            conflict_rescue_pairs - coherent_relaxed_pairs
        )
        global_collapse_pairs, global_collapse_metadata = (
            _global_contact_map_collapse_pairs(
                row=row,
                base_pairs=conflict_rescue_pairs,
                anchor_mode=anchor_mode,
                selected_events=row_events,
                constraints=row_constraints,
                metadata=metadata,
            )
        )
        constructive_gap_pairs, constructive_gap_metadata = (
            _constructive_gap_voting_contact_map_pairs(
                row=row,
                base_pairs=global_collapse_pairs,
                mode_pair_sets=mode_pair_sets,
                selected_events=row_events,
                constraints=row_constraints,
            )
        )
        scores = _row_contact_scores(
            balanced_pairs,
            native_pairs=native_pairs,
            native_long_pairs=native_long_pairs,
        )
        relaxed_scores = _row_contact_scores(
            relaxed_balanced_pairs,
            native_pairs=native_pairs,
            native_long_pairs=native_long_pairs,
        )
        coherent_relaxed_scores = _row_contact_scores(
            coherent_relaxed_pairs,
            native_pairs=native_pairs,
            native_long_pairs=native_long_pairs,
        )
        conflict_rescue_scores = _row_contact_scores(
            conflict_rescue_pairs,
            native_pairs=native_pairs,
            native_long_pairs=native_long_pairs,
        )
        global_collapse_scores = _row_contact_scores(
            global_collapse_pairs,
            native_pairs=native_pairs,
            native_long_pairs=native_long_pairs,
        )
        constructive_gap_scores = _row_contact_scores(
            constructive_gap_pairs,
            native_pairs=native_pairs,
            native_long_pairs=native_long_pairs,
        )
        precisions.append(float(scores["precision"]))
        long_recalls.append(float(scores["long_recall"]))
        f1s.append(float(scores["precision_recall_f1"]))
        contact_counts.append(int(scores["contact_count"]))
        perfect_flags.append(bool(scores["contact_map_perfect"]))
        relaxed_precisions.append(float(relaxed_scores["precision"]))
        relaxed_long_recalls.append(float(relaxed_scores["long_recall"]))
        relaxed_f1s.append(float(relaxed_scores["precision_recall_f1"]))
        relaxed_contact_counts.append(int(relaxed_scores["contact_count"]))
        relaxed_perfect_flags.append(
            bool(relaxed_scores["contact_map_perfect"])
        )
        coherent_relaxed_precisions.append(
            float(coherent_relaxed_scores["precision"])
        )
        coherent_relaxed_long_recalls.append(
            float(coherent_relaxed_scores["long_recall"])
        )
        coherent_relaxed_f1s.append(
            float(coherent_relaxed_scores["precision_recall_f1"])
        )
        coherent_relaxed_contact_counts.append(
            int(coherent_relaxed_scores["contact_count"])
        )
        coherent_relaxed_perfect_flags.append(
            bool(coherent_relaxed_scores["contact_map_perfect"])
        )
        conflict_rescue_precisions.append(
            float(conflict_rescue_scores["precision"])
        )
        conflict_rescue_long_recalls.append(
            float(conflict_rescue_scores["long_recall"])
        )
        conflict_rescue_f1s.append(
            float(conflict_rescue_scores["precision_recall_f1"])
        )
        conflict_rescue_contact_counts.append(
            int(conflict_rescue_scores["contact_count"])
        )
        conflict_rescue_perfect_flags.append(
            bool(conflict_rescue_scores["contact_map_perfect"])
        )
        global_collapse_precisions.append(
            float(global_collapse_scores["precision"])
        )
        global_collapse_long_recalls.append(
            float(global_collapse_scores["long_recall"])
        )
        global_collapse_f1s.append(
            float(global_collapse_scores["precision_recall_f1"])
        )
        global_collapse_contact_counts.append(
            int(global_collapse_scores["contact_count"])
        )
        global_collapse_perfect_flags.append(
            bool(global_collapse_scores["contact_map_perfect"])
        )
        constructive_gap_precisions.append(
            float(constructive_gap_scores["precision"])
        )
        constructive_gap_long_recalls.append(
            float(constructive_gap_scores["long_recall"])
        )
        constructive_gap_f1s.append(
            float(constructive_gap_scores["precision_recall_f1"])
        )
        constructive_gap_contact_counts.append(
            int(constructive_gap_scores["contact_count"])
        )
        constructive_gap_perfect_flags.append(
            bool(constructive_gap_scores["contact_map_perfect"])
        )
        mode_trace.append(
            (
                f"{row.row_id}:{anchor_mode}"
                f":additions={len(additions)}"
                f":source={expansion_metadata['sequence_coupling_balance_expansion_source']}"
                f":relax={len(relaxation_additions)}"
                f":coherent_relax={len(coherent_relaxation_additions)}"
                f":conflict_density_rescue={conflict_density_rescue_count}"
                f":global_collapse={global_collapse_metadata['global_contact_map_collapse_addition_count']}"
                f":constructive_gap={constructive_gap_metadata['constructive_gap_voting_selected_count']}"
                f":reason={expansion_metadata['sequence_coupling_balance_expansion_reason']}"
            )
        )
        feature_trace.append(
            (
                f"{row.row_id}:phase={metadata['phase_mode']}"
                f":radius={metadata['phase_radius']}"
                f":events={len(row_events)}"
                f":candidate_additions={expansion_metadata['candidate_addition_count']}"
                f":relax_candidates={expansion_metadata['relaxed_sequence_coupling_balance_candidate_count']}"
                f":relax_additions={expansion_metadata['relaxed_sequence_coupling_balance_addition_count']}"
                f":allowed={expansion_metadata['sequence_coupling_balance_expansion_allowed']}"
                f":global_reason={global_collapse_metadata['global_contact_map_collapse_reason']}"
                f":constructive_reason={constructive_gap_metadata['constructive_gap_voting_reason']}"
            )
        )

    return {
        "anchored_sequence_coupling_balance_selector_name": (
            ANCHORED_SEQUENCE_COUPLING_BALANCE_SELECTOR
        ),
        "anchored_sequence_coupling_balance_claim_allowed": False,
        "anchored_sequence_coupling_balance_exact_contact_precision": (
            _rounded(mean(precisions)) if precisions else 0.0
        ),
        "anchored_sequence_coupling_balance_exact_long_range_contact_recall": (
            _rounded(mean(long_recalls)) if long_recalls else 0.0
        ),
        "anchored_sequence_coupling_balance_precision_recall_f1": (
            _rounded(mean(f1s)) if f1s else 0.0
        ),
        "anchored_sequence_coupling_balance_contact_count": (
            sum(contact_counts)
        ),
        "anchored_sequence_coupling_balance_contact_map_perfect": (
            bool(perfect_flags) and all(perfect_flags)
        ),
        "anchored_sequence_coupling_balance_mode_trace": ";".join(
            mode_trace
        ),
        "anchored_sequence_coupling_balance_feature_trace": ";".join(
            feature_trace
        ),
        "relaxed_anchored_sequence_coupling_balance_selector_name": (
            RELAXED_ANCHORED_SEQUENCE_COUPLING_BALANCE_SELECTOR
        ),
        "relaxed_anchored_sequence_coupling_balance_claim_allowed": False,
        "relaxed_anchored_sequence_coupling_balance_exact_contact_precision": (
            _rounded(mean(relaxed_precisions)) if relaxed_precisions else 0.0
        ),
        "relaxed_anchored_sequence_coupling_balance_exact_long_range_contact_recall": (
            _rounded(mean(relaxed_long_recalls))
            if relaxed_long_recalls
            else 0.0
        ),
        "relaxed_anchored_sequence_coupling_balance_precision_recall_f1": (
            _rounded(mean(relaxed_f1s)) if relaxed_f1s else 0.0
        ),
        "relaxed_anchored_sequence_coupling_balance_contact_count": (
            sum(relaxed_contact_counts)
        ),
        "relaxed_anchored_sequence_coupling_balance_contact_map_perfect": (
            bool(relaxed_perfect_flags) and all(relaxed_perfect_flags)
        ),
        "coherent_relaxed_sequence_coupling_balance_selector_name": (
            COHERENT_RELAXED_SEQUENCE_COUPLING_BALANCE_SELECTOR
        ),
        "coherent_relaxed_sequence_coupling_balance_claim_allowed": False,
        "coherent_relaxed_sequence_coupling_balance_exact_contact_precision": (
            _rounded(mean(coherent_relaxed_precisions))
            if coherent_relaxed_precisions
            else 0.0
        ),
        "coherent_relaxed_sequence_coupling_balance_exact_long_range_contact_recall": (
            _rounded(mean(coherent_relaxed_long_recalls))
            if coherent_relaxed_long_recalls
            else 0.0
        ),
        "coherent_relaxed_sequence_coupling_balance_precision_recall_f1": (
            _rounded(mean(coherent_relaxed_f1s))
            if coherent_relaxed_f1s
            else 0.0
        ),
        "coherent_relaxed_sequence_coupling_balance_contact_count": (
            sum(coherent_relaxed_contact_counts)
        ),
        "coherent_relaxed_sequence_coupling_balance_contact_map_perfect": (
            bool(coherent_relaxed_perfect_flags)
            and all(coherent_relaxed_perfect_flags)
        ),
        "conflict_shell_density_rescue_selector_name": (
            CONFLICT_SHELL_DENSITY_RESCUE_SELECTOR
        ),
        "conflict_shell_density_rescue_claim_allowed": False,
        "conflict_shell_density_rescue_exact_contact_precision": (
            _rounded(mean(conflict_rescue_precisions))
            if conflict_rescue_precisions
            else 0.0
        ),
        "conflict_shell_density_rescue_exact_long_range_contact_recall": (
            _rounded(mean(conflict_rescue_long_recalls))
            if conflict_rescue_long_recalls
            else 0.0
        ),
        "conflict_shell_density_rescue_precision_recall_f1": (
            _rounded(mean(conflict_rescue_f1s))
            if conflict_rescue_f1s
            else 0.0
        ),
        "conflict_shell_density_rescue_contact_count": (
            sum(conflict_rescue_contact_counts)
        ),
        "conflict_shell_density_rescue_contact_map_perfect": (
            bool(conflict_rescue_perfect_flags)
            and all(conflict_rescue_perfect_flags)
        ),
        "global_contact_map_collapse_selector_name": (
            GLOBAL_CONTACT_MAP_COLLAPSE_SELECTOR
        ),
        "global_contact_map_collapse_claim_allowed": False,
        "global_contact_map_collapse_exact_contact_precision": (
            _rounded(mean(global_collapse_precisions))
            if global_collapse_precisions
            else 0.0
        ),
        "global_contact_map_collapse_exact_long_range_contact_recall": (
            _rounded(mean(global_collapse_long_recalls))
            if global_collapse_long_recalls
            else 0.0
        ),
        "global_contact_map_collapse_precision_recall_f1": (
            _rounded(mean(global_collapse_f1s))
            if global_collapse_f1s
            else 0.0
        ),
        "global_contact_map_collapse_contact_count": (
            sum(global_collapse_contact_counts)
        ),
        "global_contact_map_collapse_contact_map_perfect": (
            bool(global_collapse_perfect_flags)
            and all(global_collapse_perfect_flags)
        ),
        "constructive_gap_voting_selector_name": (
            CONSTRUCTIVE_GAP_VOTING_CONTACT_MAP_SELECTOR
        ),
        "constructive_gap_voting_claim_allowed": False,
        "constructive_gap_voting_exact_contact_precision": (
            _rounded(mean(constructive_gap_precisions))
            if constructive_gap_precisions
            else 0.0
        ),
        "constructive_gap_voting_exact_long_range_contact_recall": (
            _rounded(mean(constructive_gap_long_recalls))
            if constructive_gap_long_recalls
            else 0.0
        ),
        "constructive_gap_voting_precision_recall_f1": (
            _rounded(mean(constructive_gap_f1s))
            if constructive_gap_f1s
            else 0.0
        ),
        "constructive_gap_voting_contact_count": sum(
            constructive_gap_contact_counts
        ),
        "constructive_gap_voting_contact_map_perfect": (
            bool(constructive_gap_perfect_flags)
            and all(constructive_gap_perfect_flags)
        ),
    }


def _cached_contact_universe_ceiling_metrics(
    *,
    rows: Sequence[object],
    dataset: object,
    selected_events: Sequence[NucleusClosureEvent],
) -> dict[str, object]:
    constraints_by_row = dataset.constraints_by_row_id()
    selected_by_row: dict[str, list[NucleusClosureEvent]] = {}
    for event in selected_events:
        selected_by_row.setdefault(event.row_id, []).append(event)

    event_precisions: list[float] = []
    event_long_recalls: list[float] = []
    event_f1s: list[float] = []
    event_contact_counts: list[int] = []
    event_full_recall_flags: list[bool] = []
    external_precisions: list[float] = []
    external_long_recalls: list[float] = []
    external_f1s: list[float] = []
    external_contact_counts: list[int] = []
    external_full_recall_flags: list[bool] = []
    row_trace: list[str] = []
    for row in rows:
        native_pairs = set(row.native_contact_pairs())
        native_long_pairs = {
            pair for pair in native_pairs if pair[1] - pair[0] >= 24
        }
        selected_region_pairs = _valid_long_range_pairs(
            row_length=row.sequence_length,
            pairs={
                pair
                for event in selected_by_row.get(row.row_id, ())
                for pair in event.candidate_region_pairs()
            },
        )
        external_long_pairs = {
            constraint.pair()
            for constraint in constraints_by_row.get(row.row_id, ())
            if constraint.sequence_separation >= 24
        }
        event_scores = _row_contact_scores(
            selected_region_pairs,
            native_pairs=native_pairs,
            native_long_pairs=native_long_pairs,
        )
        external_scores = _row_contact_scores(
            external_long_pairs,
            native_pairs=native_pairs,
            native_long_pairs=native_long_pairs,
        )
        event_precisions.append(float(event_scores["precision"]))
        event_long_recalls.append(float(event_scores["long_recall"]))
        event_f1s.append(float(event_scores["precision_recall_f1"]))
        event_contact_counts.append(int(event_scores["contact_count"]))
        event_full_recall_flags.append(
            float(event_scores["long_recall"]) == 1.0
        )
        external_precisions.append(float(external_scores["precision"]))
        external_long_recalls.append(float(external_scores["long_recall"]))
        external_f1s.append(float(external_scores["precision_recall_f1"]))
        external_contact_counts.append(int(external_scores["contact_count"]))
        external_full_recall_flags.append(
            float(external_scores["long_recall"]) == 1.0
        )
        row_trace.append(
            (
                f"{row.row_id}:event_recall={event_scores['long_recall']}"
                f":external_recall={external_scores['long_recall']}"
            )
        )

    return {
        "event_union_ceiling_exact_contact_precision": (
            _rounded(mean(event_precisions)) if event_precisions else 0.0
        ),
        "event_union_ceiling_exact_long_range_contact_recall": (
            _rounded(mean(event_long_recalls)) if event_long_recalls else 0.0
        ),
        "event_union_ceiling_precision_recall_f1": (
            _rounded(mean(event_f1s)) if event_f1s else 0.0
        ),
        "event_union_ceiling_contact_count": sum(event_contact_counts),
        "event_union_ceiling_full_long_recall": (
            bool(event_full_recall_flags) and all(event_full_recall_flags)
        ),
        "external_long_constraint_exact_contact_precision": (
            _rounded(mean(external_precisions)) if external_precisions else 0.0
        ),
        "external_long_constraint_exact_long_range_contact_recall": (
            _rounded(mean(external_long_recalls))
            if external_long_recalls
            else 0.0
        ),
        "external_long_constraint_precision_recall_f1": (
            _rounded(mean(external_f1s)) if external_f1s else 0.0
        ),
        "external_long_constraint_contact_count": sum(external_contact_counts),
        "external_long_constraint_full_long_recall": (
            bool(external_full_recall_flags)
            and all(external_full_recall_flags)
        ),
        "contact_universe_ceiling_row_trace": ";".join(row_trace),
    }


def score_cached_contact_map_modes_v0(
    *,
    target_manifest: Path,
    selected_events_csv: Path,
    target_rows_csv: Path,
    coupling_dir: Path,
    output_dir: Path,
    selector_name: str = SATURATION_BOUNDARY_SELECTOR,
) -> tuple[Path, Path]:
    manifest = _load_json(target_manifest)
    selected_rows = _read_csv_rows(selected_events_csv)
    target_rows_by_id = {
        str(row.get("target_id", "")): row for row in _read_csv_rows(target_rows_csv)
    }
    output_rows: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []

    for raw_target in manifest.get("targets", []):
        if not isinstance(raw_target, Mapping):
            continue
        target_id = str(raw_target.get("target_id", "")).strip()
        if not target_id:
            continue
        try:
            benchmark_file = _resolve_path(
                raw_target.get("benchmark_file", ""),
                base_dir=BATTERY_REPO_ROOT,
            )
            external_coupling_file = _target_coupling_file(
                raw_target,
                coupling_dir=coupling_dir,
            )
            rows = load_real_coordinate_visual_rows(benchmark_file)
            row_by_id = {row.row_id: row for row in rows}
            import_result = import_external_coupling_dataset(
                rows=rows,
                external_coupling_file=external_coupling_file,
            )
            events = _events_for_target(
                selected_rows,
                target_id=target_id,
                row_by_id=row_by_id,
                selector_name=selector_name,
            )
            metric = _cached_metric(
                target_id=target_id,
                target_row=target_rows_by_id.get(target_id),
                selected_event_count=len(events),
            )
            run = TraceLoopRun(
                selector_name=selector_name,
                dataset=import_result.dataset,
                metric=metric,
                selected_events=events,
                selected_rows=(),
                constraint_count=len(import_result.dataset.constraints),
                control_kind="cached_external_real_selected_events",
            )
            exact_metrics = _exact_contact_metrics(
                rows=rows,
                dataset=import_result.dataset,
                run=run,
            )
            scaffold_metrics = _scaffold_contact_metrics(
                rows=rows,
                dataset=import_result.dataset,
                run=run,
                exact_metrics=exact_metrics,
            )
            row_local_critical_metrics = (
                _cached_row_local_critical_switch_metrics(
                    rows=rows,
                    dataset=import_result.dataset,
                    selected_events=events,
                )
            )
            anchored_sequence_balance_metrics = (
                _cached_anchored_sequence_coupling_balance_metrics(
                    rows=rows,
                    dataset=import_result.dataset,
                    selected_events=events,
                )
            )
            contact_universe_ceiling_metrics = (
                _cached_contact_universe_ceiling_metrics(
                    rows=rows,
                    dataset=import_result.dataset,
                    selected_events=events,
                )
            )
            row_out: dict[str, object] = {
                "target_id": target_id,
                "selector_name": selector_name,
                "selected_event_count": len(events),
                "cached_selected_event_count": len(events),
                "cached_source_selected_events_csv": str(selected_events_csv),
                "cached_source_target_rows_csv": str(target_rows_csv),
                "cache_claim_allowed": False,
            }
            row_out.update(exact_metrics.to_dict())
            row_out.update(scaffold_metrics.to_dict())
            row_out.update(row_local_critical_metrics)
            row_out.update(anchored_sequence_balance_metrics)
            row_out.update(contact_universe_ceiling_metrics)
            row_out.update(_external_coupling_quality_summary(external_coupling_file))
            output_rows.append(row_out)
        except Exception as exc:  # noqa: BLE001
            failures.append(
                {
                    "target_id": target_id,
                    "failure_kind": exc.__class__.__name__,
                    "failure_message": str(exc),
                }
            )

    switch_audit = _apply_self_critical_quality_switch(output_rows)
    for row in output_rows:
        row["row_local_critical_switch_delta_vs_self_critical_quality_switch"] = (
            _rounded(
                float(
                    row.get(
                        "row_local_critical_switch_precision_recall_f1",
                        0.0,
                    )
                )
                - float(
                    row.get(
                        "self_critical_quality_switch_precision_recall_f1",
                        0.0,
                    )
                )
            )
        )
        row["anchored_sequence_coupling_balance_delta_vs_row_local_critical_switch"] = (
            _rounded(
                float(
                    row.get(
                        "anchored_sequence_coupling_balance_precision_recall_f1",
                        0.0,
                    )
                )
                - float(
                    row.get(
                        "row_local_critical_switch_precision_recall_f1",
                        0.0,
                    )
                )
            )
        )
        row["relaxed_anchored_sequence_coupling_balance_delta_vs_anchored_sequence_coupling_balance"] = (
            _rounded(
                float(
                    row.get(
                        "relaxed_anchored_sequence_coupling_balance_precision_recall_f1",
                        0.0,
                    )
                )
                - float(
                    row.get(
                        "anchored_sequence_coupling_balance_precision_recall_f1",
                        0.0,
                    )
                )
            )
        )
        row["coherent_relaxed_sequence_coupling_balance_delta_vs_anchored_sequence_coupling_balance"] = (
            _rounded(
                float(
                    row.get(
                        "coherent_relaxed_sequence_coupling_balance_precision_recall_f1",
                        0.0,
                    )
                )
                - float(
                    row.get(
                        "anchored_sequence_coupling_balance_precision_recall_f1",
                        0.0,
                    )
                )
            )
        )
        row["conflict_shell_density_rescue_delta_vs_coherent_relaxed_sequence_coupling_balance"] = (
            _rounded(
                float(
                    row.get(
                        "conflict_shell_density_rescue_precision_recall_f1",
                        0.0,
                    )
                )
                - float(
                    row.get(
                        "coherent_relaxed_sequence_coupling_balance_precision_recall_f1",
                        0.0,
                    )
                )
            )
        )
        row["global_contact_map_collapse_delta_vs_conflict_shell_density_rescue"] = (
            _rounded(
                float(
                    row.get(
                        "global_contact_map_collapse_precision_recall_f1",
                        0.0,
                    )
                )
                - float(
                    row.get(
                        "conflict_shell_density_rescue_precision_recall_f1",
                        0.0,
                    )
                )
            )
        )
        row["constructive_gap_voting_delta_vs_global_contact_map_collapse"] = (
            _rounded(
                float(
                    row.get(
                        "constructive_gap_voting_precision_recall_f1",
                        0.0,
                    )
                )
                - float(
                    row.get(
                        "global_contact_map_collapse_precision_recall_f1",
                        0.0,
                    )
                )
            )
        )
    report = {
        "report_kind": CACHED_MODE_SCORE_REPORT_KIND,
        "source_report_kind": BATTERY_REPORT_KIND,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "target_manifest": str(target_manifest),
        "selected_events_csv": str(selected_events_csv),
        "target_rows_csv": str(target_rows_csv),
        "selector_name": selector_name,
        "target_count": len(output_rows),
        "failure_count": len(failures),
        "cache_claim_allowed": False,
        "classification": "cached_mode_score_inconclusive",
        "self_critical_quality_switch": switch_audit,
        "mean_phase_coverage_scaffold_precision_recall_f1": _mean_field(
            output_rows,
            "phase_coverage_scaffold_precision_recall_f1",
        ),
        "mean_region_density_top_l_scaffold_precision_recall_f1": _mean_field(
            output_rows,
            "region_density_top_l_scaffold_precision_recall_f1",
        ),
        "mean_phase_density_conflict_shell_scaffold_precision_recall_f1": _mean_field(
            output_rows,
            "phase_density_conflict_shell_scaffold_precision_recall_f1",
        ),
        "mean_phase_density_conflict_phase_confidence_scaffold_precision_recall_f1": _mean_field(
            output_rows,
            "phase_density_conflict_phase_confidence_scaffold_precision_recall_f1",
        ),
        "mean_self_critical_quality_switch_precision_recall_f1": _mean_field(
            output_rows,
            "self_critical_quality_switch_precision_recall_f1",
        ),
        "mean_row_local_critical_switch_precision_recall_f1": _mean_field(
            output_rows,
            "row_local_critical_switch_precision_recall_f1",
        ),
        "mean_row_local_critical_switch_delta_vs_self_critical_quality_switch": _mean_field(
            output_rows,
            "row_local_critical_switch_delta_vs_self_critical_quality_switch",
        ),
        "self_critical_quality_switch_contact_map_perfect_rate": _mean_field(
            output_rows,
            "self_critical_quality_switch_contact_map_perfect",
        ),
        "row_local_critical_switch_contact_map_perfect_rate": _mean_field(
            output_rows,
            "row_local_critical_switch_contact_map_perfect",
        ),
        "mean_anchored_sequence_coupling_balance_precision_recall_f1": _mean_field(
            output_rows,
            "anchored_sequence_coupling_balance_precision_recall_f1",
        ),
        "mean_anchored_sequence_coupling_balance_delta_vs_row_local_critical_switch": _mean_field(
            output_rows,
            "anchored_sequence_coupling_balance_delta_vs_row_local_critical_switch",
        ),
        "anchored_sequence_coupling_balance_contact_map_perfect_rate": _mean_field(
            output_rows,
            "anchored_sequence_coupling_balance_contact_map_perfect",
        ),
        "mean_relaxed_anchored_sequence_coupling_balance_precision_recall_f1": _mean_field(
            output_rows,
            "relaxed_anchored_sequence_coupling_balance_precision_recall_f1",
        ),
        "mean_relaxed_anchored_sequence_coupling_balance_delta_vs_anchored_sequence_coupling_balance": _mean_field(
            output_rows,
            "relaxed_anchored_sequence_coupling_balance_delta_vs_anchored_sequence_coupling_balance",
        ),
        "relaxed_anchored_sequence_coupling_balance_contact_map_perfect_rate": _mean_field(
            output_rows,
            "relaxed_anchored_sequence_coupling_balance_contact_map_perfect",
        ),
        "mean_coherent_relaxed_sequence_coupling_balance_precision_recall_f1": _mean_field(
            output_rows,
            "coherent_relaxed_sequence_coupling_balance_precision_recall_f1",
        ),
        "mean_coherent_relaxed_sequence_coupling_balance_delta_vs_anchored_sequence_coupling_balance": _mean_field(
            output_rows,
            "coherent_relaxed_sequence_coupling_balance_delta_vs_anchored_sequence_coupling_balance",
        ),
        "coherent_relaxed_sequence_coupling_balance_contact_map_perfect_rate": _mean_field(
            output_rows,
            "coherent_relaxed_sequence_coupling_balance_contact_map_perfect",
        ),
        "mean_conflict_shell_density_rescue_precision_recall_f1": _mean_field(
            output_rows,
            "conflict_shell_density_rescue_precision_recall_f1",
        ),
        "mean_conflict_shell_density_rescue_delta_vs_coherent_relaxed_sequence_coupling_balance": _mean_field(
            output_rows,
            "conflict_shell_density_rescue_delta_vs_coherent_relaxed_sequence_coupling_balance",
        ),
        "conflict_shell_density_rescue_contact_map_perfect_rate": _mean_field(
            output_rows,
            "conflict_shell_density_rescue_contact_map_perfect",
        ),
        "mean_global_contact_map_collapse_precision_recall_f1": _mean_field(
            output_rows,
            "global_contact_map_collapse_precision_recall_f1",
        ),
        "mean_global_contact_map_collapse_delta_vs_conflict_shell_density_rescue": _mean_field(
            output_rows,
            "global_contact_map_collapse_delta_vs_conflict_shell_density_rescue",
        ),
        "global_contact_map_collapse_contact_map_perfect_rate": _mean_field(
            output_rows,
            "global_contact_map_collapse_contact_map_perfect",
        ),
        "mean_constructive_gap_voting_precision_recall_f1": _mean_field(
            output_rows,
            "constructive_gap_voting_precision_recall_f1",
        ),
        "mean_constructive_gap_voting_delta_vs_global_contact_map_collapse": _mean_field(
            output_rows,
            "constructive_gap_voting_delta_vs_global_contact_map_collapse",
        ),
        "constructive_gap_voting_contact_map_perfect_rate": _mean_field(
            output_rows,
            "constructive_gap_voting_contact_map_perfect",
        ),
        "mean_event_union_ceiling_exact_long_range_contact_recall": _mean_field(
            output_rows,
            "event_union_ceiling_exact_long_range_contact_recall",
        ),
        "mean_external_long_constraint_exact_long_range_contact_recall": _mean_field(
            output_rows,
            "external_long_constraint_exact_long_range_contact_recall",
        ),
        "event_union_ceiling_full_long_recall_rate": _mean_field(
            output_rows,
            "event_union_ceiling_full_long_recall",
        ),
        "external_long_constraint_full_long_recall_rate": _mean_field(
            output_rows,
            "external_long_constraint_full_long_recall",
        ),
        "folding_problem_solved": False,
        "failures": failures,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "cached_contact_map_mode_score_report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    rows_path = write_csv_rows(
        output_rows,
        output_dir / "cached_contact_map_mode_score_rows.csv",
    )
    return report_path, rows_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Score contact-map completion modes from a frozen selected-events "
            "CSV without rerunning external-coupling frontier selection."
        )
    )
    parser.add_argument("--target-manifest", default=str(DEFAULT_TARGET_MANIFEST))
    parser.add_argument("--selected-events", default=str(DEFAULT_SELECTED_EVENTS))
    parser.add_argument("--target-rows", default=str(DEFAULT_TARGET_ROWS))
    parser.add_argument("--coupling-dir", default=str(DEFAULT_COUPLING_DIR))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--selector-name", default=SATURATION_BOUNDARY_SELECTOR)
    args = parser.parse_args()
    for output in score_cached_contact_map_modes_v0(
        target_manifest=Path(args.target_manifest),
        selected_events_csv=Path(args.selected_events),
        target_rows_csv=Path(args.target_rows),
        coupling_dir=Path(args.coupling_dir),
        output_dir=Path(args.output_dir),
        selector_name=args.selector_name,
    ):
        print(output)


if __name__ == "__main__":
    main()
