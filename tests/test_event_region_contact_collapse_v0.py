from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

from pharmacotopology.folding_contact_law_features import (
    contact_law_feature_rows_for_row,
)
from pharmacotopology.folding_coupling_nucleus_selector import build_coupling_nucleus_context
from pharmacotopology.folding_event_region_contact_collapse import (
    EVENT_REGION_CONTACT_COLLAPSE_BOUNDARY,
    collapse_row_event_regions,
)
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    load_real_coordinate_visual_rows,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
FRONTIER_FILE = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_coupling_trace_loop_frontier.csv"


@lru_cache(maxsize=1)
def _one_cll_frontier_events():
    rows = tuple(load_real_coordinate_visual_rows(BENCHMARK_FILE))
    row = next(candidate for candidate in rows if candidate.source_accession == "1CLL:A")
    dataset = load_coupling_dataset(COUPLING_FILE)
    context = build_coupling_nucleus_context(rows=rows, coupling_dataset=dataset)
    event_ids: list[str] = []
    with FRONTIER_FILE.open(newline="", encoding="utf-8") as handle:
        for item in csv.DictReader(handle):
            if item.get("row_id") == row.row_id:
                event_ids.append(str(item["event_id"]))
    events = tuple(context.event_by_id[event_id] for event_id in event_ids)
    row_constraints = tuple(dataset.constraints_by_row_id().get(row.row_id, ()))
    row_features = contact_law_feature_rows_for_row(row)
    return row, events, row_features, row_constraints


def test_event_region_contact_collapse_boundary_is_clean():
    assert EVENT_REGION_CONTACT_COLLAPSE_BOUNDARY == (
        "sequence_and_external_coupling_only_before_native_evaluation"
    )
    row, events, row_features, row_constraints = _one_cll_frontier_events()
    result = collapse_row_event_regions(
        row=row,
        events=events,
        row_features=row_features,
        row_constraints=row_constraints,
        collapse_strategy="frontier_precision",
    )
    assert result.evaluation.native_truth_used_before_collapse_selection is False
    assert result.evaluation.coordinate_truth_used_before_collapse_selection is False
    assert all(
        pair.native_truth_used_before_collapse_selection is False
        and pair.coordinate_truth_used_before_collapse_selection is False
        and pair.raw_sequence_exposed is False
        for pair in result.collapsed_pairs
    )


def test_1cll_frontier_precision_collapse_reduces_region_without_claiming_perfect_map():
    row, events, row_features, row_constraints = _one_cll_frontier_events()
    result = collapse_row_event_regions(
        row=row,
        events=events,
        row_features=row_features,
        row_constraints=row_constraints,
        collapse_strategy="frontier_precision",
    )
    assert len(events) == 7
    assert result.evaluation.uncollapsed_region_pair_count == 432
    assert result.evaluation.uncollapsed_true_positive_contacts == 21
    assert result.evaluation.uncollapsed_region_precision == 0.048611
    assert result.evaluation.collapsed_pair_count <= 35
    assert result.evaluation.collapsed_contact_precision >= 0.30
    assert result.evaluation.precision_improvement_factor >= 6.0
    assert result.evaluation.collapsed_contact_recall < 1.0


def test_1cll_frontier_recall_collapse_preserves_long_frontier_signal():
    row, events, row_features, row_constraints = _one_cll_frontier_events()
    result = collapse_row_event_regions(
        row=row,
        events=events,
        row_features=row_features,
        row_constraints=row_constraints,
        collapse_strategy="frontier_recall",
    )
    assert result.evaluation.collapsed_pair_count < result.evaluation.uncollapsed_region_pair_count
    assert result.evaluation.frontier_native_retention >= 0.80
    assert result.evaluation.frontier_long_native_retention == 1.0
    assert (
        result.evaluation.collapsed_long_range_recall
        == result.evaluation.uncollapsed_long_range_recall
    )
    assert (
        result.evaluation.collapsed_contact_precision
        > result.evaluation.uncollapsed_region_precision
    )


def test_1cll_frontier_balanced_collapse_enters_main_tradeoff_band():
    row, events, row_features, row_constraints = _one_cll_frontier_events()
    result = collapse_row_event_regions(
        row=row,
        events=events,
        row_features=row_features,
        row_constraints=row_constraints,
        collapse_strategy="frontier_balanced",
    )
    assert len(events) == 7
    assert 35 <= result.evaluation.collapsed_pair_count <= 50
    assert result.evaluation.collapsed_contact_precision >= 0.30
    assert result.evaluation.collapsed_long_range_precision >= 0.40
    assert result.evaluation.collapsed_long_range_recall >= 0.30
    assert result.evaluation.collapsed_long_range_f1 >= 0.35
    assert result.evaluation.precision_improvement_factor >= 6.0
    assert result.evaluation.native_truth_used_before_collapse_selection is False
    assert result.evaluation.coordinate_truth_used_before_collapse_selection is False


def test_1cll_internal_gap_balanced_collapse_cracks_contact_level_tradeoff():
    row, events, row_features, row_constraints = _one_cll_frontier_events()
    result = collapse_row_event_regions(
        row=row,
        events=events,
        row_features=row_features,
        row_constraints=row_constraints,
        collapse_strategy="frontier_internal_gap_balanced",
    )
    assert len(events) == 7
    assert result.evaluation.uncollapsed_region_pair_count == 432
    assert result.evaluation.collapsed_pair_count <= 25
    assert result.evaluation.collapse_reduction_ratio >= 0.94
    assert result.evaluation.collapsed_true_positive_contacts >= 13
    assert result.evaluation.collapsed_contact_precision >= 0.55
    assert result.evaluation.collapsed_long_range_precision >= 0.65
    assert result.evaluation.collapsed_long_range_recall >= 0.36
    assert result.evaluation.collapsed_long_range_f1 >= 0.47
    assert result.evaluation.frontier_long_native_retention == 1.0
    assert result.evaluation.precision_improvement_factor >= 10.0
    assert result.evaluation.native_truth_used_before_collapse_selection is False
    assert result.evaluation.coordinate_truth_used_before_collapse_selection is False
    reasons = {summary.collapse_decision_reason for summary in result.event_summaries}
    assert "first_gap_lock" in reasons
    assert "short_sparse_region_closed" in reasons
    assert any(summary.support_completion_pair_count for summary in result.event_summaries)
    assert any(summary.edge_rescue_pair_count for summary in result.event_summaries)
