from __future__ import annotations

import csv
from pathlib import Path

from pharmacotopology.folding_coupling_nucleus_selector import (
    PRIMARY_CONTACT_COLLAPSE_STRATEGY,
    build_coupling_nucleus_context,
    contact_collapse_results_for_selected_events,
    contact_collapse_summary,
)
from pharmacotopology.folding_event_region_contact_collapse import SELF_DECIDING_STRATEGY_NAME
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    load_real_coordinate_visual_rows,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
FRONTIER_FILE = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_coupling_trace_loop_frontier.csv"


def _frontier_events(context):
    event_ids: list[str] = []
    seen: set[str] = set()
    with FRONTIER_FILE.open(newline="", encoding="utf-8") as handle:
        for item in csv.DictReader(handle):
            event_id = str(item.get("event_id", ""))
            if event_id and event_id not in seen:
                seen.add(event_id)
                event_ids.append(event_id)
    return tuple(
        context.event_by_id[event_id]
        for event_id in event_ids
        if event_id in context.event_by_id
    )


def test_frontier_contact_collapse_is_self_deciding_pipeline_not_oracle_selection():
    rows = tuple(load_real_coordinate_visual_rows(BENCHMARK_FILE))
    dataset = load_coupling_dataset(COUPLING_FILE)
    context = build_coupling_nucleus_context(rows=rows, coupling_dataset=dataset)
    results = contact_collapse_results_for_selected_events(context, _frontier_events(context))
    summary = contact_collapse_summary(results)
    one_cll = summary["contact_collapse_1cll_self_deciding"]

    assert summary["contact_collapse_integrated"] is True
    assert summary["contact_collapse_strategy"] == PRIMARY_CONTACT_COLLAPSE_STRATEGY
    assert PRIMARY_CONTACT_COLLAPSE_STRATEGY == SELF_DECIDING_STRATEGY_NAME
    assert one_cll["collapse_strategy"] == SELF_DECIDING_STRATEGY_NAME
    assert one_cll["selected_event_count"] == 7
    assert one_cll["uncollapsed_region_pair_count"] == 432
    assert one_cll["collapsed_pair_count"] <= 25
    assert one_cll["collapse_reduction_ratio"] >= 0.94
    assert one_cll["collapsed_true_positive_contacts"] >= 12
    assert one_cll["collapsed_contact_precision"] >= 0.48
    assert one_cll["collapsed_long_range_precision"] >= 0.55
    assert one_cll["collapsed_long_range_recall"] >= 0.30
    assert one_cll["collapsed_long_range_f1"] >= 0.39
    assert one_cll["native_truth_used_before_collapse_selection"] is False
    assert one_cll["coordinate_truth_used_before_collapse_selection"] is False


def test_self_deciding_event_cutoffs_are_not_a_fixed_pairs_per_event_budget():
    rows = tuple(load_real_coordinate_visual_rows(BENCHMARK_FILE))
    dataset = load_coupling_dataset(COUPLING_FILE)
    context = build_coupling_nucleus_context(rows=rows, coupling_dataset=dataset)
    results = contact_collapse_results_for_selected_events(context, _frontier_events(context))
    one_cll = next(result for result in results if result.evaluation.source_accession == "1CLL:A")
    summaries = one_cll.event_summaries

    assert {summary.collapse_strategy for summary in summaries} == {SELF_DECIDING_STRATEGY_NAME}
    cutoffs = {summary.applied_cutoff for summary in summaries}
    reasons = {summary.collapse_decision_reason for summary in summaries}
    assert len(cutoffs) >= 3
    assert 0 in cutoffs
    assert "self_closed_no_long_range_candidate_space" in reasons
    assert "self_closed_partial_strip_without_direct_root" in reasons
    assert any(reason.startswith("self_") for reason in reasons)
    assert all(summary.max_pairs_per_event == summary.applied_cutoff for summary in summaries)
    assert all(summary.native_truth_used_before_collapse_selection is False for summary in summaries)
    assert all(summary.coordinate_truth_used_before_collapse_selection is False for summary in summaries)
