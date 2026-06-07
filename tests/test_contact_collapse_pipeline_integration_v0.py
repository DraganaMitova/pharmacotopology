from __future__ import annotations

import csv
from pathlib import Path

from pharmacotopology.folding_coupling_nucleus_selector import (
    PRIMARY_CONTACT_COLLAPSE_STRATEGY,
    build_coupling_nucleus_context,
    contact_collapse_results_for_selected_events,
    contact_collapse_summary,
)
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    load_real_coordinate_visual_rows,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
FRONTIER_FILE = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "external_coupling_trace_loop_frontier.csv"


def test_frontier_contact_collapse_is_pipeline_summary_not_oracle_selection():
    rows = tuple(load_real_coordinate_visual_rows(BENCHMARK_FILE))
    dataset = load_coupling_dataset(COUPLING_FILE)
    context = build_coupling_nucleus_context(rows=rows, coupling_dataset=dataset)
    event_ids: list[str] = []
    seen: set[str] = set()
    with FRONTIER_FILE.open(newline="", encoding="utf-8") as handle:
        for item in csv.DictReader(handle):
            event_id = str(item.get("event_id", ""))
            if event_id and event_id not in seen:
                seen.add(event_id)
                event_ids.append(event_id)
    events = tuple(
        context.event_by_id[event_id]
        for event_id in event_ids
        if event_id in context.event_by_id
    )
    results = contact_collapse_results_for_selected_events(context, events)
    summary = contact_collapse_summary(results)
    one_cll = summary["contact_collapse_1cll_internal_gap"]

    assert summary["contact_collapse_integrated"] is True
    assert summary["contact_collapse_strategy"] == PRIMARY_CONTACT_COLLAPSE_STRATEGY
    assert PRIMARY_CONTACT_COLLAPSE_STRATEGY == "frontier_internal_gap_balanced"
    assert one_cll["selected_event_count"] == 7
    assert one_cll["uncollapsed_region_pair_count"] == 432
    assert one_cll["collapsed_pair_count"] == 23
    assert one_cll["collapse_reduction_ratio"] >= 0.94
    assert one_cll["collapsed_true_positive_contacts"] >= 13
    assert one_cll["collapsed_contact_precision"] >= 0.55
    assert one_cll["collapsed_long_range_precision"] >= 0.65
    assert one_cll["collapsed_long_range_recall"] >= 0.36
    assert one_cll["collapsed_long_range_f1"] >= 0.47
    assert one_cll["native_truth_used_before_collapse_selection"] is False
    assert one_cll["coordinate_truth_used_before_collapse_selection"] is False
