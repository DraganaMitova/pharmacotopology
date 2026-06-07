from __future__ import annotations

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_coupling_nucleus_selector import (
    build_coupling_nucleus_context,
    select_coupling_trace_loop_aggressive_relative_frontier_events,
)
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows
from pharmacotopology.folding_self_consistent_contact_loop import run_self_consistent_contact_loop



def test_self_consistent_loop_respects_native_free_self_decision_without_static_threshold() -> None:
    rows = load_real_coordinate_visual_rows(REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json")
    row = next(item for item in rows if item.source_accession == "4AKE:A")
    coupling_dataset = load_coupling_dataset(
        REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
    )
    context = build_coupling_nucleus_context(rows=rows, coupling_dataset=coupling_dataset)
    events = tuple(
        event
        for event in select_coupling_trace_loop_aggressive_relative_frontier_events(context)
        if event.row_id == row.row_id
    )
    constraints = coupling_dataset.constraints_by_row_id()[row.row_id]

    packet = run_self_consistent_contact_loop(
        row=row,
        events=events,
        constraints=constraints,
        event_source="aggressive_relative_frontier",
    )
    report = packet.report.to_dict()

    assert report["kind"] == "self_consistent_contact_loop_v0"
    assert report["source_accession"] == "4AKE:A"
    assert report["event_source"] == "aggressive_relative_frontier"
    assert report["candidate_event_count"] == 1078
    assert report["candidate_pair_count"] == 19824
    assert report["external_coupling_pair_count"] == 197
    assert report["seed_pair_count"] == 197
    assert report["final_pair_count"] == 197
    assert report["decision_rule"].endswith("no_static_confidence_threshold")
    assert report["negative_control_count"] == 3
    assert report["real_rank_among_controls"] == 1
    assert report["self_consistency_status"] == "self_consistency_survived_matched_negative_controls"
    assert report["self_consistent_internal_claim_allowed"] is True
    assert report["external_independent_claim_allowed"] is False
    assert report["global_folding_claim_allowed"] is False
    assert report["folding_problem_solved"] is False

    assert report["coordinate_truth_used_before_selection"] is False
    assert report["native_truth_used_before_selection"] is False
    assert report["raw_sequence_exposed"] is False
    assert report["native_truth_attached_after_selection_for_evaluation"] is True

    assert report["contact_precision_after_native_audit"] == 0.147208
    assert report["long_range_recall_after_native_audit"] == 0.093264
    assert report["contact_precision_delta_vs_seed_after_native_audit"] == -0.030457
    assert report["long_range_recall_delta_vs_seed_after_native_audit"] == -0.005182
    assert report["final_contact_map_hash"] == (
        "sha256:877c5b19e38f3a702d1932ae2376ff210957debb5ef9c336b49b1815c6451887"
    )
