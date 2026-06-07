from __future__ import annotations

from pathlib import Path

from pharmacotopology.folding_coupling_nucleus_selector import (
    build_coupling_nucleus_context,
    select_coupling_trace_loop_aggressive_relative_frontier_events,
    select_coupling_trace_loop_events,
)
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    load_real_coordinate_visual_rows,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
COUPLING_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"


def _long_native_pairs(row):
    return {pair for pair in row.native_contact_pairs() if pair[1] - pair[0] >= 24}


def _region_ceiling(row, events):
    native_long = _long_native_pairs(row)
    region_pairs: set[tuple[int, int]] = set()
    for event in events:
        region_pairs.update(event.candidate_region_pairs())
    return len(region_pairs & native_long) / len(native_long)


def test_aggressive_relative_frontier_opens_4ake_candidate_pool_without_oracle_inputs():
    rows = tuple(load_real_coordinate_visual_rows(BENCHMARK_FILE))
    dataset = load_coupling_dataset(COUPLING_FILE)
    context = build_coupling_nucleus_context(rows=rows, coupling_dataset=dataset)
    target = next(row for row in rows if row.source_accession == "4AKE:A")

    seed_all = select_coupling_trace_loop_events(context)
    aggressive_all = select_coupling_trace_loop_aggressive_relative_frontier_events(
        context,
        seed_events=seed_all,
    )
    seed = tuple(event for event in seed_all if event.row_id == target.row_id)
    aggressive = tuple(event for event in aggressive_all if event.row_id == target.row_id)

    assert rows.index(target) == 6
    assert len(aggressive) > len(seed) * 1.2
    assert _region_ceiling(target, aggressive) >= _region_ceiling(target, seed)
    assert dataset.coordinate_truth_tainted is False
    assert dataset.native_truth_tainted is False
    assert dataset.structure_model_tainted is False
