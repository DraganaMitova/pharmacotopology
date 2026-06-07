"""
1CLL precision check.

The old 7/7 statement is event-level: all seven frontier segment events pass
adaptive gating. That does not automatically mean an exact/perfect contact map.
This test keeps the 7/7 gate check, then verifies the contact-map precision
claim separately.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


def _candidate_region_pairs(row: pd.Series) -> set[tuple[int, int]]:
    return {
        (left, right)
        for left in range(int(row["segment_a_start"]), int(row["segment_a_end"]) + 1)
        for right in range(int(row["segment_b_start"]), int(row["segment_b_end"]) + 1)
        if right - left >= 3
    }


def test_1cll_has_7_frontier_events_but_not_perfect_contact_map() -> None:
    manifest_path = (
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "external_coupling_target_manifest_v0.json"
    )
    with manifest_path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)

    frontier_path = (
        ROOT
        / "first_contact_clean_pharmacotopology_layer_run"
        / "external_coupling_trace_loop_frontier.csv"
    )
    frontier = pd.read_csv(frontier_path)
    target_ids_1cll = [entry["row_id"] for entry in manifest if entry.get("pdb_id") == "1CLL"]
    frontier_1cll = frontier[frontier["row_id"].isin(target_ids_1cll)].sort_values(
        "future_preservation_score",
        ascending=False,
    )

    assert len(frontier_1cll) == 7

    # Same event-level adaptive frontier check as before.
    pass_count = int(
        (
            (frontier_1cll["future_preservation_score"] >= 0.48)
            & (frontier_1cll["blocked_future_pressure"] <= 0.20)
        ).sum()
    )
    assert pass_count == 7

    benchmark_rows = load_real_coordinate_visual_rows(
        ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
    )
    cll_row = next(row for row in benchmark_rows if row.source_accession == "1CLL:A")
    native_pairs = set(cll_row.native_contact_pairs())
    predicted_pairs = set().union(
        *(_candidate_region_pairs(row) for _, row in frontier_1cll.iterrows())
    )
    true_positive_pairs = predicted_pairs & native_pairs
    precision = len(true_positive_pairs) / len(predicted_pairs)
    recall = len(true_positive_pairs) / len(native_pairs)

    assert precision < 1.0
    assert recall < 1.0
    assert not (precision == 1.0 and recall == 1.0)
