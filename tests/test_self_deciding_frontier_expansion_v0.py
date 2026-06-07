from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DIAGNOSTICS_FILE = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "contact_collapse_diagnostics_v0.json"
)
EXPANSION_ROWS_FILE = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "self_deciding_frontier_expansion_rows_v0.csv"
)


def _diagnostics() -> dict[str, object]:
    return json.loads(DIAGNOSTICS_FILE.read_text(encoding="utf-8"))


def test_self_deciding_frontier_expansion_is_native_free_and_recorded():
    report = _diagnostics()
    rows = report["self_deciding_frontier_expansion_rows"]
    assert EXPANSION_ROWS_FILE.exists()
    assert rows
    assert all(row["native_truth_used_before_frontier_expansion"] is False for row in rows)
    assert all(row["coordinate_truth_used_before_frontier_expansion"] is False for row in rows)
    assert any(row["self_deciding_frontier_expansion_selected"] is True for row in rows)


def test_4ake_expansion_exposes_frontier_recall_but_controller_keeps_it_as_probe():
    report = _diagnostics()
    probes = report["hard_target_rescue_probe"]["4AKE:A"]
    base = probes["frontier_self_deciding"]
    expanded = probes["self_deciding_frontier_expansion_merged"]

    assert expanded["selected_event_count"] > base["selected_event_count"]
    assert expanded["uncollapsed_long_range_recall"] > base["uncollapsed_long_range_recall"]
    assert expanded["collapsed_true_positive_contacts"] >= base["collapsed_true_positive_contacts"]
    assert expanded["collapsed_long_range_recall"] >= base["collapsed_long_range_recall"]
    assert expanded["collapsed_contact_precision"] < base["collapsed_contact_precision"]
    assert base["native_truth_used_before_collapse_selection"] is False
    assert expanded["native_truth_used_before_collapse_selection"] is False
