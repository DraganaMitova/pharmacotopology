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


def test_4ake_self_verified_expansion_adds_recall_without_precision_collapse():
    report = _diagnostics()
    probes = report["hard_target_rescue_probe"]["4AKE:A"]
    base = probes["frontier_self_deciding"]
    expanded = probes["self_deciding_frontier_expansion_merged"]

    assert expanded["selected_event_count"] > base["selected_event_count"]
    assert expanded["uncollapsed_long_range_recall"] > base["uncollapsed_long_range_recall"]
    assert expanded["collapsed_true_positive_contacts"] > base["collapsed_true_positive_contacts"]
    assert expanded["collapsed_long_range_recall"] > base["collapsed_long_range_recall"]
    assert expanded["collapsed_contact_precision"] >= 0.40
    assert expanded["collapsed_long_range_precision"] >= 0.40
    assert expanded["collapsed_contact_precision"] < base["collapsed_contact_precision"]
    assert base["native_truth_used_before_collapse_selection"] is False
    assert expanded["native_truth_used_before_collapse_selection"] is False


def test_4ake_expansion_rows_are_self_collapse_verified_not_raw_low_floor():
    report = _diagnostics()
    rows = report["self_deciding_frontier_expansion_rows"]
    accepted = [
        row for row in rows
        if row["source_accession"] == "4AKE:A"
        and row["self_verified_frontier_expansion_selected"] is True
        and row["seed_event"] is False
    ]

    assert len(accepted) == 2
    assert all(row["self_collapse_profile"] == "direct_ridge_trace" for row in accepted)
    assert all(float(row["self_collapse_confidence"]) > 0.0 for row in accepted)
    assert all(row["native_truth_used_before_frontier_expansion"] is False for row in accepted)
    assert all(row["coordinate_truth_used_before_frontier_expansion"] is False for row in accepted)


def test_1mbn_expansion_controller_rejects_precision_destroying_expansion():
    report = _diagnostics()
    probes = report["hard_target_rescue_probe"]["1MBN:A"]
    base = probes["frontier_self_deciding"]
    expanded = probes["self_deciding_frontier_expansion_merged"]

    assert expanded["selected_event_count"] == base["selected_event_count"]
    assert expanded["collapsed_pair_count"] == base["collapsed_pair_count"]
    assert expanded["collapsed_true_positive_contacts"] == base["collapsed_true_positive_contacts"]
    assert expanded["collapsed_contact_precision"] == base["collapsed_contact_precision"]
    assert expanded["native_truth_used_before_collapse_selection"] is False
