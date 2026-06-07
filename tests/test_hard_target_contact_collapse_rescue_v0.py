from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DIAGNOSTICS_FILE = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "contact_collapse_diagnostics_v0.json"
)


def _hard_target_probe(source_accession: str, strategy: str) -> dict[str, object]:
    report = json.loads(DIAGNOSTICS_FILE.read_text(encoding="utf-8"))
    return report["hard_target_rescue_probe"][source_accession][strategy]


def test_4ake_ridge_rescue_finds_precise_but_low_recall_contacts():
    result = _hard_target_probe("4AKE:A", "ridge_coupling")
    assert result["collapsed_pair_count"] == 4
    assert result["collapsed_true_positive_contacts"] == 3
    assert result["collapsed_contact_precision"] >= 0.75
    assert result["collapsed_long_range_precision"] == 1.0
    assert result["collapsed_long_range_recall"] < 0.02
    assert result["native_truth_used_before_collapse_selection"] is False
    assert result["coordinate_truth_used_before_collapse_selection"] is False


def test_1mbn_recall_probe_recovers_signal_but_remains_noisy():
    result = _hard_target_probe("1MBN:A", "frontier_recall")
    assert result["collapsed_pair_count"] == 80
    assert result["collapsed_true_positive_contacts"] >= 9
    assert result["collapsed_long_range_precision"] >= 0.50
    assert result["frontier_long_native_retention"] >= 0.75
    assert result["collapsed_contact_precision"] < 0.20
    assert result["native_truth_used_before_collapse_selection"] is False
    assert result["coordinate_truth_used_before_collapse_selection"] is False
