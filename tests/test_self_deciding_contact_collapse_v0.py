from __future__ import annotations

import json
from pathlib import Path

from pharmacotopology.folding_event_region_contact_collapse import SELF_DECIDING_STRATEGY_NAME

REPO_ROOT = Path(__file__).resolve().parents[1]
DIAGNOSTICS_FILE = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "contact_collapse_diagnostics_v0.json"
)


def _diagnostics() -> dict[str, object]:
    return json.loads(DIAGNOSTICS_FILE.read_text(encoding="utf-8"))


def test_1cll_self_deciding_cracks_contact_level_without_static_event_budget():
    report = _diagnostics()
    result = report["one_cll_self_deciding_report"]
    assert result["collapse_strategy"] == SELF_DECIDING_STRATEGY_NAME
    assert result["uncollapsed_region_pair_count"] == 432
    assert result["collapsed_pair_count"] == 25
    assert result["collapsed_true_positive_contacts"] >= 12
    assert result["collapsed_contact_precision"] >= 0.48
    assert result["collapsed_long_range_precision"] >= 0.55
    assert result["collapsed_long_range_recall"] >= 0.30
    assert result["collapsed_long_range_f1"] >= 0.39
    assert result["precision_improvement_factor"] >= 9.0
    assert result["native_truth_used_before_collapse_selection"] is False
    assert result["coordinate_truth_used_before_collapse_selection"] is False


def test_1mbn_self_deciding_keeps_recall_signal_but_removes_recall_probe_noise():
    report = _diagnostics()
    self_result = report["hard_target_rescue_probe"]["1MBN:A"][SELF_DECIDING_STRATEGY_NAME]
    recall_probe = report["hard_target_rescue_probe"]["1MBN:A"]["frontier_recall"]
    assert self_result["collapsed_pair_count"] == 16
    assert recall_probe["collapsed_pair_count"] == 80
    assert self_result["collapsed_true_positive_contacts"] >= recall_probe["collapsed_true_positive_contacts"]
    assert self_result["collapsed_long_range_recall"] >= recall_probe["collapsed_long_range_recall"]
    assert self_result["collapsed_contact_precision"] >= 0.55
    assert self_result["collapsed_contact_precision"] > recall_probe["collapsed_contact_precision"] * 4.0
    assert self_result["frontier_long_native_retention"] >= 0.75
    assert self_result["native_truth_used_before_collapse_selection"] is False
    assert self_result["coordinate_truth_used_before_collapse_selection"] is False


def test_4ake_phase_aware_self_deciding_rescues_precision_but_not_full_recall():
    report = _diagnostics()
    ridge_probe = report["hard_target_rescue_probe"]["4AKE:A"]["ridge_coupling"]
    self_result = report["hard_target_rescue_probe"]["4AKE:A"][SELF_DECIDING_STRATEGY_NAME]
    expansion = report["hard_target_rescue_probe"]["4AKE:A"][
        "self_deciding_frontier_expansion_merged"
    ]

    assert ridge_probe["collapsed_contact_precision"] >= 0.75
    assert ridge_probe["collapsed_long_range_precision"] == 1.0
    assert ridge_probe["collapsed_long_range_recall"] < 0.02

    # The self-deciding scorer now chooses a direct-ridge surface and then widens
    # only broad ridge traces through an internal low-score tier.  This improves
    # recall without accepting the whole expanded frontier, but the row is still
    # not solved because total long-range recall remains below 0.10.
    assert self_result["collapsed_pair_count"] == 34
    assert self_result["collapsed_true_positive_contacts"] >= 16
    assert self_result["collapsed_contact_precision"] >= 0.45
    assert self_result["collapsed_long_range_precision"] >= 0.45
    assert self_result["collapsed_long_range_recall"] >= 0.08
    assert self_result["collapsed_long_range_recall"] < 0.10
    assert self_result["collapsed_long_range_f1"] >= 0.14

    # Native-free frontier expansion exposes more uncollapsed recall, but the
    # collapse controller does not accept it as a solved contact map because
    # contact precision drops sharply and collapsed recall does not materially
    # beat the self-deciding seed frontier.
    assert expansion["uncollapsed_long_range_recall"] > self_result["uncollapsed_long_range_recall"]
    assert expansion["collapsed_long_range_recall"] >= self_result["collapsed_long_range_recall"]
    assert expansion["collapsed_contact_precision"] < self_result["collapsed_contact_precision"]

    assert self_result["native_truth_used_before_collapse_selection"] is False
    assert self_result["coordinate_truth_used_before_collapse_selection"] is False
