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
    assert expanded["collapsed_contact_precision"] >= 0.35
    assert expanded["collapsed_long_range_precision"] >= 0.35
    assert expanded["collapsed_long_range_recall"] >= 0.12
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

    assert len(accepted) == 3
    assert all(row["self_collapse_profile"] == "direct_ridge_trace" for row in accepted)
    assert any(
        row["self_collapse_tier_mode"] == "direct_ridge_trace_low_score_tier"
        for row in accepted
    )
    assert all(float(row["self_collapse_confidence"]) > 0.0 for row in accepted)
    assert all(float(row["self_collapse_phase_specific_confidence"]) > 0.0 for row in accepted)
    assert all(float(row["self_collapse_degree_consistency"]) > 0.0 for row in accepted)
    assert all(row["native_truth_used_before_frontier_expansion"] is False for row in accepted)
    assert all(row["coordinate_truth_used_before_frontier_expansion"] is False for row in accepted)


def test_4ake_expansion_cutoff_is_gap_based_not_fixed_confidence_threshold():
    report = _diagnostics()
    rows = report["self_deciding_frontier_expansion_rows"]
    candidates = [
        row for row in rows
        if row["source_accession"] == "4AKE:A"
        and row["seed_event"] is False
        and row["self_verified_candidate_profile_allowed"] is True
    ]

    assert candidates
    assert {row["self_verified_frontier_expansion_cutoff_reason"] for row in candidates} == {
        "self_collapse_verified_internal_gap_frontier"
    }
    assert all("distribution" not in row["self_verified_frontier_expansion_cutoff_reason"] for row in candidates)
    assert all(
        row["self_verified_frontier_expansion_selected"] is row["self_verified_frontier_expansion_gap_selected"]
        for row in candidates
    )
    accepted = [row for row in candidates if row["self_verified_frontier_expansion_selected"] is True]
    rejected = [row for row in candidates if row["self_verified_frontier_expansion_selected"] is False]
    natural_boundary = accepted[-1]["self_verified_frontier_expansion_natural_boundary"]
    assert float(natural_boundary) == min(
        float(row["self_verified_frontier_expansion_identity_normalized_score"])
        for row in accepted
    )
    assert max(
        float(row["self_verified_frontier_expansion_identity_normalized_score"])
        for row in rejected
    ) < float(natural_boundary)
    assert all(
        row["self_verified_frontier_expansion_normalization_reason"]
        == "identity_normalized_seed_baseline"
        for row in candidates
    )


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


def test_expansion_uses_identity_normalized_scores_not_raw_global_confidence():
    report = _diagnostics()
    rows = report["self_deciding_frontier_expansion_rows"]
    candidates = [
        row for row in rows
        if row["source_accession"] == "4AKE:A"
        and row["seed_event"] is False
        and row["self_verified_candidate_profile_allowed"] is True
    ]

    assert candidates
    assert all(float(row["self_verified_frontier_expansion_identity_baseline"]) > 0.0 for row in candidates)
    assert all(float(row["self_verified_frontier_expansion_identity_normalized_score"]) > 0.0 for row in candidates)
    assert all(row["native_truth_used_before_frontier_expansion"] is False for row in candidates)
    assert all(row["coordinate_truth_used_before_frontier_expansion"] is False for row in candidates)


def test_4ake_frontier_ceiling_audit_blocks_false_fully_solved_claim():
    report = _diagnostics()
    audit = report["frontier_ceiling_audit"]["4AKE:A"]
    expansion = report["hard_target_rescue_probe"]["4AKE:A"][
        "self_deciding_frontier_expansion_merged"
    ]

    assert audit["native_truth_used_before_frontier_ceiling_audit"] is False
    assert audit["native_truth_attached_after_frontier_ceiling_for_evaluation"] is True
    assert audit["expanded_region_long_range_recall_ceiling"] >= expansion["collapsed_long_range_recall"]
    assert audit["competitive_region_long_range_recall_ceiling"] < 0.40
    assert audit["candidate_generator_long_range_recall_ceiling"] == 1.0
    assert audit["frontier_generation_bottleneck_for_0_40_recall"] is True
