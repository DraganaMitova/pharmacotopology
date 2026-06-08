from pathlib import Path

from pharmacotopology.folding_five_axis_physics import (
    build_five_axis_contact_decisions,
    matched_control_pairs,
    run_five_axis_challenge,
    selected_pairs_from_decisions,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    load_real_coordinate_visual_rows,
)


DATA = Path("data/folding_real_coordinate_visual_8.locked.json")


def _rows():
    return load_real_coordinate_visual_rows(DATA)


def test_five_axis_decisions_cover_missing_axes_without_leakage():
    row = _rows()[0]
    decisions = build_five_axis_contact_decisions(row, max_sequence_separation=48)

    assert decisions
    assert any(decision.selected for decision in decisions)
    first = decisions[0]
    assert first.energy_support_score >= 0.0
    assert first.entropy_retention_score >= 0.0
    assert first.cooperative_neighbour_score >= 0.0
    assert first.environmental_context_score >= 0.0
    assert first.dynamic_ensemble_score >= 0.0
    assert first.free_energy_support_score >= 0.0
    assert first.coordinate_truth_used_before_selection is False
    assert first.native_truth_used_before_selection is False
    assert first.learned_prior_used_before_selection is False
    assert first.msa_used_before_selection is False
    assert first.raw_sequence_exposed is False


def test_matched_controls_preserve_selected_separation_profile():
    row = _rows()[0]
    decisions = build_five_axis_contact_decisions(row, max_sequence_separation=48)
    selected = selected_pairs_from_decisions(decisions)
    candidates = tuple(decision.pair() for decision in decisions)
    control = matched_control_pairs(
        row=row,
        selected_pairs=selected,
        candidate_pairs=candidates,
        control_index=1,
    )

    assert control
    assert len(control) == len(selected)
    assert sorted(pair[1] - pair[0] for pair in control) == sorted(
        pair[1] - pair[0] for pair in selected
    )
    assert not (set(control) & set(selected))


def test_five_axis_challenge_refuses_universal_claim_on_locked_subset():
    packet = run_five_axis_challenge(_rows()[:2], max_sequence_separation=48)

    assert packet.entropy_axis_included is True
    assert packet.cooperativity_axis_included is True
    assert packet.dynamics_axis_included is True
    assert packet.context_axis_included is True
    assert packet.independent_physics_axis_included is True
    assert packet.universal_physical_law_claim_allowed is False
    assert packet.folding_problem_solved is False
    assert packet.claim_rejection_reason.startswith(
        "global_folding_claim_rejected_five_axis_gate_failed_for_rows:"
    )
    assert packet.coordinate_truth_used_before_selection is False
    assert packet.native_truth_used_before_selection is False
    assert packet.learned_prior_used_before_selection is False
    assert packet.msa_used_before_selection is False
