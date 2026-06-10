from __future__ import annotations

from pathlib import Path

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_global_verifier_annealing import (
    GLOBAL_VERIFIER_MODE,
    degree_coherence_score,
    global_structure_score,
    globally_filter_contacts,
    loop_patch_coherence_score,
    run_global_verifier_annealing_packet,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows

ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
EXTERNAL_COUPLINGS = ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"


def _rows():
    rows = load_real_coordinate_visual_rows(BENCHMARK)
    return tuple(row for row in rows if row.source_accession in {"1PGA:A", "1CSP:A"})


def test_global_coherence_scores_are_native_free_and_bounded() -> None:
    patch = [(10, 60), (11, 61), (12, 62), (13, 63)]
    isolated = [(10, 60), (25, 90), (40, 130)]
    assert 0.0 <= degree_coherence_score(patch, 90) <= 1.0
    assert loop_patch_coherence_score(patch) > loop_patch_coherence_score(isolated)


def test_global_contact_filter_enforces_degree_caps_without_truth() -> None:
    row = _rows()[0]
    n = row.sequence_length
    coords = [[float(i) * 2.0, 0.0, 0.0] for i in range(n)]
    raw = [(1, j) for j in range(4, min(n, 40))]
    scores = {pair: 0.9 for pair in raw}
    selected, selected_scores = globally_filter_contacts(row, coords, raw, scores, restraints=())
    assert selected
    assert len(selected) < len(raw)
    assert all(0.0 <= score <= 1.0 for score in selected_scores.values())


def test_global_verifier_packet_is_bounded_and_claim_locked_for_small_control() -> None:
    constraints = load_coupling_dataset(EXTERNAL_COUPLINGS).constraints
    packet = run_global_verifier_annealing_packet(
        _rows(),
        constraints=constraints,
        md_steps=8,
        candidate_count=2,
        max_direct=16,
        max_sequence_closure=12,
        max_geodesic=12,
    )

    assert packet.source_mode == GLOBAL_VERIFIER_MODE
    assert packet.global_verifier_included is True
    assert packet.statistical_potential_surrogate_included is True
    assert packet.degree_distribution_coherence_included is True
    assert packet.loop_patch_coherence_included is True
    assert packet.monte_carlo_global_candidate_selection_included is True
    assert packet.contacts_predicted_before_structure is False
    assert packet.coordinate_truth_used_before_selection is False
    assert packet.native_truth_used_before_selection is False
    assert packet.structure_model_used_before_selection is False
    assert packet.learned_geometry_prior_used_before_selection is False
    assert packet.global_verifier_claim_allowed is False
    assert packet.universal_physical_law_claim_allowed is False
    assert packet.folding_problem_solved is False
    assert packet.decisions
