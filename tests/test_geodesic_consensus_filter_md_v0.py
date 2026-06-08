from __future__ import annotations

from pathlib import Path

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_geodesic_consensus_filter_md import (
    GEODESIC_CONSENSUS_FILTER_MODE,
    build_geodesic_consensus_filter_constraints,
    consensus_geodesic_contacts,
    estimate_bending_energy,
    estimate_geodesic_distance,
    filter_geodesic_anchor_pair,
    run_geodesic_consensus_filter_packet,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows

ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
EXTERNAL_COUPLINGS = ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"


def _rows():
    rows = load_real_coordinate_visual_rows(BENCHMARK)
    return tuple(row for row in rows if row.source_accession in {"1PGA:A", "1CSP:A"})


def test_geodesic_filters_are_native_free_and_reject_bad_paths() -> None:
    assert estimate_geodesic_distance((10, 100), (20, 91)) == 10
    assert estimate_bending_energy((10, 100), (10, 90)) >= 0.5
    assert filter_geodesic_anchor_pair(
        (10, 100),
        (45, 130),
        max_geodesic_distance=30,
        max_bending_energy=0.5,
    ) == (False, "distance")
    assert filter_geodesic_anchor_pair(
        (10, 100),
        (10, 90),
        max_geodesic_distance=30,
        max_bending_energy=0.5,
    ) == (False, "bending")
    assert filter_geodesic_anchor_pair(
        (10, 100),
        (20, 90),
        max_geodesic_distance=30,
        max_bending_energy=0.5,
    ) == (True, "accepted")
    consensus = consensus_geodesic_contacts((10, 100), (20, 90), sequence_length=140, vote_threshold=3)
    assert consensus
    assert all(1 <= i < j <= 140 for i, j in consensus)


def test_consensus_filter_builder_expands_but_filters_geodesics_without_truth_taint() -> None:
    row = _rows()[0]
    constraints = load_coupling_dataset(EXTERNAL_COUPLINGS).constraints
    restraints, summary, records = build_geodesic_consensus_filter_constraints(
        row,
        constraints,
        max_direct=24,
        max_sequence_closure=24,
        max_geodesic=24,
        max_geodesic_distance=30,
        max_bending_energy=0.5,
        consensus_vote_threshold=3,
    )

    assert restraints
    assert records
    assert summary.total_restraint_count >= summary.direct_dca_anchor_count
    assert summary.candidate_anchor_count_before_geodesic >= summary.direct_dca_anchor_count
    assert summary.accepted_geodesic_contact_count <= 24
    assert summary.rejected_by_distance_filter_count >= 0
    assert summary.rejected_by_bending_filter_count >= 0
    assert summary.rejected_by_consensus_filter_count >= 0
    assert summary.coordinate_truth_used_before_selection is False
    assert summary.native_truth_used_before_selection is False
    assert all(not c.coordinate_truth_used_to_build_constraint for c in restraints)
    assert all(not c.native_truth_used_before_coupling_selection for c in restraints)


def test_geodesic_consensus_filter_packet_is_bounded_and_claim_locked_for_small_control() -> None:
    constraints = load_coupling_dataset(EXTERNAL_COUPLINGS).constraints
    packet = run_geodesic_consensus_filter_packet(
        _rows(),
        constraints=constraints,
        md_steps=8,
        restarts=1,
        max_direct=16,
        max_sequence_closure=12,
        max_geodesic=12,
        max_geodesic_distance=30,
        max_bending_energy=0.5,
        consensus_vote_threshold=3,
    )

    assert packet.source_mode == GEODESIC_CONSENSUS_FILTER_MODE
    assert packet.direct_global_structure_generation_included is True
    assert packet.contacts_predicted_before_structure is False
    assert packet.geodesic_distance_filter_included is True
    assert packet.bending_energy_filter_included is True
    assert packet.consensus_geodesic_filter_included is True
    assert packet.coordinate_truth_used_before_selection is False
    assert packet.native_truth_used_before_selection is False
    assert packet.structure_model_used_before_selection is False
    assert packet.learned_geometry_prior_used_before_selection is False
    assert packet.geodesic_consensus_filter_claim_allowed is False
    assert packet.universal_physical_law_claim_allowed is False
    assert packet.folding_problem_solved is False
    assert packet.md_packet.decisions
