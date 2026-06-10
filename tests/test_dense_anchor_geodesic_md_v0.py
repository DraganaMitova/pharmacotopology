from __future__ import annotations

from pathlib import Path

from pharmacotopology.folding_dense_anchor_geodesic_md import (
    DENSE_ANCHOR_GEODESIC_MODE,
    build_dense_anchor_geodesic_constraints,
    run_dense_anchor_geodesic_packet,
)
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows

ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
EXTERNAL_COUPLINGS = ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"


def _rows():
    rows = load_real_coordinate_visual_rows(BENCHMARK)
    return tuple(row for row in rows if row.source_accession in {"1PGA:A", "1CSP:A"})


def test_dense_anchor_builder_expands_safe_dca_without_truth_taint() -> None:
    row = _rows()[0]
    constraints = load_coupling_dataset(EXTERNAL_COUPLINGS).constraints
    dense, summary, records = build_dense_anchor_geodesic_constraints(
        row,
        constraints,
        max_direct=24,
        max_sequence_closure=24,
        max_geodesic=48,
    )

    assert summary.total_dense_anchor_count >= summary.direct_dca_anchor_count
    assert summary.geodesic_anchor_count >= 0
    assert dense
    assert records
    assert summary.coordinate_truth_used_before_selection is False
    assert summary.native_truth_used_before_selection is False
    assert summary.structure_model_used_before_selection is False
    assert summary.learned_geometry_prior_used_before_selection is False
    assert all(not c.coordinate_truth_used_to_build_constraint for c in dense)
    assert all(not c.native_truth_used_before_coupling_selection for c in dense)


def test_dense_anchor_geodesic_packet_is_structure_first_and_bounded() -> None:
    constraints = load_coupling_dataset(EXTERNAL_COUPLINGS).constraints
    packet = run_dense_anchor_geodesic_packet(
        _rows(),
        constraints=constraints,
        md_steps=10,
        restarts=1,
        max_direct=20,
        max_sequence_closure=18,
        max_geodesic=24,
    )

    assert packet.source_mode == DENSE_ANCHOR_GEODESIC_MODE
    assert packet.direct_global_structure_generation_included is True
    assert packet.contacts_predicted_before_structure is False
    assert packet.dense_anchor_expansion_included is True
    assert packet.anchor_chaining_included is True
    assert packet.geodesic_interpolation_included is True
    assert packet.coordinate_truth_used_before_selection is False
    assert packet.native_truth_used_before_selection is False
    assert packet.structure_model_used_before_selection is False
    assert packet.learned_geometry_prior_used_before_selection is False
    assert packet.universal_physical_law_claim_allowed is False
    assert packet.md_packet.decisions


def test_dense_anchor_geodesic_claim_gate_stays_locked_for_small_control() -> None:
    constraints = load_coupling_dataset(EXTERNAL_COUPLINGS).constraints
    packet = run_dense_anchor_geodesic_packet(
        _rows(),
        constraints=constraints,
        md_steps=8,
        restarts=1,
        max_direct=16,
        max_sequence_closure=12,
        max_geodesic=12,
    )

    assert packet.dense_anchor_geodesic_claim_allowed is False
    assert packet.universal_physical_law_claim_allowed is False
    assert packet.folding_problem_solved is False
    assert packet.claim_rejection_reason.startswith("dense_anchor_geodesic_claim_rejected")
