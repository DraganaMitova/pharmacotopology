from __future__ import annotations

from pathlib import Path

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_long_range_calibrated_verifier import (
    POTENTIAL_SEP_FLOOR,
    build_leave_one_out_long_range_potential,
    long_range_contact_score,
    run_long_range_calibrated_verifier_packet,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows

ROOT = Path(__file__).resolve().parents[1]
ROWS = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
COUPLINGS = ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"


def test_leave_one_out_long_range_potential_excludes_target() -> None:
    rows = load_real_coordinate_visual_rows(ROWS)
    target = rows[0]
    potential = build_leave_one_out_long_range_potential(rows, target.row_id)
    assert potential.target_row_id == target.row_id
    assert potential.reference_row_count == len(rows) - 1
    assert potential.target_native_excluded_from_potential is True
    assert potential.possible_long_range_pair_count > 0
    assert 0.0 < potential.base_probability < 1.0


def test_long_range_contact_score_is_sequence_only_and_bounded() -> None:
    rows = load_real_coordinate_visual_rows(ROWS)
    target = rows[0]
    potential = build_leave_one_out_long_range_potential(rows, target.row_id)
    score = long_range_contact_score((1, min(target.sequence_length, POTENTIAL_SEP_FLOOR + 2)), target.sequence, potential)
    assert 0.0 <= score <= 1.0
    short_score = long_range_contact_score((1, 10), target.sequence, potential)
    assert short_score == 0.0


def test_long_range_calibrated_packet_keeps_truth_boundary() -> None:
    rows = load_real_coordinate_visual_rows(ROWS)
    constraints = load_coupling_dataset(COUPLINGS).constraints
    packet = run_long_range_calibrated_verifier_packet(
        rows,
        constraints=constraints,
        evaluation_source_accessions=("4AKE:A",),
        md_steps=8,
        candidate_count=1,
        max_direct=10,
        max_sequence_closure=12,
        max_geodesic=20,
    )
    assert packet.row_count == 1
    assert packet.long_range_reference_potential_included is True
    assert packet.leave_one_target_out_long_range_calibration is True
    assert packet.target_native_excluded_before_selection_for_all_rows is True
    assert packet.universal_physical_law_claim_allowed is False
    row = packet.rows[0]
    assert row.source_accession == "4AKE:A"
    assert row.native_truth_used_before_target_selection is False
    assert row.coordinate_truth_used_before_target_selection is False
    assert row.target_native_excluded_from_long_range_potential is True
