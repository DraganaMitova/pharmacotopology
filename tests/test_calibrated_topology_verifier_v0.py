from __future__ import annotations

from pathlib import Path

from pharmacotopology.folding_calibrated_topology_verifier import (
    CALIBRATED_VERIFIER_MODE,
    calibrated_filter_contacts,
    optional_energy_backend_status,
    run_calibrated_topology_verifier_packet,
)
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows

ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
EXTERNAL_COUPLINGS = ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"


def _rows():
    rows = load_real_coordinate_visual_rows(BENCHMARK)
    return tuple(row for row in rows if row.source_accession in {"1PGA:A", "1CSP:A", "1UBQ:A"})


def test_optional_backend_probe_is_explicit_and_non_fatal() -> None:
    status = optional_energy_backend_status()
    assert "pyrosetta_available" in status
    assert "dfire_available" in status
    assert status["backend_used"] == "transparent_leave_one_out_feature_discriminator"


def test_calibrated_filter_is_bounded_and_degree_capped() -> None:
    row = _rows()[0]
    coords = [[float(i) * 2.0, 0.0, 0.0] for i in range(row.sequence_length)]
    raw = [(1, j) for j in range(4, min(row.sequence_length, 60))]
    scores = {pair: 0.8 for pair in raw}
    selected, selected_scores = calibrated_filter_contacts(row, coords, raw, scores, restraints=(), weights=(0.1,) * 9, probability_threshold=0.1)
    assert selected
    assert len(selected) < len(raw)
    assert all(0.0 <= score <= 1.0 for score in selected_scores.values())


def test_calibrated_packet_is_leave_one_out_and_claim_locked_for_small_control() -> None:
    constraints = load_coupling_dataset(EXTERNAL_COUPLINGS).constraints
    packet = run_calibrated_topology_verifier_packet(
        _rows(),
        constraints=constraints,
        md_steps=8,
        candidate_count=1,
        max_direct=10,
        max_sequence_closure=8,
        max_geodesic=8,
    )
    assert packet.source_mode == CALIBRATED_VERIFIER_MODE
    assert packet.calibrated_structural_feature_discriminator_included is True
    assert packet.leave_one_target_out_calibration is True
    assert packet.target_native_excluded_before_selection_for_all_rows is True
    assert packet.pdb_calibration_native_contacts_used_from_other_rows is True
    assert packet.universal_physical_law_claim_allowed is False
    assert packet.folding_problem_solved is False
    assert packet.decisions
