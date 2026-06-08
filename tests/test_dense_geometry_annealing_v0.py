from __future__ import annotations

from pathlib import Path

from pharmacotopology.folding_dense_geometry_annealing import (
    EXTERNAL_DCA_MULTIFIELD_ANNEALING_MODE,
    PURE_SEQUENCE_COARSE_ANNEALING_MODE,
    run_dense_geometry_annealing_packet,
    run_dense_geometry_annealing_row,
)
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows

ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
EXTERNAL_COUPLINGS = ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"


def _small_rows():
    rows = load_real_coordinate_visual_rows(BENCHMARK)
    return tuple(row for row in rows if row.source_accession in {"1PGA:A", "1CSP:A"})


def test_pure_sequence_coarse_annealing_is_native_free_and_nonempty() -> None:
    row = _small_rows()[0]
    report, factors, solutions, decisions = run_dense_geometry_annealing_row(
        row=row,
        constraints=(),
        source_mode=PURE_SEQUENCE_COARSE_ANNEALING_MODE,
        relaxation_steps=10,
        ensemble_size=2,
        pair_pool_multiplier=3.0,
        max_selected_contact_multiplier=1.3,
    )

    assert report.selected_contact_count > 0
    assert factors
    assert solutions
    assert decisions
    assert report.coordinate_truth_used_before_selection is False
    assert report.native_truth_used_before_selection is False
    assert report.structure_model_used_before_selection is False
    assert report.learned_geometry_prior_used_before_selection is False
    assert report.msa_dca_used_before_selection is False
    assert report.row_universal_physical_law_claim_allowed is False


def test_external_multifield_annealing_runs_without_coordinate_native_or_learned_taint() -> None:
    rows = _small_rows()
    constraints = load_coupling_dataset(EXTERNAL_COUPLINGS).constraints
    packet = run_dense_geometry_annealing_packet(
        rows,
        constraints=constraints,
        source_mode=EXTERNAL_DCA_MULTIFIELD_ANNEALING_MODE,
        relaxation_steps=10,
        ensemble_size=2,
        pair_pool_multiplier=3.0,
        max_selected_contact_multiplier=1.3,
    )

    assert packet.row_count == len(rows)
    assert packet.coordinate_truth_used_before_selection is False
    assert packet.native_truth_used_before_selection is False
    assert packet.structure_model_used_before_selection is False
    assert packet.learned_geometry_prior_used_before_selection is False
    assert packet.msa_dca_used_before_selection is True
    assert packet.universal_physical_law_claim_allowed is False
    assert packet.folding_problem_solved is False
    assert packet.decisions


def test_external_multifield_records_both_selection_channels() -> None:
    row = _small_rows()[0]
    constraints = load_coupling_dataset(EXTERNAL_COUPLINGS).constraints
    report, _factors, _solutions, decisions = run_dense_geometry_annealing_row(
        row=row,
        constraints=constraints,
        source_mode=EXTERNAL_DCA_MULTIFIELD_ANNEALING_MODE,
        relaxation_steps=8,
        ensemble_size=2,
        pair_pool_multiplier=3.0,
        max_selected_contact_multiplier=1.3,
    )

    assert report.compact_channel_selected_count > 0
    assert report.long_range_channel_selected_count > 0
    assert any(decision.selected_by_compact_channel for decision in decisions)
    assert any(decision.selected_by_long_range_channel for decision in decisions)
