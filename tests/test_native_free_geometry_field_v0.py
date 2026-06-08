from __future__ import annotations

from pathlib import Path

from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_native_free_geometry_field import (
    EXTERNAL_DCA_GEOMETRY_FIELD_MODE,
    ORACLE_ANCHOR_GEOMETRY_FIELD_CONTROL_MODE,
    PURE_SEQUENCE_SYMBOLIC_GEOMETRY_MODE,
    build_native_free_geometry_scores,
    run_native_free_geometry_field_packet,
    run_native_free_geometry_field_row,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows


ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
EXTERNAL_COUPLINGS = ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"
ORACLE_COUPLINGS = ROOT / "data" / "folding_real_coordinate_visual_8_couplings.locked.json"


def _small_rows():
    rows = load_real_coordinate_visual_rows(BENCHMARK)
    return tuple(row for row in rows if row.source_accession in {"1PGA:A", "1CSP:A"})


def test_pure_sequence_geometry_scores_are_native_free_and_nonempty() -> None:
    row = _small_rows()[0]
    scores, components, anchor_count, coordinate_taint, native_taint, structure_model, dca_used = build_native_free_geometry_scores(
        row=row,
        constraints=(),
        source_mode=PURE_SEQUENCE_SYMBOLIC_GEOMETRY_MODE,
    )

    assert scores
    assert components
    assert anchor_count == 0
    assert coordinate_taint is False
    assert native_taint is False
    assert structure_model is False
    assert dca_used is False
    assert max(scores.values()) > 0.0


def test_external_geometry_field_runs_without_coordinate_or_native_taint() -> None:
    rows = _small_rows()
    constraints = load_coupling_dataset(EXTERNAL_COUPLINGS).constraints
    packet = run_native_free_geometry_field_packet(
        rows,
        constraints=constraints,
        source_mode=EXTERNAL_DCA_GEOMETRY_FIELD_MODE,
        ensemble_size=3,
        pair_pool_multiplier=4.0,
        max_selected_contact_multiplier=1.4,
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


def test_oracle_geometry_field_control_is_tainted_and_never_claimable() -> None:
    row = _small_rows()[0]
    constraints = load_coupling_dataset(ORACLE_COUPLINGS).constraints
    report, _factors, _solutions, _decisions = run_native_free_geometry_field_row(
        row=row,
        constraints=constraints,
        source_mode=ORACLE_ANCHOR_GEOMETRY_FIELD_CONTROL_MODE,
        ensemble_size=2,
        pair_pool_multiplier=3.0,
        max_selected_contact_multiplier=1.2,
    )

    assert report.coordinate_truth_used_before_selection is True
    assert report.native_truth_used_before_selection is True
    assert report.row_geometry_field_claim_allowed is False
    assert report.row_universal_physical_law_claim_allowed is False
    assert report.row_claim_rejection_reason == "oracle_or_coordinate_native_tainted_anchor_control_not_claimable"
