from __future__ import annotations

from pathlib import Path

from pharmacotopology.folding_coarse_grain_md_geometry import (
    MULTISTART_DCA_MD_MODE,
    PURE_SEQUENCE_MD_MODE,
    run_coarse_grain_md_geometry_packet,
    run_coarse_grain_md_geometry_row,
)
from pharmacotopology.folding_evolutionary_constraints import load_coupling_dataset
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows

ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
EXTERNAL_COUPLINGS = ROOT / "data" / "folding_real_coordinate_visual_8_external_couplings.v0.locked.json"


def _small_rows():
    rows = load_real_coordinate_visual_rows(BENCHMARK)
    return tuple(row for row in rows if row.source_accession in {"1PGA:A", "1CSP:A"})


def test_pure_sequence_md_geometry_is_structure_first_and_native_free() -> None:
    row = _small_rows()[0]
    report, decisions = run_coarse_grain_md_geometry_row(
        row=row,
        constraints=(),
        source_mode=PURE_SEQUENCE_MD_MODE,
        md_steps=16,
        restarts=1,
    )

    assert report.extracted_contact_count > 0
    assert decisions
    assert report.coordinate_truth_used_before_selection is False
    assert report.native_truth_used_before_selection is False
    assert report.structure_model_used_before_selection is False
    assert report.learned_geometry_prior_used_before_selection is False
    assert report.msa_dca_used_before_selection is False
    assert report.row_universal_physical_law_claim_allowed is False
    assert all(decision.selected_from_final_structure for decision in decisions)


def test_multistart_dca_md_geometry_runs_without_coordinate_native_or_learned_taint() -> None:
    rows = _small_rows()
    constraints = load_coupling_dataset(EXTERNAL_COUPLINGS).constraints
    packet = run_coarse_grain_md_geometry_packet(
        rows,
        constraints=constraints,
        source_mode=MULTISTART_DCA_MD_MODE,
        md_steps=14,
        restarts=2,
        max_restraints=48,
    )

    assert packet.row_count == len(rows)
    assert packet.direct_global_structure_generation_included is True
    assert packet.contacts_predicted_before_structure is False
    assert packet.coordinate_truth_used_before_selection is False
    assert packet.native_truth_used_before_selection is False
    assert packet.structure_model_used_before_selection is False
    assert packet.learned_geometry_prior_used_before_selection is False
    assert packet.msa_dca_used_before_selection is True
    assert packet.universal_physical_law_claim_allowed is False
    assert packet.folding_problem_solved is False
    assert packet.decisions


def test_md_geometry_packet_claim_gate_stays_locked_for_small_control() -> None:
    rows = _small_rows()
    packet = run_coarse_grain_md_geometry_packet(
        rows,
        constraints=(),
        source_mode=PURE_SEQUENCE_MD_MODE,
        md_steps=12,
        restarts=1,
    )

    assert packet.md_geometry_claim_allowed is False
    assert packet.universal_physical_law_claim_allowed is False
    assert packet.folding_problem_solved is False
    assert packet.claim_rejection_reason.startswith("coarse_grain_md_geometry_claim_rejected")
