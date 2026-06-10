from __future__ import annotations

from pathlib import Path

from pharmacotopology.folding_global_factor_graph_ensemble import (
    LEARNED_GEOMETRY_MODE,
    SEQUENCE_ONLY_MODE,
    GLOBAL_FACTOR_GRAPH_KIND,
    GlobalFactorGraphModelSpec,
    run_global_factor_graph_packet,
    run_global_factor_graph_row,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows


REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
ESMFOLD_4AKE_PDB = REPO_ROOT / "external_msa_free_predictors" / "tryhard_runs_live" / "4ake_esmfold_api.pdb"


def _rows():
    return load_real_coordinate_visual_rows(BENCHMARK_FILE)


def _row_4ake():
    return next(row for row in _rows() if row.source_accession == "4AKE:A")


def test_sequence_only_factor_graph_builds_global_groups_without_native_leakage():
    row = _rows()[0]
    report, factors, mutex_groups, solutions, decisions, accepted_specs = run_global_factor_graph_row(
        row=row,
        source_mode=SEQUENCE_ONLY_MODE,
        max_sequence_separation=48,
        ensemble_size=3,
    )

    assert report.source_mode == SEQUENCE_ONLY_MODE
    assert report.all_or_none_factor_count == len(factors)
    assert factors
    assert all(factor.factor_type == "all_or_none" for factor in factors)
    assert solutions
    assert decisions
    assert accepted_specs == ()
    assert report.coordinate_truth_used_before_selection is False
    assert report.native_truth_used_before_selection is False
    assert report.learned_geometry_prior_used_before_selection is False
    assert all(not factor.coordinate_truth_used_before_selection for factor in factors)
    assert all(not decision.native_truth_used_before_selection for decision in decisions)
    # Mutex constraints are generated from candidates, not from native truth.  A
    # small row/subset may have zero explicit residue over-degree groups, but the
    # row still records the mutex gate and the greedy solver enforces it.
    assert report.mutex_group_count == len(mutex_groups)


def test_sequence_only_global_factor_graph_refuses_universal_claim_on_locked_subset():
    packet = run_global_factor_graph_packet(
        _rows()[:2],
        source_mode=SEQUENCE_ONLY_MODE,
        max_sequence_separation=48,
        ensemble_size=3,
    )

    assert packet.kind == GLOBAL_FACTOR_GRAPH_KIND
    assert packet.factor_graph_included is True
    assert packet.all_or_none_factors_included is True
    assert packet.mutex_constraints_included is True
    assert packet.top_k_ensemble_included is True
    assert packet.learned_geometry_prior_used_before_selection is False
    assert packet.factor_graph_ensemble_claim_allowed is False
    assert packet.universal_physical_law_claim_allowed is False
    assert packet.folding_problem_solved is False
    assert packet.claim_rejection_reason.startswith(
        "sequence_only_global_factor_graph_claim_rejected_for_rows:"
    )
    assert packet.coordinate_truth_used_before_selection is False
    assert packet.native_truth_used_before_selection is False
    assert packet.alphafold_used_before_selection is False
    assert packet.msa_used_before_selection is False
    assert packet.template_used_before_selection is False


def test_factor_graph_solves_4ake_when_given_real_non_af_esmfold_geometry_prior():
    row = _row_4ake()
    assert ESMFOLD_4AKE_PDB.exists()
    packet = run_global_factor_graph_packet(
        [row],
        source_mode=LEARNED_GEOMETRY_MODE,
        model_specs_by_accession={
            "4AKE:A": [
                GlobalFactorGraphModelSpec(
                    source_id="esmfold_single_sequence_4ake",
                    pdb_path=str(ESMFOLD_4AKE_PDB),
                    chain_id="A",
                )
            ]
        },
        ensemble_size=4,
    )

    assert packet.kind == GLOBAL_FACTOR_GRAPH_KIND
    assert packet.source_mode == LEARNED_GEOMETRY_MODE
    assert packet.learned_geometry_prior_used_before_selection is True
    assert packet.factor_graph_ensemble_claim_allowed is True
    assert packet.universal_physical_law_claim_allowed is False
    assert packet.folding_problem_solved is True
    assert packet.mean_native_contact_precision_after_audit >= 0.70
    assert packet.mean_native_contact_recall_after_audit >= 0.70
    assert packet.mean_long_range_contact_recall_after_audit >= 0.70
    assert packet.rows[0].row_claim_rejection_reason == (
        "target_solved_by_global_factor_graph_over_non_alphafold_learned_geometry_prior"
    )
    assert packet.coordinate_truth_used_before_selection is False
    assert packet.native_truth_used_before_selection is False
    assert packet.alphafold_used_before_selection is False
    assert packet.msa_used_before_selection is False
    assert packet.template_used_before_selection is False
    assert packet.raw_sequence_exposed is False
