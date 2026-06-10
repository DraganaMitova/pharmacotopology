from __future__ import annotations

from pathlib import Path

from pharmacotopology.folding_msa_free_model_ensemble import (
    MSA_FREE_MODEL_ENSEMBLE_REPORT_KIND,
    MSAFreeStructureModelSpec,
    run_msa_free_model_ensemble,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows


REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"
ESMFOLD_4AKE_PDB = REPO_ROOT / "external_msa_free_predictors" / "tryhard_runs_live" / "4ake_esmfold_api.pdb"
ALPHAFOLD_4AKE_PDB = REPO_ROOT / "data" / "independent_contact_sources" / "AF-P69441-F1-model_v4.pdb"


def _row_4ake():
    rows = load_real_coordinate_visual_rows(BENCHMARK_FILE)
    return next(item for item in rows if item.source_accession == "4AKE:A")


def test_msa_free_model_ensemble_solves_4ake_from_real_esmfold_without_native_leakage() -> None:
    row = _row_4ake()
    assert ESMFOLD_4AKE_PDB.exists(), "this package should include the successful Mac-safe ESMFold output"

    packet = run_msa_free_model_ensemble(
        row=row,
        model_specs=(
            MSAFreeStructureModelSpec(
                source_id="esmfold_single_sequence_4ake",
                pdb_path=str(ESMFOLD_4AKE_PDB),
                chain_id="A",
            ),
        ),
        iteration_count=3,
    )

    assert packet.kind == MSA_FREE_MODEL_ENSEMBLE_REPORT_KIND
    assert packet.usable_model_count == 1
    assert packet.folding_problem_solved is True
    assert packet.direct_structure_solved is True
    assert packet.consensus_contact_collapse_solved is True
    assert packet.folding_solution_mode == "single_model_consensus_physics_refined"
    assert packet.direct_best_model_metric.native_contact_precision >= 0.70
    assert packet.direct_best_model_metric.native_contact_recall >= 0.70
    assert packet.consensus_metric.native_contact_precision >= 0.70
    assert packet.consensus_metric.native_contact_recall >= 0.70
    assert packet.native_truth_used_before_selection is False
    assert packet.coordinate_truth_used_before_selection is False
    assert packet.alphafold_used_before_selection is False
    assert packet.msa_used_before_selection is False
    assert packet.template_used_before_selection is False
    assert packet.raw_sequence_exposed is False
    assert all(not decision.native_truth_used_before_selection for decision in packet.decisions)
    assert all(not decision.coordinate_truth_used_before_selection for decision in packet.decisions)


def test_msa_free_model_ensemble_supports_multi_model_consensus_with_duplicate_non_af_sources(tmp_path: Path) -> None:
    row = _row_4ake()
    copy_a = tmp_path / "esmfold_model_a.pdb"
    copy_b = tmp_path / "omegafold_like_model_b.pdb"
    copy_a.write_text(ESMFOLD_4AKE_PDB.read_text(encoding="utf-8"), encoding="utf-8")
    copy_b.write_text(ESMFOLD_4AKE_PDB.read_text(encoding="utf-8"), encoding="utf-8")

    packet = run_msa_free_model_ensemble(
        row=row,
        model_specs=(
            MSAFreeStructureModelSpec(source_id="esmfold_single_sequence_4ake", pdb_path=str(copy_a), chain_id="A"),
            MSAFreeStructureModelSpec(source_id="omegafold_single_sequence_4ake", pdb_path=str(copy_b), chain_id="A"),
        ),
        iteration_count=2,
    )

    assert packet.usable_model_count == 2
    assert packet.folding_problem_solved is True
    assert packet.folding_solution_mode == "multi_model_consensus_physics_refined"
    selected = [decision for decision in packet.decisions if decision.selected]
    assert selected
    assert max(decision.model_vote_count for decision in selected) == 2


def test_msa_free_model_ensemble_rejects_alphafold_like_sources_by_default() -> None:
    row = _row_4ake()
    packet = run_msa_free_model_ensemble(
        row=row,
        model_specs=(
            MSAFreeStructureModelSpec(
                source_id="AF-P69441-F1-model_v4",
                pdb_path=str(ALPHAFOLD_4AKE_PDB),
                chain_id="A",
            ),
        ),
        iteration_count=1,
    )

    assert packet.usable_model_count == 0
    assert packet.rejected_model_count == 1
    assert packet.folding_problem_solved is False
    assert packet.claim_rejection_reason == "missing_usable_msa_free_model_output"
    assert packet.models[0].status == "rejected_alphafold_like_source_for_msa_free_ensemble"
    assert packet.native_truth_used_before_selection is False
    assert packet.coordinate_truth_used_before_selection is False
