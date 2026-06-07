from __future__ import annotations

import sys
from pathlib import Path

from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows
from pharmacotopology.folding_single_sequence_structure_source import (
    is_alphafold_like_source_id,
    run_single_sequence_prediction_source,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def _row_4ake():
    rows = load_real_coordinate_visual_rows(REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json")
    return next(item for item in rows if item.source_accession == "4AKE:A")


def test_msa_free_adapter_rejects_alphafold_like_source_id_without_opt_in(tmp_path: Path) -> None:
    row = _row_4ake()
    assert is_alphafold_like_source_id("AF-P69441-F1-model_v4") is True
    run = run_single_sequence_prediction_source(
        row=row,
        predicted_pdb_path=tmp_path / "prediction.pdb",
        predictor_source_id="AF-P69441-F1-model_v4",
        prediction_command=f"{sys.executable} -c 'raise SystemExit(99)'",
        timeout_seconds=2,
    )

    assert run.source_rejected is True
    assert run.prediction_command_executed is False
    assert run.status == "rejected_alphafold_like_source_id_for_msa_free_probe"
    assert run.output_pdb_exists is False


def test_msa_free_adapter_times_out_and_does_not_persist_query_fasta_by_default(tmp_path: Path) -> None:
    row = _row_4ake()
    run = run_single_sequence_prediction_source(
        row=row,
        predicted_pdb_path=tmp_path / "prediction.pdb",
        predictor_source_id="custom_single_sequence_model",
        prediction_command=f"{sys.executable} -c 'import time; time.sleep(5)'",
        timeout_seconds=1,
        working_directory=tmp_path,
    )

    assert run.prediction_command_executed is True
    assert run.prediction_command_timed_out is True
    assert run.status == "prediction_command_timed_out"
    assert run.output_pdb_exists is False
    assert run.query_fasta_persisted is False
    assert run.raw_sequence_exposed_in_persisted_artifacts is False
    assert not (tmp_path / "query.fasta").exists()
