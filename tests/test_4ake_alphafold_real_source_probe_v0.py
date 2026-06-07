from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

from pharmacotopology.folding_independent_contact_evidence import (
    contact_evidence_from_predicted_pdb,
    parse_ca_pdb_points,
)
from pharmacotopology.folding_real_coordinate_visual_benchmark import (
    load_real_coordinate_visual_rows,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
PDB_PATH = REPO_ROOT / "data" / "independent_contact_sources" / "AF-P69441-F1-model_v4.pdb"
PROVENANCE_PATH = (
    REPO_ROOT
    / "data"
    / "independent_contact_sources"
    / "AF-P69441-F1-model_v4.provenance.json"
)
REPORT_PATH = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "independent_contact_ensemble_4ake_alphafold_v0.json"
)
BENCHMARK_PATH = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"


def _row_4ake():
    rows = load_real_coordinate_visual_rows(BENCHMARK_PATH)
    return next(row for row in rows if row.source_accession == "4AKE:A")


def test_real_alphafold_4ake_pdb_is_present_and_hash_locked() -> None:
    raw = PDB_PATH.read_bytes()
    provenance = json.loads(PROVENANCE_PATH.read_text(encoding="utf-8"))

    assert PDB_PATH.exists()
    assert provenance["model_identity"]["alphafold_model_id"] == "AF-P69441-F1-model_v4"
    assert provenance["file"]["sha256"] == sha256(raw).hexdigest()
    assert provenance["file"]["size_bytes"] == len(raw)
    assert provenance["safety_and_leakage"]["coordinate_truth_used_before_selection"] is False
    assert provenance["safety_and_leakage"]["native_truth_used_before_selection"] is False


def test_real_alphafold_4ake_pdb_parses_exact_chain_length_without_native_truth() -> None:
    row = _row_4ake()
    points = parse_ca_pdb_points(PDB_PATH, chain_id="A")
    evidence = contact_evidence_from_predicted_pdb(
        row=row,
        pdb_path=PDB_PATH,
        chain_id="A",
        source_id="alphafold_db_AF-P69441-F1-model_v4",
    )

    assert row.sequence_length == 214
    assert len(points) == row.sequence_length
    assert len(evidence) == 600
    assert all(not item.coordinate_truth_used_before_selection for item in evidence)
    assert all(not item.native_truth_used_before_selection for item in evidence)
    assert all(not item.raw_sequence_exposed for item in evidence)


def test_real_alphafold_4ake_probe_cracks_precision_without_leakage() -> None:
    payload = json.loads(REPORT_PATH.read_text(encoding="utf-8"))
    report = payload["report"]
    safety = payload["safety"]

    assert safety["benchmark_claim_allowed"] is True
    assert safety["claim_rejection_reason"] == "none"
    assert safety["coordinate_truth_used_before_selection"] is False
    assert safety["native_truth_used_before_selection"] is False
    assert safety["raw_sequence_exposed"] is False

    assert report["independent_structure_pair_count"] == 600
    assert report["final_long_range_pair_count"] == 235
    assert report["true_positive_long_range_contacts"] == 175
    assert report["long_range_precision"] >= 0.70
    assert report["long_range_recall"] >= 0.90
