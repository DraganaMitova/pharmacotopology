from __future__ import annotations

import json
import subprocess
import sys
from hashlib import sha256
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PDB_PATH = REPO_ROOT / "data" / "independent_contact_sources" / "AF-P69441-F1-model_v4.pdb"
PROVENANCE_PATH = (
    REPO_ROOT
    / "data"
    / "independent_contact_sources"
    / "AF-P69441-F1-model_v4.provenance.json"
)
DONE_MODE_RUNNER = REPO_ROOT / "scripts" / "run_done_mode_contact_prediction_v0.py"


def test_4ake_done_mode_cracks_with_visual_proof_without_leakage_or_timeout(tmp_path: Path) -> None:
    """The one active regression test: 4AKE real AF source -> claim + GIF.

    It runs the actual DONE-mode wrapper, which calls the ensemble probe and then
    renders the visual proof GIF.  The timeout keeps the suite from returning to
    the old slow/irrelevant-test state.
    """

    raw_pdb = PDB_PATH.read_bytes()
    provenance = json.loads(PROVENANCE_PATH.read_text(encoding="utf-8"))

    assert provenance["model_identity"]["alphafold_model_id"] == "AF-P69441-F1-model_v4"
    assert provenance["file"]["sha256"] == sha256(raw_pdb).hexdigest()
    assert provenance["file"]["size_bytes"] == len(raw_pdb)
    assert provenance["safety_and_leakage"]["coordinate_truth_used_before_selection"] is False
    assert provenance["safety_and_leakage"]["native_truth_used_before_selection"] is False

    subprocess.run(
        [
            sys.executable,
            str(DONE_MODE_RUNNER),
            "--source-accession",
            "4AKE:A",
            "--predicted-pdb",
            str(PDB_PATH.relative_to(REPO_ROOT)),
            "--predicted-pdb-chain",
            "A",
            "--predicted-source-id",
            "alphafold_db_AF-P69441-F1-model_v4",
            "--event-source",
            "aggressive_relative_frontier",
            "--min-votes",
            "2",
            "--out-dir",
            str(tmp_path),
            "--render-visual-proof",
        ],
        cwd=REPO_ROOT,
        check=True,
        timeout=30,
        capture_output=True,
        text=True,
    )

    done_manifest = json.loads((tmp_path / "done_mode_manifest.json").read_text(encoding="utf-8"))
    payload = json.loads((tmp_path / "contact_prediction_report.json").read_text(encoding="utf-8"))
    gif_manifest = json.loads((tmp_path / "4ake_visual_proof_manifest.json").read_text(encoding="utf-8"))
    gif_path = tmp_path / "4ake_visual_proof.gif"

    assert gif_path.exists()
    assert gif_path.stat().st_size > 50_000
    assert done_manifest["kind"] == "done_mode_contact_prediction_v0"
    assert done_manifest["status"] == "claim_allowed"
    assert done_manifest["done_contract"] == {
        "abstains_instead_of_guessing_without_independent_evidence": True,
        "claim_requires_independent_evidence": True,
        "native_or_coordinate_truth_before_selection_allowed": False,
        "predictor_runs_for_every_locked_target": True,
        "raw_sequence_exposure_allowed": False,
    }

    report = payload["report"]
    safety = payload["safety"]
    assert payload["row_id"] == "coord_007_pdb_4AKE_A_adenylate_kinase"
    assert payload["source_accession"] == "4AKE:A"
    assert payload["event_source"] == "aggressive_relative_frontier"
    assert payload["selected_event_count"] == 1078
    assert payload["source_family_pair_counts"] == {
        "candidate_region": 19824,
        "external_coupling": 214,
        "independent_structure": 600,
    }
    assert payload["selected_contact_map_hash"] == (
        "sha256:51ab5218c8f3dfc587f17014912cec649897382f8f960051b457531fcfa19ff5"
    )

    assert safety == {
        "benchmark_claim_allowed": True,
        "claim_rejection_reason": "none",
        "coordinate_truth_used_before_selection": False,
        "native_truth_attached_after_selection_for_evaluation": True,
        "native_truth_used_before_selection": False,
        "raw_sequence_exposed": False,
    }

    assert report["benchmark_claim_allowed"] is True
    assert report["claim_rejection_reason"] == "none"
    assert report["coordinate_truth_used_before_selection"] is False
    assert report["native_truth_used_before_selection"] is False
    assert report["raw_sequence_exposed"] is False

    assert report["independent_structure_pair_count"] == 600
    assert report["external_coupling_pair_count"] == 214
    assert report["candidate_region_pair_count"] == 19824
    assert report["final_pair_count"] == 302
    assert report["final_long_range_pair_count"] == 235
    assert report["true_positive_long_range_contacts"] == 175

    assert report["contact_precision"] == 0.774834
    assert report["contact_recall"] == 0.418605
    assert report["long_range_precision"] == 0.744681
    assert report["long_range_recall"] == 0.906736

    assert gif_manifest["kind"] == "4ake_alphafold_visual_proof_gif_v0"
    assert gif_manifest["visual_boundary"] == {
        "coordinate_truth_used_before_selection": False,
        "native_truth_used_before_selection": False,
        "native_truth_used_for_rendering_only_after_selection": True,
    }
    assert gif_manifest["pair_counts"] == {
        "candidate_region": 19824,
        "external_coupling": 214,
        "false_positive_long_range_pairs": 60,
        "independent_structure": 600,
        "missed_native_long_range_pairs": 18,
        "selected_long_range_pairs": 235,
        "selected_pairs": 302,
        "true_positive_long_range_pairs": 175,
    }
