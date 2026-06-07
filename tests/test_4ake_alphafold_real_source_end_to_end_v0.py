from __future__ import annotations

import json
import os
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
RUNNER = REPO_ROOT / "scripts" / "run_independent_contact_ensemble_probe_v0.py"


def test_4ake_alphafold_real_source_cracks_without_leakage_or_timeout(tmp_path: Path) -> None:
    """The one benchmark that matters now: real AlphaFold source rescues 4AKE.

    This is intentionally the only active regression test in the cleaned suite.
    It executes the actual probe into a temporary output directory, so it does not
    merely trust committed JSON artifacts.  The subprocess timeout prevents the
    old slow-test behavior from returning silently.
    """

    raw_pdb = PDB_PATH.read_bytes()
    provenance = json.loads(PROVENANCE_PATH.read_text(encoding="utf-8"))

    assert provenance["model_identity"]["alphafold_model_id"] == "AF-P69441-F1-model_v4"
    assert provenance["file"]["sha256"] == sha256(raw_pdb).hexdigest()
    assert provenance["file"]["size_bytes"] == len(raw_pdb)
    assert provenance["safety_and_leakage"]["coordinate_truth_used_before_selection"] is False
    assert provenance["safety_and_leakage"]["native_truth_used_before_selection"] is False

    report_path = tmp_path / "report.json"
    decisions_path = tmp_path / "decisions.csv"
    selected_path = tmp_path / "selected_contacts.csv"
    evidence_path = tmp_path / "evidence.json"

    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT / "src")

    subprocess.run(
        [
            sys.executable,
            str(RUNNER),
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
            "--report-output",
            str(report_path),
            "--decisions-output",
            str(decisions_path),
            "--selected-output",
            str(selected_path),
            "--evidence-output",
            str(evidence_path),
        ],
        cwd=REPO_ROOT,
        env=env,
        check=True,
        timeout=30,
        capture_output=True,
        text=True,
    )

    payload = json.loads(report_path.read_text(encoding="utf-8"))
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
