from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V83_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V83"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v83_certificate_passes_complexity_token_heuristic_scaling_audit() -> None:
    cert = _read(V83_ROOT / "v83_complexity_token_heuristic_scaling_audit_certificate.json")

    assert cert["status"] == "V83_COMPLEXITY_TOKEN_HEURISTIC_SCALING_AUDIT_PASSED"
    assert cert["engine_version_used"] == "E77"
    assert cert["baseline_engine_version"] == "E76"
    assert cert["targets_total"] == 4000
    assert cert["row_family_counts"] == {
        "abstention_taxonomy_control": 500,
        "blind_candidate": 2000,
        "heuristic_ablation_control": 500,
        "sentence_composition_sentinel": 500,
        "token_soup_adversarial_control": 500,
    }
    assert cert["failed_accepted_count"] == 0
    assert cert["token_only_acceptance_count"] == 0
    assert cert["sentinels_preserved"] == 500
    assert cert["overconservative_abstains_reported"] == 100
    assert cert["registry_export_deterministic"] is True
    assert cert["registry_hashes"] == cert["repeated_registry_hashes"]
    assert cert["repeated_run_certificate_core_hash_stable"] is True
    assert cert["heuristic_ablation_unsupported_flip_count"] == 0
    assert cert["coordinate_native_truth_stays_sealed"] is True
    assert cert["no_static_thresholds_used"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["failed_controls"] == []


def test_v83_rows_expose_abstention_and_token_soup_without_accepting() -> None:
    rows = _read(V83_ROOT / "v83_complexity_token_heuristic_scaling_audit_report.json")["rows"]

    assert len(rows) == 4000
    token_rows = [row for row in rows if row["row_family"] == "token_soup_adversarial_control"]
    ablation_rows = [row for row in rows if row["row_family"] == "heuristic_ablation_control"]
    abstain_rows = [row for row in rows if row["row_family"] == "abstention_taxonomy_control"]

    assert len(token_rows) == 500
    assert len(ablation_rows) == 500
    assert len(abstain_rows) == 500
    assert all(row["token_hit_role"] == "evidence_proposal_only" for row in token_rows)
    assert all(row["cannot_directly_accept"] is True for row in token_rows)
    assert all(row["token_only_acceptance"] is False for row in rows)
    assert all(row["failed_accepted"] is False for row in rows)
    assert all(row["heuristic_ablation_flipped_supported_claim"] is False for row in rows)
    assert any(row.get("overconservative_abstain_reported") is True for row in abstain_rows)
    assert all(row["coordinate_native_leakage"] is False for row in rows)
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["protein_folding_solved"] is False for row in rows)
