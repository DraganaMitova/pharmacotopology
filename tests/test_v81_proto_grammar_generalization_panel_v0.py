from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V81_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V81"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v81_certificate_passes_generalization_without_lexicon_inflation() -> None:
    cert = _read(V81_ROOT / "v81_proto_grammar_generalization_panel_certificate.json")

    assert cert["status"] == "V81_PROTO_GRAMMAR_GENERALIZATION_PANEL_PASSED"
    assert cert["engine_version_used"] == "E75"
    assert cert["baseline_engine_version"] == "E74"
    assert cert["proto_grammar_total"] == 45
    assert cert["rows_per_proto"] == 24
    assert cert["panel_rows_total"] == 1080
    assert cert["row_role_counts"] == {
        "enemy_grammar_candidate": 225,
        "fresh_usage_context_candidate": 225,
        "merge_negative_candidate": 225,
        "metadata_masked_control": 135,
        "sentinel_neighbor_control": 135,
        "sequence_counterfactual_control": 135,
    }
    assert cert["classification_counts"] == {
        "crystallized_grammar": 8,
        "keep_as_proto_unknown": 21,
        "merge_into_existing_word": 15,
        "retire_as_context_artifact": 1,
    }
    assert cert["crystallized_grammar_count"] == 8
    assert cert["merged_after_generalization_count"] == 15
    assert cert["kept_proto_unknown_count"] == 21
    assert cert["retired_context_artifact_count"] == 1
    assert cert["failed_accepted_count"] == 0
    assert cert["sentinel_regressions"] == 0
    assert cert["static_thresholds_used"] is False
    assert cert["no_static_thresholds_used"] is True
    assert cert["no_forced_expected_labels_used"] is True
    assert cert["v79_v80_context_used_as_seed_not_proof"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["failed_controls"] == []


def test_v81_rows_are_fresh_generalization_controls_without_failed_accepts() -> None:
    report = _read(V81_ROOT / "v81_proto_grammar_generalization_panel_report.json")
    rows = report["rows"]
    crystallization_rows = report["crystallization_rows"]

    assert len(crystallization_rows) == 45
    assert len(rows) == 1080
    assert all(row["fresh_target_context"] is True for row in rows)
    assert all(row["v79_v80_context_reused_as_proof"] is False for row in rows)
    assert all(row["failed_accepted"] is False for row in rows)
    assert all(row["sentinel_regression"] is False for row in rows)
    assert all(row["uses_static_observable_thresholds"] is False for row in rows)
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["protein_folding_solved"] is False for row in rows)
