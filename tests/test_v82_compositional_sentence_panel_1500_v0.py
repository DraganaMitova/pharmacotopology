from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V82_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V82"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v82_certificate_passes_compositional_sentence_panel() -> None:
    cert = _read(V82_ROOT / "v82_compositional_sentence_panel_1500_certificate.json")

    assert cert["status"] == "V82_COMPOSITIONAL_SENTENCE_PANEL_PASSED"
    assert cert["engine_version_used"] == "E76"
    assert cert["baseline_engine_version"] == "E75"
    assert cert["targets_total"] == 1500
    assert cert["row_family_counts"] == {
        "bag_of_words_control": 150,
        "enemy_order_control": 200,
        "four_word_sentence": 200,
        "metadata_sequence_masked_hard_abstain_control": 150,
        "old_single_word_sentinel": 200,
        "three_word_sentence": 300,
        "two_word_sentence": 300,
    }
    assert cert["sentence_supported_count"] == 800
    assert cert["single_word_supported_count"] == 200
    assert cert["composition_required_sentence_supported_count"] == 800
    assert cert["sentence_supported_count"] > cert["single_word_supported_count"]
    assert cert["old_single_word_sentinels_preserved"] == 200
    assert cert["bag_of_words_control_beats_bound_sentence_count"] == 0
    assert cert["wrong_word_order_controls_failed"] == 200
    assert cert["masked_hard_abstain_controls_supported"] == 150
    assert cert["failed_accepted_count"] == 0
    assert cert["sentinel_regressions"] == 0
    assert cert["withheld_context_leakage_detected"] is False
    assert cert["coordinate_native_truth_stays_sealed"] is True
    assert cert["static_thresholds_used"] is False
    assert cert["no_static_thresholds_used"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["failed_controls"] == []


def test_v82_rows_preserve_sentences_sentinels_and_controls() -> None:
    rows = _read(V82_ROOT / "v82_compositional_sentence_panel_report.json")["rows"]

    assert len(rows) == 1500
    sentence_rows = [row for row in rows if row["row_family"] in {"two_word_sentence", "three_word_sentence", "four_word_sentence"}]
    sentinel_rows = [row for row in rows if row["row_family"] == "old_single_word_sentinel"]
    bag_rows = [row for row in rows if row["row_family"] == "bag_of_words_control"]
    enemy_rows = [row for row in rows if row["row_family"] == "enemy_order_control"]
    masked_rows = [row for row in rows if row["row_family"] == "metadata_sequence_masked_hard_abstain_control"]

    assert len(sentence_rows) == 800
    assert len(sentinel_rows) == 200
    assert all(row["sentence_supported"] for row in sentence_rows)
    assert all(row["single_word_supported"] for row in sentinel_rows)
    assert all(not row.get("bag_of_words_control_beats_bound_sentence", False) for row in bag_rows)
    assert all(row["wrong_word_order_fails"] for row in enemy_rows)
    assert all(row["clean_abstain_supported"] for row in masked_rows)
    assert all(row["failed_accepted"] is False for row in rows)
    assert all(row["sentinel_regression"] is False for row in rows)
    assert all(row["coordinate_native_truth_used"] is False for row in rows)
    assert all(row["uses_static_observable_thresholds"] is False for row in rows)
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["protein_folding_solved"] is False for row in rows)
