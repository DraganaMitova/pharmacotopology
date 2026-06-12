from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V80_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V80"


ALLOWED_CLASSIFICATIONS = {
    "merge_into_existing_word",
    "proto_grammar",
    "keep_as_unknown_clean_abstain",
    "reject_as_noise",
    "retire_due_to_regression",
}


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v80_certificate_passes_candidate_word_lexicon_delta_triage() -> None:
    cert = _read(V80_ROOT / "v80_candidate_word_lexicon_delta_triage_77_certificate.json")

    assert cert["status"] == "V80_CANDIDATE_WORD_LEXICON_DELTA_TRIAGE_PASSED"
    assert cert["engine_version_used"] == "E74"
    assert cert["baseline_engine_version"] == "E73"
    assert cert["input_candidate_word_count"] == 77
    assert cert["classified_candidate_word_count"] == 77
    assert sum(cert["classification_counts"].values()) == 77
    assert cert["failed_accepted_count"] == 0
    assert cert["candidate_words_classified_reproducibly"] is True
    assert cert["classification_hash"]
    assert cert["candidate_word_graph_node_count"] == 77
    assert cert["candidate_word_graph_edge_count"] > 0
    assert cert["paired_control_rows_total"] == 539
    assert cert["proto_grammars_matched_control_dominance_passed"] is True
    assert cert["merged_words_do_not_reduce_existing_grammar_performance"] is True
    assert cert["rejected_or_retired_words_record_pressure_reason"] is True
    assert cert["sentinel_regressions"] == 0
    assert cert["withheld_context_leakage_detected"] is False
    assert cert["coordinate_native_truth_stays_sealed"] is True
    assert cert["no_forced_expected_labels_used"] is True
    assert cert["no_static_thresholds_used"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["failed_controls"] == []


def test_v80_rows_record_pressure_ledger_and_matched_controls() -> None:
    rows = _read(V80_ROOT / "v80_candidate_word_lexicon_delta_triage_report.json")["rows"]

    assert len(rows) == 77
    assert {row["classification"] for row in rows}.issubset(ALLOWED_CLASSIFICATIONS)
    assert all(row["pressure_ledger"] for row in rows)
    assert all(row["compression_gain"] for row in rows)
    assert all(row["matched_control_dominance"]["uses_static_observable_thresholds"] is False for row in rows)
    assert all(
        len(row["matched_control_dominance"]["control_rows"]) == 7
        for row in rows
    )
    assert all(
        row["matched_control_dominance"]["matched_control_dominance_passed"] is True
        for row in rows
        if row["classification"] == "proto_grammar"
    )
    assert all(
        row["merge_candidate"]
        for row in rows
        if row["classification"] == "merge_into_existing_word"
    )
    assert all(
        row["rejected_or_retired_reason"]
        for row in rows
        if row["classification"] in {"reject_as_noise", "retire_due_to_regression"}
    )
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["folding_problem_solved"] is False for row in rows)
