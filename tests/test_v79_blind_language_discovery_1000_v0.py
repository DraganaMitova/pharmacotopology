from __future__ import annotations

import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
V79_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V79"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v79_certificate_passes_blind_language_discovery() -> None:
    cert = _read(V79_ROOT / "v79_blind_language_discovery_1000_certificate.json")

    assert cert["status"] == "V79_BLIND_LANGUAGE_DISCOVERY_PASSED"
    assert cert["engine_version_used"] == "E73"
    assert cert["baseline_engine_version"] == "E72"
    assert cert["targets_total"] == 1000
    assert cert["failed_accepted_count"] == 0
    assert cert["accepted_accuracy"] == pytest.approx(1.0)
    assert cert["no_known_missing_word_queue_used"] is True
    assert cert["no_target_specific_expected_mechanism_labels_used_for_prediction"] is True
    assert cert["candidate_word_count"] >= 1
    assert cert["candidate_words_ranked_by_endogenous_pressure_support"] is True
    assert cert["withheld_context_leakage_detected"] is False
    assert cert["sentinel_regressions"] == 0
    assert cert["controls_passed"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False


def test_v79_rows_use_matched_controls_and_pressure_observations() -> None:
    scoring = _read(V79_ROOT / "v79_blind_language_discovery_scoring_report.json")
    rows = scoring["rows"]
    accepted = [row for row in rows if row["acceptance_decision"] == "accepted"]

    assert len(rows) == 1000
    assert all(row["matched_control_dominance"]["uses_static_observable_thresholds"] is False for row in rows)
    assert all(row["matched_control_dominance_passed"] for row in accepted)
    assert all(row["failed_accepted"] is False for row in accepted)
    assert any(row["active_pressure_channels"] for row in rows)
    assert all(row["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["folding_problem_solved"] is False for row in rows)


def test_v79_lexicon_delta_is_reproducible_and_nonpromissory() -> None:
    cert = _read(V79_ROOT / "v79_blind_language_discovery_1000_certificate.json")
    lexicon = _read(V79_ROOT / "v79_lexicon_delta_report.json")

    assert lexicon["candidate_word_proposal_hash"] == cert["candidate_word_proposal_hash"]
    assert lexicon["candidate_word_count"] == cert["candidate_word_count"]
    assert all(
        row["promotion_outcome"].startswith("cleanly_abstained")
        or row["promotion_outcome"] == "promoted_through_proto_grammar_tests"
        for row in lexicon["candidate_words"]
    )
    assert all(row["physical_execution_expectations"] for row in lexicon["candidate_words"])
    assert lexicon["physical_basis_claim_allowed"] is False
    assert lexicon["folding_problem_solved"] is False
