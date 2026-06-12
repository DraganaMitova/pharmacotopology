from __future__ import annotations

import json
from pathlib import Path

from pharmacotopology.protein_esperanto_engine import (
    E74_CANDIDATE_WORD_CLASSIFICATIONS,
    candidate_word_compression_cortex,
)


ROOT = Path(__file__).resolve().parents[1]
V79_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V79"
V79P_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V79P"
E74_ROOT = ROOT / "data" / "protein_esperanto_engine" / "E74"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_e74_compression_cortex_classifies_all_v79_candidate_words() -> None:
    lexicon = _read(V79_ROOT / "v79_lexicon_delta_report.json")
    scoring = _read(V79_ROOT / "v79_blind_language_discovery_scoring_report.json")["rows"]
    physical = _read(V79P_ROOT / "v79p_physical_falsification_128_rows.json")["rows"]
    sentinel = _read(V79_ROOT / "v79_old_grammar_sentinel_replay.json")

    cortex = candidate_word_compression_cortex(
        v79_lexicon_delta=lexicon,
        v79_scoring_rows=scoring,
        v79p_physical_rows=physical,
        sentinel_regression_count=sentinel["sentinel_regressions"],
    )

    classifications = {row["classification"] for row in cortex["compression_rows"]}

    assert cortex["engine_revision"] == "E74"
    assert cortex["baseline_engine_revision"] == "E73"
    assert cortex["input_candidate_word_count"] == 77
    assert cortex["classified_candidate_word_count"] == 77
    assert classifications.issubset(set(E74_CANDIDATE_WORD_CLASSIFICATIONS))
    assert cortex["candidate_word_graph"]["node_count"] == 77
    assert cortex["candidate_word_graph"]["edge_count"] > 0
    assert cortex["compression_hash"]
    assert cortex["candidate_words_classified_reproducibly"] is True
    assert cortex["failed_accepted_count"] == 0
    assert cortex["sentinel_regressions"] == 0
    assert cortex["no_static_thresholds_used"] is True
    assert cortex["physical_basis_claim_allowed"] is False
    assert cortex["folding_problem_solved"] is False
    assert all(row["compression_gain"]["uses_static_threshold"] is False for row in cortex["compression_rows"])


def test_e74_certificate_records_compression_without_physical_claim() -> None:
    cert = _read(E74_ROOT / "e74_candidate_word_compression_cortex_certificate.json")

    assert cert["kind"] == "E74_CANDIDATE_WORD_COMPRESSION_CORTEX_CERTIFICATE_v0"
    assert cert["engine_revision"] == "E74"
    assert cert["baseline_engine_revision"] == "E73"
    assert cert["input_candidate_word_count"] == 77
    assert sum(cert["classification_counts"].values()) == 77
    assert cert["classification_hash"]
    assert cert["candidate_word_graph_node_count"] == 77
    assert cert["candidate_word_graph_edge_count"] > 0
    assert cert["no_static_thresholds_used"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
