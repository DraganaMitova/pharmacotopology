from __future__ import annotations

import json
from pathlib import Path

from pharmacotopology.protein_esperanto_engine import (
    E75_PROTO_GRAMMAR_CLASSIFICATIONS,
    proto_grammar_crystallization_cortex,
)


ROOT = Path(__file__).resolve().parents[1]
V80_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V80"
V80P_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V80P"
E75_ROOT = ROOT / "data" / "protein_esperanto_engine" / "E75"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_e75_crystallization_cortex_splits_proto_grammars_without_inflation() -> None:
    triage = _read(V80_ROOT / "v80_candidate_word_lexicon_delta_triage_report.json")
    physical_rows = _read(V80P_ROOT / "v80p_independent_physical_holdout_gate_256_rows.json")["rows"]
    cortex = proto_grammar_crystallization_cortex(
        v80_triage_report=triage,
        v80p_holdout_rows=physical_rows,
        sentinel_regression_count=0,
    )

    classifications = {row["classification"] for row in cortex["crystallization_rows"]}

    assert cortex["engine_revision"] == "E75"
    assert cortex["baseline_engine_revision"] == "E74"
    assert cortex["input_proto_grammar_count"] == 45
    assert cortex["input_merged_word_count"] == 32
    assert len(cortex["crystallization_rows"]) == 45
    assert classifications.issubset(set(E75_PROTO_GRAMMAR_CLASSIFICATIONS))
    assert cortex["classification_counts"] == {
        "crystallized_grammar": 8,
        "keep_as_proto_unknown": 21,
        "merge_into_existing_word": 15,
        "retire_as_context_artifact": 1,
    }
    assert cortex["classification_counts"]["crystallized_grammar"] < 45
    assert cortex["crystallization_hash"]
    assert cortex["proto_grammars_classified_reproducibly"] is True
    assert cortex["sentinel_regressions"] == 0
    assert cortex["no_static_thresholds_used"] is True
    assert cortex["no_forced_expected_labels_used"] is True
    assert cortex["v79_v80_context_used_as_seed_not_proof"] is True
    assert cortex["physical_basis_claim_allowed"] is False
    assert cortex["folding_problem_solved"] is False
    assert all(row["compression_gain"]["uses_static_threshold"] is False for row in cortex["crystallization_rows"])


def test_e75_certificate_records_crystallization_without_folding_claim() -> None:
    cert = _read(E75_ROOT / "e75_proto_grammar_crystallization_cortex_certificate.json")

    assert cert["kind"] == "E75_PROTO_GRAMMAR_CRYSTALLIZATION_CORTEX_CERTIFICATE_v0"
    assert cert["engine_revision"] == "E75"
    assert cert["baseline_engine_revision"] == "E74"
    assert cert["input_proto_grammar_count"] == 45
    assert cert["input_merged_word_count"] == 32
    assert cert["classification_counts"] == {
        "crystallized_grammar": 8,
        "keep_as_proto_unknown": 21,
        "merge_into_existing_word": 15,
        "retire_as_context_artifact": 1,
    }
    assert cert["crystallization_hash"]
    assert cert["sentinel_regressions"] == 0
    assert cert["no_static_thresholds_used"] is True
    assert cert["no_forced_expected_labels_used"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
