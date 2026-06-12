from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    E73_WORD_LIFECYCLE,
    NEGATIVE_EVIDENCE_PRESSURE_CHANNELS,
    build_sealed_operator_state_packet,
    language_acquisition_observation_from_packet,
    protein_language_acquisition_cortex,
)


def _source(statement: str) -> dict[str, object]:
    return {
        "source_id": "E73_TEST_SOURCE",
        "source_class": "pure_non_coordinate",
        "source_role": "prediction_input",
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "spatial_proxy": False,
        "metadata_context_marks": [],
        "evidence_statement": statement,
    }


def test_e73_packet_exposes_epistemological_status_and_operator_state_packet() -> None:
    packet = build_sealed_operator_state_packet(
        target_id="E73_EPISTEMOLOGY_TEST",
        target_name="E73 epistemology test",
        sequence="ACDEFGHIKLMNPQRSTVWY" * 6,
        sources=[_source("raw blind protein context without a learned mechanism")],
        perturbations=[],
    )

    status = packet["epistemological_status"]
    assert status["language_layer"] == "protein_esperanto_mechanism_language"
    assert status["physical_execution_performed"] is False
    assert status["atomistic_md_performed"] is False
    assert status["operator_state_api"] == "propagate_operator_state"
    assert status["physical_basis_claim_allowed"] is False
    assert status["folding_problem_solved"] is False
    assert packet["kind"] == "V52_COARSE_OPERATOR_STATE_PROPAGATION_PACKET_v0"
    assert packet["operator_state_propagation_summary"]["kind"] == "PROTEIN_ESPERANTO_OPERATOR_STATE_PROPAGATION_v0"
    assert packet["hypothesized_interaction_language_map"] == packet["predicted_contact_interaction_probability_map"]


def test_e73_language_acquisition_cortex_builds_proto_words_from_pressure() -> None:
    packets = [
        build_sealed_operator_state_packet(
            target_id=f"E73_PRESSURE_TEST_{index}",
            target_name="E73 pressure test",
            sequence="GGGSSSQPNY" * 12,
            sources=[_source("unseen extracellular enzyme context not membrane no cofactor")],
            perturbations=[],
        )
        for index in range(2)
    ]
    observations = [
        language_acquisition_observation_from_packet(
            packet,
            visible_context_text="unseen extracellular enzyme context not membrane no cofactor",
            matched_control_dominance_passed=True,
        )
        for packet in packets
    ]
    cortex = protein_language_acquisition_cortex(observations)

    assert E73_WORD_LIFECYCLE == [
        "unseen_pattern",
        "pressure_cluster",
        "candidate_word",
        "proto_grammar",
        "learned_grammar",
        "rejected_or_merged_word",
    ]
    assert "not_membrane_pressure" in NEGATIVE_EVIDENCE_PRESSURE_CHANNELS
    assert cortex["engine_revision"] == "E73"
    assert cortex["candidate_word_count"] >= 1
    assert cortex["candidate_words_ranked_by_endogenous_pressure_support"] is True
    assert cortex["candidate_words"][0]["lifecycle_state"] in {"proto_grammar", "learned_grammar"}
    assert cortex["candidate_words"][0]["promotion_outcome"] in {
        "cleanly_abstained_pressure_support_not_yet_learned",
        "promoted_through_proto_grammar_tests",
    }
    assert cortex["physical_basis_claim_allowed"] is False
    assert cortex["folding_problem_solved"] is False


def test_e73_certificate_records_language_acquisition_layer() -> None:
    cert = json.loads(
        (
            ROOT
            / "data"
            / "protein_esperanto_engine"
            / "E73"
            / "e73_protein_language_acquisition_cortex_certificate.json"
        ).read_text(encoding="utf-8")
    )

    assert cert["engine_revision"] == "E73"
    assert cert["known_missing_word_queue_status"] == "closed_before_v79"
    assert cert["candidate_grammars_remaining"] == []
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
