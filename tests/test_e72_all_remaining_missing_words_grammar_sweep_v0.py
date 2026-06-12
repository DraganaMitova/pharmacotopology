from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    MECHANISM_CLASSES,
    SELF_DECISION_CANDIDATE_GRAMMARS,
    STATE_VARIABLES,
    UNIVERSAL_OPERATORS,
    build_sealed_simulation_packet,
)


def _source(statement: str) -> dict[str, object]:
    return {
        "source_id": "E72_TEST_SOURCE",
        "source_class": "pure_non_coordinate",
        "source_role": "prediction_input",
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "spatial_proxy": False,
        "metadata_context_marks": statement.split(),
        "evidence_statement": statement,
    }


def _packet(statement: str, sequence: str) -> dict[str, object]:
    return build_sealed_simulation_packet(
        target_id="E72_ALL_WORDS_TEST",
        target_name="E72 all remaining words test",
        sequence=sequence,
        sources=[_source(statement)],
        perturbations=[],
    )


def test_e72_promotes_all_remaining_candidate_words() -> None:
    assert SELF_DECISION_CANDIDATE_GRAMMARS == {}
    for mechanism in [
        "coiled_coil_register_topology",
        "repeat_solenoid_topology",
        "knotted_topology",
    ]:
        assert mechanism in MECHANISM_CLASSES
    for operator in [
        "heptad_register_operator",
        "coiled_coil_interface_operator",
        "repeat_phase_operator",
        "solenoid_axis_operator",
        "threading_operator",
        "topological_closure_operator",
    ]:
        assert operator in UNIVERSAL_OPERATORS
    for state in [
        "heptad_register_context",
        "global_repeat_topology",
        "topological_closure_constraint",
    ]:
        assert state in STATE_VARIABLES


def test_e72_accepts_promoted_words_and_exposes_known_word_fields() -> None:
    cases = [
        (
            "coiled_coil_register_topology",
            "coiled_coil_register",
            "known_coiled_coil_word",
            "coiled_coil_register heptad_repeat register_alignment hydrophobic_repeat_phase oligomeric_coiled_coil_core",
            "LEKLAAL" * 24,
        ),
        (
            "repeat_solenoid_topology",
            "repeat_solenoid_topology",
            "known_repeat_solenoid_word",
            "repeat_solenoid_topology repeat_unit solenoid_axis curved_repeat_stack local_repeat_closure global_repeat_topology",
            "TPRAGLYVPGSTNQ" * 18,
        ),
        (
            "knotted_topology",
            "knotted_topology",
            "known_knotted_topology_word",
            "knotted_topology knot_core_context threading_loop_context slipknot topological_closure_constraint long_range_threading_dependency",
            "AVILGSPNTQYFWDEKRH" * 16,
        ),
    ]
    for expected, word, known_field, statement, sequence in cases:
        packet = _packet(statement, sequence)
        judge = packet["self_decision_judge"]
        assert packet["selected_mechanism_grammar"]["mechanism_class"] == expected
        assert judge["final_self_decision"] == "accepted"
        assert judge[known_field] == word
        assert judge["missing_word_candidate"] is None
        assert packet["acceptance_firewall"]["acceptance_decision"] == "accepted"


def test_e72_context_free_sequence_only_hydrophobicity_abstains() -> None:
    packet = _packet("raw sequence only no biological grammar context", "LLLLLAAAVVVIIILL" * 12)

    assert packet["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain"
    assert packet["acceptance_firewall"]["acceptance_decision"] == "abstain_recommended"
    assert packet["self_decision_judge"]["final_self_decision"] == "clean_abstain_low_internal_consensus"


def test_e72_certificate_records_saturation_but_blocks_folding_claim() -> None:
    cert = json.loads(
        (
            ROOT
            / "data"
            / "protein_esperanto_engine"
            / "E72"
            / "e72_all_remaining_missing_words_grammar_sweep_certificate.json"
        ).read_text(encoding="utf-8")
    )

    assert cert["engine_revision"] == "E72"
    assert cert["promoted_words"] == [
        "coiled_coil_register",
        "repeat_solenoid_topology",
        "knotted_topology",
    ]
    assert cert["candidate_grammars_remaining"] == []
    assert cert["protein_esperanto_language_saturation_status"] == "current_known_missing_words_implemented"
    assert cert["protein_folding_solved"] is False
    assert cert["claim_allowed"] is False
    assert cert["next_required_batch"] == "V79_BLIND_RCSB_DISCOVERY_500_NO_KNOWN_MISSING_WORD_QUEUE"
