from __future__ import annotations

import json
from pathlib import Path

from pharmacotopology.protein_esperanto_engine import (
    E76_COMPOSITION_LAWS,
    MECHANISM_CLASSES,
    compositional_protein_sentence_grammar,
    protein_sentence_packet,
)


ROOT = Path(__file__).resolve().parents[1]
E76_ROOT = ROOT / "data" / "protein_esperanto_engine" / "E76"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_e76_sentence_packet_binds_words_instead_of_accepting_bag_of_words() -> None:
    packet = protein_sentence_packet(
        words=[
            "membrane_multidomain_folding_proteostasis",
            "coiled_coil_register_topology",
        ],
        known_words=set(MECHANISM_CLASSES),
        sentence_id="E76_TEST_SENTENCE",
    )

    assert packet["kind"] == "E76_PROTEIN_SENTENCE_PACKET_v0"
    assert packet["head_word"] == "coiled_coil_register_topology"
    assert packet["modifier_words"] == ["membrane_multidomain_folding_proteostasis"]
    assert packet["word_order"] == [
        "membrane_multidomain_folding_proteostasis",
        "coiled_coil_register_topology",
    ]
    assert packet["dependency_edges"] == [
        {
            "source": "membrane_multidomain_folding_proteostasis",
            "target": "coiled_coil_register_topology",
            "relation": "routing_clause_constrains_topology_clause",
            "source_clause": "routing_clause",
            "target_clause": "topology_clause",
        }
    ]
    assert packet["controls"]["selected_sentence_beats_bag_of_words"] is True
    assert packet["controls"]["selected_sentence_beats_wrong_order"] is True
    assert packet["controls"]["selected_sentence_beats_wrong_head"] is True
    assert packet["controls"]["modifier_cannot_steal_head"] is True
    assert packet["controls"]["uses_static_threshold"] is False
    assert packet["sentence_acceptance_decision"] == "accepted_supported"
    assert packet["physical_basis_claim_allowed"] is False
    assert packet["protein_folding_solved"] is False


def test_e76_proto_unknown_words_stay_cleanly_abstainable() -> None:
    packet = protein_sentence_packet(
        words=["not_globular_4250a3073b24", "not_globular_479cb9dc15ce"],
        known_words=set(MECHANISM_CLASSES),
        proto_unknown_words={"not_globular_4250a3073b24", "not_globular_479cb9dc15ce"},
        sentence_id="E76_TEST_PROTO_UNKNOWN",
    )

    assert packet["sentence_acceptance_decision"] == "clean_abstain_supported"
    assert packet["dominant_phrase"] is None
    assert packet["abstained_phrases"]
    assert all(row["status"] == "proto_unknown_clean_abstain" for row in packet["word_entries"])


def test_e76_compositional_grammar_and_certificate_block_physical_claim() -> None:
    grammar = compositional_protein_sentence_grammar(
        sentence_word_lists=[
            ["signal_peptide_vs_true_tm_routing", "secretory_disulfide_redox_topology"],
            ["assembly_required_folding", "oligomerization_controlled_folding"],
        ],
        known_words=set(MECHANISM_CLASSES),
    )
    cert = _read(E76_ROOT / "e76_compositional_protein_sentence_grammar_certificate.json")

    assert grammar["engine_revision"] == "E76"
    assert grammar["composition_laws"] == E76_COMPOSITION_LAWS
    assert grammar["accepted_sentence_count"] == 2
    assert grammar["no_static_thresholds_used"] is True
    assert grammar["physical_basis_claim_allowed"] is False
    assert grammar["protein_folding_solved"] is False
    assert cert["kind"] == "E76_COMPOSITIONAL_PROTEIN_SENTENCE_GRAMMAR_CERTIFICATE_v0"
    assert cert["engine_revision"] == "E76"
    assert cert["composition_laws"] == E76_COMPOSITION_LAWS
    assert cert["no_static_thresholds_used"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
