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
    PROCESS_CLASSES,
    SELF_DECISION_CANDIDATE_GRAMMARS,
    STATE_VARIABLES,
    UNIVERSAL_OPERATORS,
    build_sealed_operator_state_packet,
)


SIGNAL_SEQUENCE = "MKKLLLLLLLLLLLLLLLLAAASA" + ("STNQDEKRASTNQDEKR" * 12)


def _source(statement: str) -> dict[str, object]:
    return {
        "source_id": "E71_TEST_SOURCE",
        "source_class": "pure_non_coordinate",
        "source_role": "prediction_input",
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "spatial_proxy": False,
        "metadata_context_marks": statement.split(),
        "evidence_statement": statement,
    }


def _packet(statement: str, *, sequence: str = SIGNAL_SEQUENCE, perturbations: list[dict[str, object]] | None = None) -> dict[str, object]:
    return build_sealed_operator_state_packet(
        target_id="E71_SIGNAL_TM_TEST",
        target_name="E71 signal peptide / true TM test",
        sequence=sequence,
        sources=[_source(statement)],
        perturbations=perturbations or [],
    )


def test_e71_promotes_signal_tm_boundary_to_learned_grammar_space() -> None:
    assert "signal_peptide_vs_true_tm_routing" in MECHANISM_CLASSES
    assert "signal_peptide_tm_boundary" in PROCESS_CLASSES
    assert "signal_peptide_vs_true_TM" not in SELF_DECISION_CANDIDATE_GRAMMARS
    for operator in [
        "signal_peptide_routing_operator",
        "tm_insertion_operator",
        "cleavage_context_operator",
        "secretory_routing_operator",
        "membrane_pressure_operator",
        "frustration_operator",
    ]:
        assert operator in UNIVERSAL_OPERATORS
    for word in [
        "signal_peptide_routing_context",
        "cleavage_site_context",
        "n_terminal_secretory_hydrophobic_patch",
        "true_transmembrane_span_context",
        "single_pass_tm_conflict",
        "multi_pass_tm_conflict",
        "secretory_lumenal_routing",
        "membrane_insertion_routing",
        "signal_anchor_ambiguity",
    ]:
        assert word in STATE_VARIABLES

    cert = json.loads(
        (
            ROOT
            / "data"
            / "protein_esperanto_engine"
            / "E71"
            / "e71_signal_peptide_true_tm_routing_grammar_certificate.json"
        ).read_text(encoding="utf-8")
    )
    assert cert["engine_revision"] == "E71"
    assert cert["candidate_word_promoted_to_learned"] == "signal_peptide_vs_true_TM"
    assert cert["proof_batch"] == "V77_SIGNAL_TM_BOUNDARY_PANEL_300"
    assert cert["v77_failed_accepted_count"] == 0
    assert cert["v77_uses_static_observable_thresholds"] is False


def test_e71_accepts_signal_peptide_only_when_views_and_controls_bind() -> None:
    packet = _packet(
        "signal_peptide_vs_true_TM signal_peptide_routing_context cleavage_site_context "
        "secretory_lumenal_routing cleaved signal peptide n-terminal signal peptide",
        perturbations=[
            {
                "perturbation_id": "E71_N_TERMINAL_SIGNAL_MASK",
                "description": "mask N-terminal signal and cleavage context",
                "n_terminal_mask": 1.0,
                "cleavage_loss": 1.0,
                "metric": "signal_peptide_routing_context",
                "expected_direction": "decrease",
            },
            {
                "perturbation_id": "E71_TRUE_TM_CONTROL",
                "description": "paired true-TM boundary control",
                "true_tm_decoy": 1.0,
                "metric": "signal_peptide_routing_context",
                "expected_direction": "decrease",
            },
        ],
    )
    mechanism = packet["selected_mechanism_grammar"]
    judge = packet["self_decision_judge"]
    final = packet["operator_state_propagation_summary"]["final_state_summary"]

    assert mechanism["mechanism_class"] == "signal_peptide_vs_true_tm_routing"
    assert mechanism["selected_signal_peptide_word"] == "signal_peptide_vs_true_TM"
    assert judge["final_self_decision"] == "accepted"
    assert judge["known_signal_peptide_word"] == "signal_peptide_vs_true_TM"
    assert judge["coefficient_probe_mode"] == "endogenous_observed_operator_permutations_no_static_scale_range"
    assert judge["coefficient_perturbation_probe"]["coefficient_scale_values_used"] == []
    assert judge["cross_view_binding_probe"]["missing_view_families"] == []
    assert judge["temporal_binding"] == "selected_observables_temporally_coherent"
    assert packet["acceptance_firewall"]["acceptance_decision"] == "accepted"
    assert final["signal_peptide_routing_context"] > final["membrane_insertion_routing"]
    assert final["secretory_lumenal_routing"] > final["signal_anchor_ambiguity"]
    assert any(
        row["interaction_type"] == "signal_peptide_secretory_routing_contact"
        for row in packet["predicted_contact_interaction_probability_map"]
    )
    signal_delta = packet["predicted_perturbation_table"][0]["baseline_value"] - packet["predicted_perturbation_table"][0]["perturbed_value"]
    tm_delta = packet["predicted_perturbation_table"][1]["baseline_value"] - packet["predicted_perturbation_table"][1]["perturbed_value"]
    assert signal_delta > tm_delta


def test_e71_negative_gates_preserve_true_tm_disulfide_and_missing_words() -> None:
    true_tm = _packet(
        "membrane_context_strong transmembrane_context topology_evidence true_transmembrane_span_context",
        sequence=("MSTNQDEKR" * 5) + "LLLLLLLLLLLLLLLLLLLL" + ("STNQDEKR" * 10),
    )
    anchor = _packet("signal_peptide_vs_true_TM signal_anchor_ambiguity signal anchor uncleaved signal anchor")
    disulfide = _packet(
        "disulfide_secretory_redox_context disulfide_bond_topology secretory_redox_context "
        "extracellular_stabilized_fold glycosylation_context cysteine_pairing_constraint",
        sequence=("ACDECGHIKCDEFGCNPQRSTCAVCDEF" * 10)[:280],
    )
    coiled = _packet(
        "coiled_coil_register heptad_repeat register_alignment hydrophobic_repeat_phase oligomeric_coiled_coil_core",
        sequence="LEKLAAL" * 24,
    )

    assert true_tm["selected_mechanism_grammar"]["mechanism_class"] == "membrane_multidomain_folding_proteostasis"
    assert true_tm["self_decision_judge"]["final_self_decision"] == "accepted"

    assert anchor["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain"
    assert anchor["acceptance_firewall"]["acceptance_decision"] == "abstain_recommended"

    assert disulfide["selected_mechanism_grammar"]["mechanism_class"] == "secretory_disulfide_redox_topology"
    assert disulfide["self_decision_judge"]["final_self_decision"] == "accepted"

    assert coiled["selected_mechanism_grammar"]["mechanism_class"] == "coiled_coil_register_topology"
    assert coiled["self_decision_judge"]["final_self_decision"] == "accepted"
    assert coiled["self_decision_judge"]["known_coiled_coil_word"] == "coiled_coil_register"
    assert coiled["self_decision_judge"]["missing_word_candidate"] is None
