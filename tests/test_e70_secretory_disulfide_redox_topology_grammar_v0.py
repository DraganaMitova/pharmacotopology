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
    build_sealed_operator_state_packet,
)


def _source(statement: str, *, withheld: list[str] | None = None) -> dict[str, object]:
    return {
        "source_id": "E70_TEST_SOURCE",
        "source_class": "pure_non_coordinate",
        "source_role": "prediction_input",
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "spatial_proxy": False,
        "metadata_context_marks": statement.split(),
        "withheld_context_marks": withheld or [],
        "evidence_statement": statement,
    }


def _packet(statement: str, *, sequence: str | None = None, perturbations: list[dict[str, object]] | None = None) -> dict[str, object]:
    return build_sealed_operator_state_packet(
        target_id="E70_SECRETORY_DISULFIDE_TEST",
        target_name="E70 secretory disulfide test",
        sequence=sequence or ("ACDECGHIKCDEFGCNPQRSTCAVCDEF" * 10)[:280],
        sources=[_source(statement)],
        perturbations=perturbations or [],
    )


def test_e70_promotes_disulfide_word_to_learned_grammar_space() -> None:
    assert "secretory_disulfide_redox_topology" in MECHANISM_CLASSES
    assert "disulfide_secretory_redox_context" not in SELF_DECISION_CANDIDATE_GRAMMARS
    assert "disulfide_pairing_operator" in UNIVERSAL_OPERATORS
    assert "secretory_redox_operator" in UNIVERSAL_OPERATORS
    for word in [
        "secretory_redox_context",
        "disulfide_pairing_topology",
        "cysteine_pairing_constraint",
        "extracellular_stabilized_fold",
        "glycosylation_context",
        "redox_mispaired_frustration",
        "signal_peptide_removed_context",
        "secretory_quality_control",
    ]:
        assert word in STATE_VARIABLES

    cert = json.loads(
        (
            ROOT
            / "data"
            / "protein_esperanto_engine"
            / "E70"
            / "e70_secretory_disulfide_redox_topology_grammar_certificate.json"
        ).read_text(encoding="utf-8")
    )
    assert cert["engine_revision"] == "E70"
    assert cert["candidate_word_promoted_to_learned"] == "disulfide_secretory_redox_context"
    assert cert["new_mechanism_class"] == "secretory_disulfide_redox_topology"
    assert cert["proof_batch"] == "V76_SECRETORY_DISULFIDE_REPAIR_PANEL_200"
    assert cert["v76_failed_accepted_count"] == 0
    assert cert["v76_fixed_accepted_support_thresholds_used"] is False


def test_e70_accepts_secretory_disulfide_only_when_views_bind() -> None:
    packet = _packet(
        "disulfide_secretory_redox_context disulfide_bond_topology secretory_redox_context "
        "extracellular_stabilized_fold glycosylation_context cysteine_pairing_constraint",
        perturbations=[
            {
                "perturbation_id": "E70_DISULFIDE_DAMAGE",
                "description": "reduce or mispair secretory disulfide topology",
                "operator_scales": {
                    "disulfide_pairing_operator": 0.24,
                    "secretory_redox_operator": 0.50,
                },
                "disulfide_damage": 0.58,
                "redox_shift": 0.25,
                "metric": "disulfide_pairing_topology",
                "expected_direction": "decrease",
            }
        ],
    )
    mechanism = packet["selected_mechanism_grammar"]
    judge = packet["self_decision_judge"]
    final = packet["operator_state_propagation_summary"]["final_state_summary"]

    assert mechanism["mechanism_class"] == "secretory_disulfide_redox_topology"
    assert mechanism["selected_secretory_disulfide_word"] == "disulfide_secretory_redox_context"
    assert judge["final_self_decision"] == "accepted"
    assert judge["known_secretory_disulfide_word"] == "disulfide_secretory_redox_context"
    assert judge["coefficient_probe_mode"] == "endogenous_observed_operator_permutations_no_static_scale_range"
    assert judge["coefficient_perturbation_probe"]["coefficient_scale_values_used"] == []
    assert judge["cross_view_binding_probe"]["missing_view_families"] == []
    assert judge["temporal_binding"] == "selected_observables_temporally_coherent"
    assert packet["acceptance_firewall"]["acceptance_decision"] == "accepted"
    assert final["disulfide_pairing_topology"] > final["redox_mispaired_frustration"]
    assert final["secretory_redox_context"] > 0.0
    assert final["extracellular_stabilized_fold"] > 0.0
    assert any(
        row["interaction_type"] == "secretory_disulfide_pairing_contact"
        for row in packet["predicted_contact_interaction_probability_map"]
    )
    assert packet["predicted_perturbation_table"][0]["perturbed_value"] < packet["predicted_perturbation_table"][0]["baseline_value"]


def test_e70_negative_gates_prevent_disulfide_overcapture() -> None:
    metal = _packet("metal_cluster_geometry coordination_shell_integrity cys-his coordination zinc-binding disulfide_secretory_redox_context")
    signal_only = _packet("signal peptide without disulfide signal_peptide_vs_true_TM")
    coiled = _packet("coiled_coil_register heptad_repeat register_alignment disulfide_secretory_redox_context")
    membrane = _packet(
        "membrane_context_strong transmembrane_context topology_evidence explicit_tm_evidence disulfide_secretory_redox_context",
        sequence=("M" + "LVVVVVVVVVVVVVVVVVVVVV" + "ACDECGHIKCDEFGCNPQRSTCAVCDEF" * 8)[:260],
    )

    assert metal["selected_mechanism_grammar"]["mechanism_class"] == "metal_cluster_and_ligand_locked_basin"
    assert metal["self_decision_judge"]["final_self_decision"] == "accepted"

    assert signal_only["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain"
    assert signal_only["self_decision_judge"]["final_self_decision"] == "clean_abstain_low_internal_consensus"
    assert signal_only["self_decision_judge"]["missing_word_candidate"] is None
    assert signal_only["acceptance_firewall"]["acceptance_decision"] == "abstain_recommended"

    assert coiled["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain"
    assert coiled["self_decision_judge"]["final_self_decision"] in {
        "clean_abstain_low_internal_consensus",
        "clean_abstain_conflict",
    }
    assert coiled["self_decision_judge"]["missing_word_candidate"] is None
    assert coiled["acceptance_firewall"]["acceptance_decision"] == "abstain_recommended"

    assert membrane["selected_mechanism_grammar"]["mechanism_class"] == "membrane_multidomain_folding_proteostasis"
    assert membrane["self_decision_judge"]["missing_word_candidate"] != "disulfide_secretory_redox_context"
