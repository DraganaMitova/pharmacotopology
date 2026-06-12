from __future__ import annotations

import pharmacotopology.protein_esperanto_engine as engine
from pharmacotopology.protein_esperanto_engine import (
    KNOTTED_TOPOLOGY_LANGUAGE_CLAIM_SCOPE,
    KNOTTED_TOPOLOGY_STRUCTURAL_PROOF_REQUIREMENT,
    build_sealed_operator_state_packet,
    make_openmm_bridge_spec,
)


def _source(statement: str) -> dict[str, object]:
    return {
        "source_id": "KNOT_BOUNDARY_TEST_SOURCE",
        "source_class": "pure_non_coordinate",
        "source_role": "prediction_input",
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "spatial_proxy": False,
        "metadata_context_marks": [],
        "evidence_statement": statement,
    }


def test_engine_uses_current_operator_state_api_names_only() -> None:
    assert hasattr(engine, "propagate_operator_state")
    assert hasattr(engine, "build_sealed_operator_state_packet")
    assert not hasattr(engine, "simulate_operator_" + "trajectory")
    assert not hasattr(engine, "build_sealed_" + "simulation_packet")


def test_knot_grammar_is_language_scoped_and_structural_claim_blocked() -> None:
    packet = build_sealed_operator_state_packet(
        target_id="KNOT_LANGUAGE_SCOPE_TEST",
        target_name="Knot language scope test",
        sequence="ACDEFGHIKLMNPQRSTVWY" * 12,
        sources=[
            _source(
                "knotted_topology knot_core_context threading_loop_context "
                "slipknot_intermediate_context topological_closure_constraint "
                "long_range_threading_dependency"
            )
        ],
        perturbations=[],
    )

    mechanism = packet["selected_mechanism_grammar"]
    operator_state = packet["operator_state_propagation_summary"]
    judge = packet["self_decision_judge"]

    assert packet["kind"] == "V52_COARSE_OPERATOR_STATE_PROPAGATION_PACKET_v0"
    assert "trajectory_" + "summary" not in packet
    assert mechanism["mechanism_class"] == "knotted_topology"
    assert mechanism["selected_knotted_topology_claim_scope"] == KNOTTED_TOPOLOGY_LANGUAGE_CLAIM_SCOPE
    assert mechanism["structural_knot_detection_performed"] is False
    assert mechanism["structural_knot_claim_requires"] == KNOTTED_TOPOLOGY_STRUCTURAL_PROOF_REQUIREMENT
    assert operator_state["kind"] == "PROTEIN_ESPERANTO_OPERATOR_STATE_PROPAGATION_v0"
    assert operator_state["knot_claim_scope"] == KNOTTED_TOPOLOGY_LANGUAGE_CLAIM_SCOPE
    assert operator_state["structural_knot_detection_performed"] is False
    assert judge["known_knotted_topology_claim_scope"] == KNOTTED_TOPOLOGY_LANGUAGE_CLAIM_SCOPE
    assert judge["structural_knot_detection_performed"] is False
    assert judge["structural_knot_claim_requires"] == KNOTTED_TOPOLOGY_STRUCTURAL_PROOF_REQUIREMENT
    assert packet["acceptance_firewall"]["kind"] == "PROTEIN_ESPERANTO_SELF_DECISION_ACCEPTANCE_VIEW_v0"


def test_openmm_bridge_spec_cannot_unlock_physical_claim_by_itself() -> None:
    spec = make_openmm_bridge_spec()

    assert spec["kind"] == "V56_OPERATOR_TO_CUSTOM_FORCE_BRIDGE_SPEC_v0"
    assert spec["execution_status"] == "specification_only_no_openmm_execution"
    assert spec["openmm_executed"] is False
    assert spec["validated_coarse_execution_executed"] is False
    assert spec["physical_claim_allowed_from_spec"] is False
    assert spec["claim_ceiling"] == "language_or_proxy_only_until_real_execution"
