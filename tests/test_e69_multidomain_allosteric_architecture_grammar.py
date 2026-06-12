from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.protein_esperanto_engine import MECHANISM_CLASSES, STATE_VARIABLES, build_sealed_operator_state_packet


def _source(statement: str, *, withheld: list[str] | None = None) -> dict[str, object]:
    return {
        "source_id": "E69_TEST_SOURCE",
        "source_class": "pure_non_coordinate",
        "source_role": "prediction_input",
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "spatial_proxy": False,
        "metadata_context_marks": statement.split(),
        "withheld_context_marks": withheld or [],
        "evidence_statement": statement,
    }


def _packet(
    statement: str,
    *,
    withheld: list[str] | None = None,
    perturbations: list[dict[str, object]] | None = None,
    sequence: str | None = None,
) -> dict[str, object]:
    sequence = sequence or ("VIFYWTERKLGHADP" * 18)[:270]
    return build_sealed_operator_state_packet(
        target_id="E69_MULTIDOMAIN_TEST",
        target_name="E69 multidomain allostery test",
        sequence=sequence,
        sources=[_source(statement, withheld=withheld)],
        perturbations=perturbations or [],
    )


def test_e69_adds_multidomain_mechanism_and_state_words() -> None:
    assert "multidomain_allosteric_architecture" in MECHANISM_CLASSES
    for word in [
        "multidomain_allostery",
        "domain_boundary",
        "hinge_region",
        "interdomain_lock",
        "allosteric_basin_shift",
        "domain_reorientation",
        "modular_architecture",
        "domain_swapping",
    ]:
        assert word in STATE_VARIABLES

    cert = json.loads(
        (
            ROOT
            / "data"
            / "protein_esperanto_engine"
            / "E69"
            / "e69_multidomain_allosteric_architecture_grammar_certificate.json"
        ).read_text(encoding="utf-8")
    )
    assert cert["engine_revision"] == "E69"
    assert cert["trigger_batch"] == "V74_RCSB_NONREDUNDANT_200_DISCOVERY_E68"
    assert cert["trigger_failure_count"] == 33
    assert cert["new_mechanism_class"] == "multidomain_allosteric_architecture"


def test_e69_selects_multidomain_allostery_and_outputs_readouts() -> None:
    packet = _packet("multidomain_allostery domain_boundary hinge_region interdomain_lock allosteric_basin_shift domain_reorientation modular_architecture")
    mechanism = packet["selected_mechanism_grammar"]
    final = packet["operator_state_propagation_summary"]["final_state_summary"]

    assert mechanism["mechanism_class"] == "multidomain_allosteric_architecture"
    assert mechanism["selected_multidomain_word"] == "allosteric_basin_shift"
    assert packet["self_decision_judge"]["zero_failed_accepted_required"] is True
    assert packet["self_decision_judge"]["final_self_decision"] == "accepted"
    assert packet["self_decision_judge"]["dominance_law"] == "single_dominant_learned_mechanism_bound_across_views"
    assert packet["self_decision_judge"]["cross_view_binding_probe"]["missing_view_families"] == []
    assert packet["self_decision_judge"]["operator_basis_stability"] in {
        "stable_under_endogenous_operator_basis_probe",
        "definition_sensitive_under_semantic_operator_basis_probe",
    }
    assert packet["self_decision_judge"]["coefficient_probe_mode"] == "endogenous_observed_operator_permutations_no_static_scale_range"
    assert packet["self_decision_judge"]["physical_basis_claim_allowed"] is False
    assert packet["acceptance_firewall"]["acceptance_decision"] == "accepted"
    assert final["multidomain_allostery"] > 0.0
    assert final["domain_boundary"] > 0.0
    assert final["hinge_region"] > 0.0
    assert final["interdomain_lock"] > 0.0
    assert final["allosteric_basin_shift"] > 0.0
    assert final["domain_reorientation"] > 0.0
    assert final["modular_architecture"] > 0.0
    assert any(row["interaction_type"] == "interdomain_allosteric_lock" for row in packet["predicted_contact_interaction_probability_map"])


def test_e69_selects_domain_swapping_subtype() -> None:
    packet = _packet("domain_swapping domain-swapped interdomain_lock modular_architecture")
    mechanism = packet["selected_mechanism_grammar"]
    final = packet["operator_state_propagation_summary"]["final_state_summary"]

    assert mechanism["mechanism_class"] == "multidomain_allosteric_architecture"
    assert mechanism["selected_multidomain_word"] == "domain_swapping"
    assert final["domain_swapping"] > 0.0


def test_e69_interdomain_lock_perturbation_decreases_lock() -> None:
    packet = _packet(
        "multidomain_allostery domain_boundary hinge_region interdomain_lock allosteric_basin_shift",
        perturbations=[
            {
                "perturbation_id": "E69_LOCK_DAMAGE",
                "description": "damage the interdomain lock",
                "operator_scales": {"interface_operator": 0.45},
                "lock_damage": 0.45,
                "metric": "interdomain_lock",
                "expected_direction": "decrease",
            }
        ],
    )
    row = packet["predicted_perturbation_table"][0]
    assert row["observed_direction"] == "decrease"
    assert row["direction_passed"] is True


def test_e69_complex_word_boundary_preserves_low_complexity_disorder() -> None:
    sequence = "MASASSSQRGRSGSGNFGGGRGGGFGGNDNFGRGGNFSGRGGFGGSRGGGGYGGSGDGYNGFGNDGGYGGGGPGY"
    packet = build_sealed_operator_state_packet(
        target_id="E69_LOW_COMPLEXITY_BOUNDARY_TEST",
        target_name="E69 low complexity boundary test",
        sequence=sequence,
        sources=[
            _source(
                "Low complexity prion-like LCD with aromatic glycine marks, phase separation tendency, and disorder ensemble behavior."
            )
        ],
        perturbations=[],
    )

    assert packet["selected_mechanism_grammar"]["mechanism_class"] == "intrinsic_disorder_phase_separation"
    assert packet["selected_mechanism_grammar"]["selection_reason"] == "low_complexity_disorder_phase_evidence"


def test_e69_former_missing_word_is_now_promoted_by_e72() -> None:
    packet = _packet(
        "coiled_coil_register heptad_repeat register_alignment hydrophobic_repeat_phase oligomeric_coiled_coil_core",
        sequence="LEKLAAL" * 24,
    )
    judge = packet["self_decision_judge"]
    acceptance_view = packet["acceptance_firewall"]

    assert judge["zero_failed_accepted_required"] is True
    assert packet["selected_mechanism_grammar"]["mechanism_class"] == "coiled_coil_register_topology"
    assert judge["final_self_decision"] == "accepted"
    assert judge["known_coiled_coil_word"] == "coiled_coil_register"
    assert judge["missing_word_candidate"] is None
    assert judge["wrong_grammar_separation"] == "wrong_grammars_fail"
    assert acceptance_view["unknown_word_signals"] == []
    assert acceptance_view["acceptance_decision"] == "accepted"


def test_e69_preserves_prior_repair_priorities_and_withheld_isolation() -> None:
    membrane = _packet("membrane_context_strong transmembrane_context topology_evidence multidomain_allostery")
    assembly = _packet("assembly_required_core assembly_required_folding partner_completed_core multidomain_allostery")
    metal = _packet("metal_cluster_geometry ligand_locked_basin multidomain_allostery")
    disorder = _packet("disorder_context fold-upon-binding region multidomain_allostery")
    beta = _packet("closed_beta_topology strand_register beta_sheet_closure multidomain_allostery")
    withheld = _packet(
        "multidomain_allostery domain_boundary hinge_region",
        withheld=["assembly_required_core", "membrane_beta_barrel", "metal_cluster_geometry"],
    )

    assert membrane["selected_mechanism_grammar"]["mechanism_class"] == "membrane_multidomain_folding_proteostasis"
    assert assembly["selected_mechanism_grammar"]["mechanism_class"] == "assembly_required_folding"
    assert metal["selected_mechanism_grammar"]["mechanism_class"] == "metal_cluster_and_ligand_locked_basin"
    assert disorder["selected_mechanism_grammar"]["mechanism_class"] == "disorder_boundary_and_fold_upon_binding"
    assert beta["selected_mechanism_grammar"]["mechanism_class"] == "beta_closure_topology"
    assert withheld["selected_mechanism_grammar"]["mechanism_class"] == "multidomain_allosteric_architecture"
