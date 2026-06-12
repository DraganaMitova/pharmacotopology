from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.protein_esperanto_engine import build_sealed_operator_state_packet


def test_e61_cofactor_context_activates_ligand_stabilization_grammar() -> None:
    packet = build_sealed_operator_state_packet(
        target_id="E61_MYOGLOBIN_CONTEXT_PROBE",
        target_name="myoglobin heme context probe",
        sequence="VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG",
        sources=[
            {
                "source_id": "E61_HEME_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "evidence_statement": "RCSB non-coordinate metadata exposes heme_context cofactor_context ligand_context for an oxygen-binding globin.",
            }
        ],
        perturbations=[
            {
                "perturbation_id": "E61_COFATOR_REMOVAL",
                "description": "remove cofactor stabilization",
                "cofactor_loss": 0.60,
                "metric": "basin:ligand_stabilized_basin",
                "expected_direction": "decrease",
            }
        ],
    )
    assert packet["selected_mechanism_grammar"]["mechanism_class"] == "cofactor_ligand_assisted_stabilization"
    final = packet["operator_state_propagation_summary"]["final_state_summary"]
    assert final["state_basin_occupancy"]["ligand_stabilized_basin"] > final["state_basin_occupancy"]["apo_weak_basin"]
    assert packet["predicted_perturbation_table"][0]["direction_passed"] is True
    assert packet["folding_problem_solved"] is False


def test_e61_oligomer_context_activates_assembly_grammar() -> None:
    packet = build_sealed_operator_state_packet(
        target_id="E61_HOMOMERIC_CONTEXT_PROBE",
        target_name="homomeric assembly context probe",
        sequence="MKKVVVLGLGLVGLAAGAAADAFKKAGVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
        sources=[
            {
                "source_id": "E61_ASSEMBLY_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "evidence_statement": "RCSB non-coordinate metadata exposes homomeric_context oligomer_context assembly_context partner_copy_context.",
            }
        ],
        perturbations=[
            {
                "perturbation_id": "E61_INTERFACE_DAMAGE",
                "description": "damage assembly interface",
                "interface_disruption": 0.60,
                "metric": "basin:assembly_stabilized_basin",
                "expected_direction": "decrease",
            }
        ],
    )
    assert packet["selected_mechanism_grammar"]["mechanism_class"] == "oligomerization_controlled_folding"
    final = packet["operator_state_propagation_summary"]["final_state_summary"]
    assert final["state_basin_occupancy"]["assembly_stabilized_basin"] > final["state_basin_occupancy"]["interface_rejected_basin"]
    assert any(row["interaction_type"] == "assembly_stabilized_interface" for row in packet["predicted_contact_interaction_probability_map"])
    assert packet["predicted_perturbation_table"][0]["direction_passed"] is True
