from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.protein_esperanto_engine import MECHANISM_CLASSES, build_sealed_simulation_packet


SEQUENCE = "ACDEFGHIKLMNPQRSTVWY" * 8
MEMBRANE_SEQUENCE = "M" + ("LLLLIIVVFF" * 12) + "GSGS" + ("AKTQ" * 10)


def _packet(statement: str, sequence: str = SEQUENCE) -> dict[str, object]:
    return build_sealed_simulation_packet(
        target_id="E66_METAL_PROBE",
        target_name="E66 metal probe",
        sequence=sequence,
        sources=[
            {
                "source_id": "E66_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "evidence_statement": statement,
            }
        ],
    )


def _selected(statement: str, sequence: str = SEQUENCE) -> str:
    return _packet(statement, sequence)["selected_mechanism_grammar"]["mechanism_class"]


def test_e66_adds_metal_ligand_locked_basin_class() -> None:
    assert "metal_cluster_and_ligand_locked_basin" in MECHANISM_CLASSES


def test_e66_routes_metal_and_ligand_locked_evidence_to_new_class() -> None:
    metal = _packet("cofactor_context ligand_context metal_context heme_context metal_cluster_geometry")
    assert metal["selected_mechanism_grammar"]["mechanism_class"] == "metal_cluster_and_ligand_locked_basin"
    assert metal["trajectory_summary"]["final_state_summary"]["metal_cluster_geometry"] > 0.0
    ligand = _packet("cofactor_context ligand_context ligand_locked_basin apo_holo_basin_shift")
    assert ligand["selected_mechanism_grammar"]["mechanism_class"] == "metal_cluster_and_ligand_locked_basin"
    assert ligand["trajectory_summary"]["final_state_summary"]["ligand_locked_basin"] > 0.0


def test_e66_keeps_generic_cofactor_generic() -> None:
    assert (
        _selected("cofactor_context ligand_context generic cofactor evidence without metal or locked basin")
        == "cofactor_ligand_assisted_stabilization"
    )


def test_e66_preserves_membrane_and_assembly_priority() -> None:
    assert (
        _selected("membrane_context_strong transmembrane_context topology_evidence metal_context", MEMBRANE_SEQUENCE)
        == "membrane_multidomain_folding_proteostasis"
    )
    assert (
        _selected("assembly_required_core partner_completed_core biological_oligomer_context metal_context")
        == "assembly_required_folding"
    )


def test_e66_certificate_records_v69_dominant_failure() -> None:
    cert = json.loads(
        (
            ROOT
            / "data"
            / "protein_esperanto_engine"
            / "E66"
            / "e66_metal_cluster_ligand_locked_basin_grammar_certificate.json"
        ).read_text(encoding="utf-8")
    )
    assert cert["engine_version"] == "E66"
    assert cert["source_batch_trigger"] == "V69_RCSB_NONREDUNDANT_200_DISCOVERY_E65"
    assert cert["failure_mode_addressed"] == "metal_cluster_geometry"
    assert cert["secondary_failure_mode_addressed"] == "ligand_locked_basin"
    assert cert["new_mechanism_class"] == "metal_cluster_and_ligand_locked_basin"
    assert cert["next_required_batch"] == "V70_E66_METAL_CLUSTER_LIGAND_REPAIR_PANEL_200"
