from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.protein_esperanto_engine import build_sealed_operator_state_packet


MEMBRANE_RICH_SEQUENCE = "M" + ("LLLLIIVVFF" * 12) + "GSGS" + ("AKTQ" * 10)


def _packet(statement: str) -> dict[str, object]:
    return build_sealed_operator_state_packet(
        target_id="E63_MEMBRANE_TOPOLOGY_PROBE",
        target_name="E63 topology probe",
        sequence=MEMBRANE_RICH_SEQUENCE,
        sources=[
            {
                "source_id": "E63_TOPOLOGY_SOURCE",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "evidence_statement": statement,
            }
        ],
    )


def _selected(statement: str) -> str:
    packet = _packet(statement)
    return packet["selected_mechanism_grammar"]["mechanism_class"]


def test_e63_preserves_true_transmembrane_topology_routing() -> None:
    assert (
        _selected("membrane_context_strong transmembrane_context topology_evidence transmembrane helix")
        == "membrane_multidomain_folding_proteostasis"
    )


def test_e63_abstains_on_hydrophobicity_without_topology() -> None:
    assert (
        _selected(
            "soluble_hydrophobic_core_context hydrophobicity-alone signal; "
            "no membrane topology evidence and no OPM/PDBTM/MemProtMD transmembrane assignment."
        )
        == "insufficient_evidence_clean_abstain"
    )


def test_e63_abstains_on_peripheral_or_monotopic_membrane_context() -> None:
    assert (
        _selected(
            "OPM monotopic/peripheral source: peripheral membrane association, "
            "amphipathic peripheral helix, not transmembrane, no bilayer-spanning topology evidence."
        )
        == "insufficient_evidence_clean_abstain"
    )


def test_e63_keeps_cofactor_and_oligomer_explanations_for_hydrophobicity() -> None:
    assert (
        _selected("cofactor_context ligand_context cofactor-buried hydrophobic pocket; no membrane topology evidence.")
        == "cofactor_ligand_assisted_stabilization"
    )
    assert (
        _selected("oligomer_context assembly_context partner_copy_context oligomeric interface hydrophobicity; no membrane topology evidence.")
        == "oligomerization_controlled_folding"
    )


def test_e63_revision_certificate_records_v65_failure_source() -> None:
    cert_path = ROOT / "data" / "protein_esperanto_engine" / "E63" / "e63_v65_membrane_topology_grammar_revision_certificate.json"
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    assert cert["engine_version"] == "E63"
    assert cert["derived_from_batch"] == "V65_MEMBRANE_TOPOLOGY_ESPERANTO_PANEL"
    assert cert["failure_modes_addressed"] == [
        "soluble_hydrophobic_false_membrane",
        "peripheral_misread_as_transmembrane",
    ]
    assert cert["engine_modified_after_v65"] is True
