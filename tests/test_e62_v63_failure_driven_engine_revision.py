from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.protein_esperanto_engine import build_sealed_simulation_packet


MEMBRANE_RICH_SEQUENCE = "M" + ("LLLLIIVVFF" * 12) + "GSGS" + ("AKTQ" * 10)
SOLUBLE_SEQUENCE = (
    "MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWAEYFHQKALDDFKNAKPNKTDPNKLAEQLKQLE"
    "EQLKQLEQAQKQLEQAKQ"
)


def test_e62_membrane_context_precedes_incidental_oligomer_context() -> None:
    packet = build_sealed_simulation_packet(
        target_id="E62_MEMBRANE_OLIGOMER_PRIORITY_PROBE",
        target_name="membrane plus incidental assembly context",
        sequence=MEMBRANE_RICH_SEQUENCE,
        sources=[
            {
                "source_id": "E62_MEMBRANE_AND_ASSEMBLY_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "evidence_statement": "Public metadata exposes membrane_context_strong transmembrane_context oligomer_context assembly_context partner_copy_context homomeric_context.",
            }
        ],
    )
    selected = packet["selected_mechanism_grammar"]
    assert selected["mechanism_class"] == "membrane_multidomain_folding_proteostasis"
    assert "prioritized" in selected["selection_reason"]


def test_e62_membrane_context_precedes_incidental_cofactor_context() -> None:
    packet = build_sealed_simulation_packet(
        target_id="E62_MEMBRANE_COFATOR_PRIORITY_PROBE",
        target_name="membrane plus incidental metal context",
        sequence=MEMBRANE_RICH_SEQUENCE,
        sources=[
            {
                "source_id": "E62_MEMBRANE_AND_METAL_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "evidence_statement": "Public metadata exposes membrane_context_strong transmembrane_context cofactor_context ligand_context metal_context.",
            }
        ],
    )
    assert packet["selected_mechanism_grammar"]["mechanism_class"] == "membrane_multidomain_folding_proteostasis"


def test_e62_preserves_pure_cofactor_and_oligomer_routing() -> None:
    cofactor = build_sealed_simulation_packet(
        target_id="E62_PURE_COFATOR_CONTEXT_PROBE",
        target_name="pure cofactor context",
        sequence=SOLUBLE_SEQUENCE,
        sources=[
            {
                "source_id": "E62_COFATOR_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "evidence_statement": "Public metadata exposes cofactor_context ligand_context metal_context with no membrane topology context.",
            }
        ],
    )
    oligomer = build_sealed_simulation_packet(
        target_id="E62_PURE_OLIGOMER_CONTEXT_PROBE",
        target_name="pure oligomer context",
        sequence=SOLUBLE_SEQUENCE,
        sources=[
            {
                "source_id": "E62_OLIGOMER_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "evidence_statement": "Public metadata exposes oligomer_context assembly_context partner_copy_context homomeric_context with no membrane topology context.",
            }
        ],
    )
    assert cofactor["selected_mechanism_grammar"]["mechanism_class"] == "cofactor_ligand_assisted_stabilization"
    assert oligomer["selected_mechanism_grammar"]["mechanism_class"] == "oligomerization_controlled_folding"


def test_e62_revision_certificate_records_v63_failure_source() -> None:
    cert_path = ROOT / "data" / "protein_esperanto_engine" / "E62" / "e62_v63_failure_driven_engine_revision_certificate.json"
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    assert cert["engine_version"] == "E62"
    assert cert["derived_from_batch"] == "V63_RCSB_500_DISCOVERY_BATCH"
    assert cert["failure_modes_addressed"] == ["membrane_misread"]
    assert cert["engine_modified_after_v63"] is True
