from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.protein_esperanto_engine import MECHANISM_CLASSES, build_sealed_operator_state_packet


SEQUENCE = "LRVQPEAQAKVDVFREDLCTKTENLLGSYFPKKISELDAFLKEPALNEANLSNLKAPLDI"
MEMBRANE_SEQUENCE = "M" + ("LLLLIIVVFF" * 12) + "GSGS" + ("AKTQ" * 10)


def _selected(statement: str, sequence: str = SEQUENCE) -> str:
    packet = build_sealed_operator_state_packet(
        target_id="E65_ASSEMBLY_PROBE",
        target_name="E65 assembly probe",
        sequence=sequence,
        sources=[
            {
                "source_id": "E65_CONTEXT",
                "source_class": "pure_non_coordinate",
                "source_role": "prediction_input",
                "coordinate_derived": False,
                "internal_runtime_source": False,
                "evidence_statement": statement,
            }
        ],
    )
    return packet["selected_mechanism_grammar"]["mechanism_class"]


def test_e65_adds_assembly_required_folding_class() -> None:
    assert "assembly_required_folding" in MECHANISM_CLASSES


def test_e65_routes_explicit_partner_completed_core_to_assembly_required() -> None:
    assert (
        _selected(
            "assembly_required_core assembly_required_folding partner_completed_core "
            "interface_buried_hydrophobicity monomer_incomplete_topology biological_oligomer_context"
        )
        == "assembly_required_folding"
    )


def test_e65_preserves_true_transmembrane_priority_over_assembly() -> None:
    assert (
        _selected(
            "membrane_context_strong transmembrane_context topology_evidence channel_context "
            "assembly_required_core partner_completed_core biological_oligomer_context",
            MEMBRANE_SEQUENCE,
        )
        == "membrane_multidomain_folding_proteostasis"
    )


def test_e65_generic_complex_only_abstains() -> None:
    assert (
        _selected("generic_complex_only generic complex only; complex annotation alone is not obligate assembly")
        == "insufficient_evidence_clean_abstain"
    )


def test_e65_negated_assembly_tokens_do_not_trigger_false_assembly() -> None:
    assert (
        _selected("cofactor_context ligand_context cofactor-stabilized soluble core; not assembly_required_core")
        == "cofactor_ligand_assisted_stabilization"
    )
    assert (
        _selected(
            "soluble_monomeric_core_context complete soluble monomer hydrophobic core; "
            "not assembly_required_folding and not biological_oligomer_context"
        )
        == "globular_closure"
    )


def test_e65_revision_certificate_records_v67_failure_source() -> None:
    cert = json.loads(
        (
            ROOT
            / "data"
            / "protein_esperanto_engine"
            / "E65"
            / "e65_assembly_required_folding_grammar_certificate.json"
        ).read_text(encoding="utf-8")
    )
    assert cert["engine_version"] == "E65"
    assert cert["source_batch_trigger"] == "V67_RCSB_NONREDUNDANT_200_DISCOVERY_E64"
    assert cert["failure_mode_addressed"] == "assembly_required_core_vs_membrane_topology"
    assert cert["missing_esperanto_word"] == "assembly_required_core_vs_topology_provider"
    assert cert["next_required_batch"] == "V68_OLIGOMER_ASSEMBLY_PANEL_200"
