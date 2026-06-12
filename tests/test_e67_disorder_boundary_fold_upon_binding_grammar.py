from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.protein_esperanto_engine import MECHANISM_CLASSES, build_sealed_operator_state_packet


def _source(source_id: str, statement: str) -> dict[str, object]:
    return {
        "source_id": source_id,
        "source_class": "pure_non_coordinate",
        "source_role": "prediction_input",
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "spatial_proxy": False,
        "evidence_statement": statement,
    }


def _packet(statement: str) -> dict[str, object]:
    sequence = "M" + "GSQNPY" * 35 + "KRDE" * 15
    return build_sealed_operator_state_packet(
        target_id="E67_DISORDER_BOUNDARY_TEST",
        target_name="E67 disorder boundary test",
        sequence=sequence,
        sources=[_source("E67_TEST_SOURCE", statement)],
        perturbations=[],
    )


def test_e67_adds_disorder_boundary_mechanism_class() -> None:
    assert "disorder_boundary_and_fold_upon_binding" in MECHANISM_CLASSES
    cert = json.loads(
        (
            ROOT
            / "data"
            / "protein_esperanto_engine"
            / "E67"
            / "e67_disorder_boundary_fold_upon_binding_grammar_certificate.json"
        ).read_text(encoding="utf-8")
    )
    assert cert["engine_revision"] == "E67"
    assert cert["new_mechanism_class"] == "disorder_boundary_and_fold_upon_binding"
    assert cert["v71_disorder_failures_repaired"] == 31


def test_e67_disorder_boundary_outranks_generic_oligomer() -> None:
    packet = _packet(
        "disorder_context phase_prone_low_complexity oligomer_context "
        "assembly_context partner_copy_context low complexity tendency"
    )
    assert packet["selected_mechanism_grammar"]["mechanism_class"] == "disorder_boundary_and_fold_upon_binding"
    final = packet["operator_state_propagation_summary"]["final_state_summary"]
    assert final["IDR_boundary"] > 0
    assert final["phase_prone_low_complexity"] > 0


def test_e67_preserves_higher_priority_contexts() -> None:
    tm = _packet("membrane_context_strong transmembrane_context topology_evidence disorder_context")
    assembly = _packet("assembly_required_core assembly_required_folding partner_completed_core disorder_context")
    metal = _packet("metal_cluster_geometry disorder_context")

    assert tm["selected_mechanism_grammar"]["mechanism_class"] == "membrane_multidomain_folding_proteostasis"
    assert assembly["selected_mechanism_grammar"]["mechanism_class"] == "assembly_required_folding"
    assert metal["selected_mechanism_grammar"]["mechanism_class"] == "metal_cluster_and_ligand_locked_basin"
