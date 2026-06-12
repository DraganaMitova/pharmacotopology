from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pharmacotopology.protein_esperanto_engine import MECHANISM_CLASSES, STATE_VARIABLES, build_sealed_simulation_packet


def _source(statement: str, *, source_id: str = "E68_TEST_SOURCE", marks: list[str] | None = None, withheld: list[str] | None = None) -> dict[str, object]:
    return {
        "source_id": source_id,
        "source_class": "pure_non_coordinate",
        "source_role": "prediction_input",
        "coordinate_derived": False,
        "internal_runtime_source": False,
        "spatial_proxy": False,
        "metadata_context_marks": marks or [],
        "withheld_context_marks": withheld or [],
        "evidence_statement": statement,
    }


def _packet(statement: str, *, marks: list[str] | None = None, withheld: list[str] | None = None) -> dict[str, object]:
    sequence = ("VIFYWT" * 45)[:240]
    return build_sealed_simulation_packet(
        target_id="E68_BETA_CLOSURE_TEST",
        target_name="E68 beta closure test",
        sequence=sequence,
        sources=[_source(statement, marks=marks, withheld=withheld)],
        perturbations=[],
    )


def test_e68_adds_beta_closure_mechanism_and_state_words() -> None:
    assert "beta_closure_topology" in MECHANISM_CLASSES
    for word in [
        "closed_beta_topology",
        "strand_register",
        "beta_sheet_closure",
        "soluble_beta_barrel",
        "membrane_beta_barrel",
        "beta_propeller_repeat_closure",
        "beta_sandwich_core",
        "jelly_roll_wrap",
        "greek_key_beta_lock",
        "beta_helix_solenoid_stack",
        "alpha_beta_barrel_distinction",
        "open_beta_sheet_ambiguous",
        "beta_topology_conflict",
    ]:
        assert word in STATE_VARIABLES

    cert = json.loads(
        (
            ROOT
            / "data"
            / "protein_esperanto_engine"
            / "E68"
            / "e68_beta_closure_topology_grammar_certificate.json"
        ).read_text(encoding="utf-8")
    )
    assert cert["engine_revision"] == "E68"
    assert cert["new_mechanism_class"] == "beta_closure_topology"
    assert cert["v72_closed_beta_failures_repaired"] == 30


def test_e68_repairs_closed_beta_vs_generic_oligomer() -> None:
    packet = _packet("closed_beta_topology strand_register beta_sheet_closure oligomer_context assembly_context partner_copy_context")
    assert packet["selected_mechanism_grammar"]["mechanism_class"] == "beta_closure_topology"
    assert packet["selected_mechanism_grammar"]["selected_beta_topology_word"] == "closed_beta_topology"
    final = packet["trajectory_summary"]["final_state_summary"]
    assert final["closed_beta_topology"] > 0.0
    assert final["strand_register"] > 0.0
    assert final["beta_sheet_closure"] > 0.0
    assert any(row["interaction_type"] == "beta_strand_register_closure" for row in packet["predicted_contact_interaction_probability_map"])


def test_e68_selects_beta_subtypes_and_abstains_on_open_beta_only() -> None:
    sandwich = _packet("beta_sandwich_core closed_beta_topology strand_register beta_sheet_closure")
    membrane_barrel = _packet("membrane_beta_barrel membrane_context_strong transmembrane_context closed_beta_topology strand_register")
    ambiguous = _packet("open_beta_sheet_ambiguous beta propensity only strand_register_insufficient")

    assert sandwich["selected_mechanism_grammar"]["selected_beta_topology_word"] == "beta_sandwich_core"
    assert sandwich["trajectory_summary"]["final_state_summary"]["beta_sandwich_core"] > 0.0
    assert membrane_barrel["selected_mechanism_grammar"]["mechanism_class"] == "beta_closure_topology"
    assert membrane_barrel["selected_mechanism_grammar"]["selected_beta_topology_word"] == "membrane_beta_barrel"
    assert membrane_barrel["trajectory_summary"]["final_state_summary"]["membrane_beta_barrel"] > 0.0
    assert ambiguous["selected_mechanism_grammar"]["mechanism_class"] == "insufficient_evidence_clean_abstain"


def test_e68_preserves_priorities_and_keeps_withheld_context_invisible() -> None:
    membrane = _packet("membrane_context_strong transmembrane_context topology_evidence closed_beta_topology")
    assembly = _packet("assembly_required_core assembly_required_folding partner_completed_core closed_beta_topology")
    metal = _packet("metal_cluster_geometry ligand_locked_basin closed_beta_topology")
    disorder = _packet("disorder_context IDR_boundary structured_domain_plus_IDR_tail closed_beta_topology")
    withheld = _packet(
        "closed_beta_topology strand_register beta_sheet_closure",
        withheld=["assembly_required_core", "membrane_beta_barrel"],
    )

    assert membrane["selected_mechanism_grammar"]["mechanism_class"] == "membrane_multidomain_folding_proteostasis"
    assert assembly["selected_mechanism_grammar"]["mechanism_class"] == "assembly_required_folding"
    assert metal["selected_mechanism_grammar"]["mechanism_class"] == "metal_cluster_and_ligand_locked_basin"
    assert disorder["selected_mechanism_grammar"]["mechanism_class"] == "disorder_boundary_and_fold_upon_binding"
    assert withheld["selected_mechanism_grammar"]["mechanism_class"] == "beta_closure_topology"
