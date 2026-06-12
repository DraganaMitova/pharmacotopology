from __future__ import annotations

import json
from pathlib import Path

from pharmacotopology.protein_esperanto_engine import (
    E78_CLAIM_LADDER,
    final_claim_ledger,
    final_closed_loop_packet,
    universal_claim_firewall,
)


ROOT = Path(__file__).resolve().parents[1]
E78_ROOT = ROOT / "data" / "protein_esperanto_engine" / "E78"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_e78_claim_ledger_unlocks_highest_honest_bounded_claim() -> None:
    ledger = final_claim_ledger(
        language_claim_count=6000,
        independent_observable_support_count=3000,
        mechanism_claim_count=0,
        topology_claim_count=1000,
        bounded_physical_claim_count=1000,
        target_fold_claim_count=0,
        general_solution_candidate_claim_allowed=False,
        universal_solution_claim_allowed=False,
        unsupported_physical_claims=0,
        unsupported_topology_claims=0,
        unsupported_fold_claims=0,
        failed_accepted_count=0,
    )

    assert ledger["claim_ladder"] == E78_CLAIM_LADDER
    assert ledger["highest_claim_tier_unlocked"] == "CLAIM_4_COARSE_PHYSICAL_SUPPORTED"
    assert ledger["final_status"] == "FINAL_LOOP_CLOSED_BOUNDED_PHYSICAL_LANGUAGE"
    assert ledger["protein_folding_solved"] is False


def test_e78_universal_firewall_blocks_without_target_fold_subset() -> None:
    firewall = universal_claim_firewall(
        broad_fresh_generalization_passed=False,
        independent_physical_observables_passed=True,
        target_fold_subset_passed=False,
        unsupported_claims=0,
        failed_accepted_count=0,
    )

    assert firewall["kind"] == "E78_UNIVERSAL_CLAIM_FIREWALL_v0"
    assert firewall["universal_solution_claim_allowed"] is False
    assert firewall["protein_folding_solved"] is False
    assert firewall["blocked_reasons"] == [
        "broad_fresh_generalization_not_sufficient",
        "target_fold_subset_not_sufficient",
    ]


def test_e78_final_packet_seals_prediction_before_holdouts_and_blocks_fold_claim() -> None:
    packet = final_closed_loop_packet(
        target_id="E78_PACKET_TEST",
        sequence_hash="abc123",
        target_family="coarse_physical_execution",
        sealed_language_prediction={"decision": "accepted", "token_only_acceptance": False},
        sealed_sentence_prediction={"decision": "accepted_supported"},
        claim_tier_unlocked="CLAIM_4_COARSE_PHYSICAL_SUPPORTED",
        claim_allowed=True,
        claim_blocked_reason=None,
        supports_selected=True,
        proxy_physical_execution=True,
    )

    assert packet["kind"] == "E78_FINAL_CLOSED_LOOP_PACKET_v0"
    assert packet["prediction_hash"]
    assert all(value is False for value in packet["preseal_allowed_sources"].values())
    assert packet["postseal_biology_holdout"]["opened_after_prediction_hash"] is True
    assert packet["postseal_biology_holdout"]["used_before_prediction"] is False
    assert packet["postseal_biology_holdout"]["observable_supports_selected_sentence"] is True
    assert packet["postseal_physical_holdout"]["kind"] == "E78_PHYSICAL_EXECUTION_BINDING_v0"
    assert packet["postseal_physical_holdout"]["proxy_physical_execution"] is True
    assert packet["postseal_physical_holdout"]["postseal_physical_holdout_is_proxy"] is True
    assert packet["postseal_physical_holdout"]["real_physical_execution_data_source"] is None
    assert packet["postseal_physical_holdout"]["openmm_execution_performed"] is False
    assert packet["postseal_physical_holdout"]["validated_coarse_execution_performed"] is False
    assert packet["postseal_physical_holdout"]["physical_claim_ceiling"] == "CLAIM_4_COARSE_PHYSICAL_SUPPORTED"
    assert packet["postseal_physical_holdout"]["physical_basis_claim_allowed"] is False
    assert packet["target_fold_support_packet"]["target_fold_claim_allowed"] is False
    assert packet["target_fold_support_packet"]["structural_knot_claim_allowed"] is False
    assert packet["physical_basis_claim_allowed"] is False
    assert packet["protein_folding_solved"] is False


def test_e78_certificate_records_final_engine_objects() -> None:
    cert = _read(E78_ROOT / "e78_closed_loop_biophysical_protein_language_engine_certificate.json")

    assert cert["kind"] == "E78_CLOSED_LOOP_BIOPHYSICAL_PROTEIN_LANGUAGE_ENGINE_CERTIFICATE_v0"
    assert cert["engine_revision"] == "E78"
    assert cert["baseline_engine_revision"] == "E77"
    assert cert["claim_ladder"] == E78_CLAIM_LADDER
    assert cert["final_status"] == "FINAL_LOOP_CLOSED_BOUNDED_PHYSICAL_LANGUAGE"
    assert cert["highest_claim_tier_unlocked"] == "CLAIM_4_COARSE_PHYSICAL_SUPPORTED"
    assert cert["final_closed_loop_packet_kind"] == "E78_FINAL_CLOSED_LOOP_PACKET_v0"
    assert cert["independent_holdout_binding_kind"] == "E78_INDEPENDENT_HOLDOUT_BINDING_v0"
    assert cert["physical_execution_binding_kind"] == "E78_PHYSICAL_EXECUTION_BINDING_v0"
    assert cert["target_fold_support_packet_kind"] == "E78_TARGET_FOLD_SUPPORT_PACKET_v0"
    assert cert["universal_claim_firewall_kind"] == "E78_UNIVERSAL_CLAIM_FIREWALL_v0"
    assert cert["failure_autopsy_packet_kind"] == "E78_FAILURE_AUTOPSY_PACKET_v0"
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
