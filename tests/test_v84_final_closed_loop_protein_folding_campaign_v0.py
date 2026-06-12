from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V84_ROOT = ROOT / "data" / "protein_esperanto_engine" / "V84"


def _read(path: Path) -> dict[str, object]:
    assert path.exists()
    return json.loads(path.read_text(encoding="utf-8"))


def test_v84_certificate_reports_highest_honest_final_status() -> None:
    cert = _read(V84_ROOT / "v84_final_closed_loop_protein_folding_campaign_certificate.json")

    assert cert["kind"] == "V84_FINAL_CLOSED_LOOP_PROTEIN_FOLDING_CAMPAIGN_CERTIFICATE_v0"
    assert cert["status"] == "FINAL_LOOP_CLOSED_BOUNDED_PHYSICAL_LANGUAGE"
    assert cert["campaign_controls_passed"] is True
    assert cert["engine_version_used"] == "E78"
    assert cert["baseline_engine_version"] == "E77"
    assert cert["targets_total"] == 9000
    assert cert["row_family_counts"] == {
        "compositional_sentence_target": 1000,
        "fresh_blind_nonredundant_target": 3000,
        "hard_abstain_no_evidence_metadata_masked_control": 1000,
        "old_sentinel_all_grammar_family": 1000,
        "physical_execution_target": 1000,
        "token_soup_adversarial_control": 1000,
        "topology_family_holdout_target": 1000,
    }
    assert cert["fresh_target_shortage"] is True
    assert cert["fresh_source_count_available"] == 1000
    assert cert["deterministic_variant_count"] == 2000
    assert cert["claim_ceiling_due_to_cache"] == "CLAIM_4_COARSE_PHYSICAL_SUPPORTED"
    assert cert["final_claim_ledger"]["highest_claim_tier_unlocked"] == "CLAIM_4_COARSE_PHYSICAL_SUPPORTED"
    assert cert["independent_observable_support_count"] == 3000
    assert cert["bounded_physical_claim_count"] == 1000
    assert cert["topology_claim_count"] == 1000
    assert cert["target_fold_claim_count"] == 0
    assert cert["general_solution_candidate_claim_allowed"] is False
    assert cert["universal_folding_solution_claim_allowed"] is False
    assert cert["protein_esperanto_language_complete_for_current_registry"] is True
    assert cert["failed_accepted_count"] == 0
    assert cert["token_only_acceptance_count"] == 0
    assert cert["coordinate_native_leakage"] is False
    assert cert["sentinels_preserved"] == 1000
    assert cert["hard_abstains_clean"] is True
    assert cert["unsupported_claims"] == 0
    assert cert["unsupported_physical_claims"] == 0
    assert cert["unsupported_topology_claims"] == 0
    assert cert["unsupported_fold_claims"] == 0
    assert cert["bounded_physical_language_claim_allowed"] is True
    assert cert["physical_basis_claim_allowed"] is False
    assert cert["protein_folding_solved"] is False
    assert cert["claim_blocked_reason"] == "fresh_target_shortage_and_proxy_physical_execution_cap_claim_ceiling_at_CLAIM_4"
    assert cert["failure_autopsy_count"] == 5500
    assert cert["failed_controls"] == []


def test_v84_report_contains_final_packets_with_sealed_holdouts_and_controls() -> None:
    report = _read(V84_ROOT / "v84_final_closed_loop_protein_folding_campaign_report.json")
    rows = report["rows"]

    assert len(rows) == 9000
    assert len(report["failure_autopsies"]) == 5500
    assert all("closed_loop_packet" in row for row in rows)
    assert all(row["closed_loop_packet"]["kind"] == "E78_FINAL_CLOSED_LOOP_PACKET_v0" for row in rows)
    assert all(row["closed_loop_packet"]["prediction_hash"] for row in rows)
    assert all(
        all(value is False for value in row["closed_loop_packet"]["preseal_allowed_sources"].values())
        for row in rows
    )
    assert all(
        row["closed_loop_packet"]["postseal_biology_holdout"]["opened_after_prediction_hash"] is True
        and row["closed_loop_packet"]["postseal_biology_holdout"]["used_before_prediction"] is False
        for row in rows
    )
    assert all(row["closed_loop_packet"]["wrong_grammar_controls"]["supports_selected"] is False for row in rows)
    assert all(row["closed_loop_packet"]["bag_of_words_controls"]["supports_selected"] is False for row in rows)
    assert all(row["closed_loop_packet"]["masked_controls"]["supports_selected"] is False for row in rows)
    assert all(row["closed_loop_packet"]["wrong_target_controls"]["supports_selected"] is False for row in rows)
    assert all(row["closed_loop_packet"]["physical_basis_claim_allowed"] is False for row in rows)
    assert all(row["closed_loop_packet"]["protein_folding_solved"] is False for row in rows)

    physical_rows = [row for row in rows if row["row_family"] == "physical_execution_target"]
    assert len(physical_rows) == 1000
    assert sum(1 for row in physical_rows if row.get("atomistic_or_validated_coarse_subset_member")) == 500
    assert all(row.get("subset_execution_mode") == "deterministic_equivalent_proxy" for row in physical_rows[:500])
    assert all(
        row["closed_loop_packet"]["claim_tier_unlocked"] == "CLAIM_4_COARSE_PHYSICAL_SUPPORTED"
        for row in physical_rows
    )
    assert all(
        row["closed_loop_packet"]["target_fold_support_packet"]["target_fold_claim_allowed"] is False
        for row in physical_rows
    )
    assert all(
        row["closed_loop_packet"]["postseal_physical_holdout"]["postseal_physical_holdout_is_proxy"] is True
        for row in physical_rows
    )
    assert all(
        row["closed_loop_packet"]["postseal_physical_holdout"]["real_physical_execution_data_source"] is None
        for row in physical_rows
    )
    assert all(
        row["closed_loop_packet"]["postseal_physical_holdout"]["openmm_execution_performed"] is False
        for row in physical_rows
    )
    assert all(
        row["closed_loop_packet"]["postseal_physical_holdout"]["validated_coarse_execution_performed"] is False
        for row in physical_rows
    )
    assert all(
        row["closed_loop_packet"]["postseal_physical_holdout"]["physical_claim_ceiling"]
        == "CLAIM_4_COARSE_PHYSICAL_SUPPORTED"
        for row in physical_rows
    )
    assert all(
        row["closed_loop_packet"]["postseal_physical_holdout"]["physical_basis_claim_allowed"] is False
        for row in physical_rows
    )
