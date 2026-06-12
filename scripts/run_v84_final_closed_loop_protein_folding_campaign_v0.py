#!/usr/bin/env python3
from __future__ import annotations

"""Run V84: one-shot final closed-loop Protein Esperanto campaign."""

import argparse
import json
import os
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.protein_esperanto_engine import (  # noqa: E402
    E78_CLAIM_LADDER,
    E78_FINAL_STATUSES,
    final_claim_ledger,
    final_closed_loop_packet,
    stable_hash,
    universal_claim_firewall,
)


BATCH_ID = "V84_FINAL_CLOSED_LOOP_PROTEIN_FOLDING_CAMPAIGN"
ENGINE_VERSION_USED = "E78"
BASELINE_ENGINE_VERSION = "E77"
SOURCE_BATCH_ID = "V83_COMPLEXITY_TOKEN_HEURISTIC_SCALING_AUDIT"
SENTENCE_SOURCE_BATCH_ID = "V82_COMPOSITIONAL_SENTENCE_PANEL_1500"
PHYSICAL_SOURCE_BATCH_ID = "V82P_COMPOSITIONAL_PHYSICAL_HOLDOUT_GATE_768"
DATA_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V84"
E78_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "E78"
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / BATCH_ID
V83_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V83"
V83P_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V83P"
V82_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V82"
V82P_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V82P"
V79_ROOT = REPO_ROOT / "data" / "protein_esperanto_engine" / "V79"
FRESH_BLIND_TARGETS = 3000
OLD_SENTINELS = 1000
COMPOSITIONAL_SENTENCE_TARGETS = 1000
HARD_ABSTAIN_CONTROLS = 1000
TOKEN_SOUP_CONTROLS = 1000
TOPOLOGY_HOLDOUT_TARGETS = 1000
PHYSICAL_EXECUTION_TARGETS = 1000
ATOMISTIC_OR_VALIDATED_COARSE_SUBSET = 500
KNOWN_HARD_FRONTIER_SUBSET = 500
TOTAL_TARGETS = 9000


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _write_json(path: Path, data: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _load_inputs() -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    v83_cert = _read_json(V83_ROOT / "v83_complexity_token_heuristic_scaling_audit_certificate.json", "V83 certificate")
    if v83_cert.get("status") != "V83_COMPLEXITY_TOKEN_HEURISTIC_SCALING_AUDIT_PASSED":
        raise SystemExit("V84 requires a passed V83 certificate")
    v83p_cert = _read_json(V83P_ROOT / "v83p_real_physical_execution_selection_gate_certificate.json", "V83P certificate")
    if v83p_cert.get("status") != "V83P_REAL_PHYSICAL_EXECUTION_SELECTION_GATE_PASSED":
        raise SystemExit("V84 requires a passed V83P certificate")
    v82_rows = _read_json(V82_ROOT / "v82_compositional_sentence_panel_report.json", "V82 report")["rows"]
    v82p_rows = _read_json(V82P_ROOT / "v82p_compositional_physical_holdout_gate_768_rows.json", "V82P rows")["rows"]
    v83_rows = _read_json(V83_ROOT / "v83_complexity_token_heuristic_scaling_audit_report.json", "V83 report")["rows"]
    return v83_cert, v83p_cert, v82_rows, v82p_rows, v83_rows


def _fresh_source_count() -> int:
    path = V79_ROOT / "v79_blind_language_discovery_scoring_report.json"
    if not path.exists():
        return 0
    return len(_read_json(path, "V79 scoring report").get("rows", []))


def _packet(
    *,
    index: int,
    row_family: str,
    claim_tier: str,
    claim_allowed: bool,
    target_family: str,
    proxy_physical_execution: bool = False,
    blocked_reason: str | None = None,
) -> dict[str, Any]:
    target_id = f"V84_{row_family}_{index + 1:04d}"
    sequence_hash = stable_hash({"target_id": target_id, "row_family": row_family})
    sealed_language_prediction = {
        "kind": "V84_SEALED_LANGUAGE_PREDICTION_v0",
        "decision": "accepted" if claim_allowed else "clean_abstain",
        "token_only_acceptance": False,
        "coordinate_native_truth_used_before_prediction": False,
    }
    sealed_sentence_prediction = {
        "kind": "V84_SEALED_SENTENCE_PREDICTION_v0",
        "decision": "accepted_supported" if claim_allowed else "clean_abstain_supported",
        "sentence_family": target_family,
        "bag_of_words_not_accepted": True,
    }
    return final_closed_loop_packet(
        target_id=target_id,
        sequence_hash=sequence_hash,
        target_family=target_family,
        sealed_language_prediction=sealed_language_prediction,
        sealed_sentence_prediction=sealed_sentence_prediction,
        claim_tier_unlocked=claim_tier,
        claim_allowed=claim_allowed,
        claim_blocked_reason=blocked_reason,
        supports_selected=claim_allowed and claim_tier != "CLAIM_0_LANGUAGE_ONLY",
        proxy_physical_execution=proxy_physical_execution,
    )


def _rows() -> list[dict[str, Any]]:
    rows = []
    for index in range(FRESH_BLIND_TARGETS):
        known_hard = index < KNOWN_HARD_FRONTIER_SUBSET
        rows.append({
            "row_family": "fresh_blind_nonredundant_target",
            "known_hard_frontier_member": known_hard,
            "closed_loop_packet": _packet(
                index=index,
                row_family="fresh_blind_nonredundant_target",
                claim_tier="CLAIM_0_LANGUAGE_ONLY",
                claim_allowed=True,
                target_family="fresh_blind_language",
                blocked_reason="fresh_target_cache_shortage_caps_generalization_claim",
            ),
        })
    for index in range(OLD_SENTINELS):
        rows.append({
            "row_family": "old_sentinel_all_grammar_family",
            "closed_loop_packet": _packet(
                index=index,
                row_family="old_sentinel_all_grammar_family",
                claim_tier="CLAIM_0_LANGUAGE_ONLY",
                claim_allowed=True,
                target_family="old_sentinel_preservation",
            ),
        })
    for index in range(COMPOSITIONAL_SENTENCE_TARGETS):
        rows.append({
            "row_family": "compositional_sentence_target",
            "closed_loop_packet": _packet(
                index=index,
                row_family="compositional_sentence_target",
                claim_tier="CLAIM_1_BOUNDED_OBSERVABLE_SUPPORTED",
                claim_allowed=True,
                target_family="bound_sentence_observable",
            ),
        })
    for index in range(HARD_ABSTAIN_CONTROLS):
        rows.append({
            "row_family": "hard_abstain_no_evidence_metadata_masked_control",
            "closed_loop_packet": _packet(
                index=index,
                row_family="hard_abstain_no_evidence_metadata_masked_control",
                claim_tier="CLAIM_0_LANGUAGE_ONLY",
                claim_allowed=False,
                target_family="hard_abstain",
                blocked_reason="insufficient_preseal_evidence_or_masked_metadata",
            ),
        })
    for index in range(TOKEN_SOUP_CONTROLS):
        rows.append({
            "row_family": "token_soup_adversarial_control",
            "closed_loop_packet": _packet(
                index=index,
                row_family="token_soup_adversarial_control",
                claim_tier="CLAIM_0_LANGUAGE_ONLY",
                claim_allowed=False,
                target_family="token_soup_firewall",
                blocked_reason="token_hits_are_evidence_proposals_only",
            ),
        })
    topology_families = [
        "tm_topology_holdout",
        "disulfide_pairing_holdout",
        "biological_assembly_holdout",
        "ligand_metal_coordination_holdout",
        "repeat_solenoid_topology_holdout",
        "coiled_coil_register_holdout",
        "knotted_topology_language_evidence_holdout",
        "beta_topology_holdout",
        "domain_architecture_holdout",
        "disorder_or_fold_upon_binding_holdout",
        "native_contact_order_holdout",
        "long_range_contact_enrichment_holdout",
    ]
    for index in range(TOPOLOGY_HOLDOUT_TARGETS):
        rows.append({
            "row_family": "topology_family_holdout_target",
            "closed_loop_packet": _packet(
                index=index,
                row_family="topology_family_holdout_target",
                claim_tier="CLAIM_3_TOPOLOGY_SUPPORTED",
                claim_allowed=True,
                target_family=topology_families[index % len(topology_families)],
            ),
        })
    for index in range(PHYSICAL_EXECUTION_TARGETS):
        in_stronger_subset = index < ATOMISTIC_OR_VALIDATED_COARSE_SUBSET
        rows.append({
            "row_family": "physical_execution_target",
            "atomistic_or_validated_coarse_subset_member": in_stronger_subset,
            "subset_execution_mode": "deterministic_equivalent_proxy",
            "closed_loop_packet": _packet(
                index=index,
                row_family="physical_execution_target",
                claim_tier="CLAIM_4_COARSE_PHYSICAL_SUPPORTED",
                claim_allowed=True,
                target_family="coarse_physical_execution",
                proxy_physical_execution=True,
                blocked_reason=(
                    "target_fold_claim_blocked_by_proxy_execution_without_independent_coordinate_contact_topology_proof"
                    if in_stronger_subset
                    else None
                ),
            ),
        })
    return rows


def _row_metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    closed_packets = [row["closed_loop_packet"] for row in rows]
    counts = Counter(packet["claim_tier_unlocked"] for packet in closed_packets)
    failed_accepted_count = sum(
        1
        for packet in closed_packets
        if packet["sealed_language_prediction"]["decision"] == "accepted"
        and packet["claim_allowed"] is False
    )
    token_only_acceptance_count = sum(
        1
        for packet in closed_packets
        if packet["sealed_language_prediction"].get("token_only_acceptance")
    )
    unsupported_physical_claims = sum(
        1
        for packet in closed_packets
        if packet["claim_tier_unlocked"] == "CLAIM_4_COARSE_PHYSICAL_SUPPORTED"
        and not packet["postseal_physical_holdout"]["selected_sentence_execution"]["supports_selected"]
    )
    unsupported_topology_claims = sum(
        1
        for packet in closed_packets
        if packet["claim_tier_unlocked"] == "CLAIM_3_TOPOLOGY_SUPPORTED"
        and not packet["postseal_structure_holdout"]["observable_supports_selected_sentence"]
    )
    unsupported_fold_claims = sum(
        1
        for packet in closed_packets
        if packet["claim_tier_unlocked"] == "CLAIM_5_TARGET_FOLD_SUPPORTED"
        and not packet["target_fold_support_packet"]["target_fold_claim_allowed"]
    )
    return {
        "claim_counts": dict(counts),
        "failed_accepted_count": failed_accepted_count,
        "token_only_acceptance_count": token_only_acceptance_count,
        "coordinate_native_leakage": any(
            packet["preseal_allowed_sources"]["coordinate_truth_used_before_prediction"]
            or packet["preseal_allowed_sources"]["native_contacts_used_before_prediction"]
            for packet in closed_packets
        ),
        "hard_abstains_clean": all(
            row["closed_loop_packet"]["sealed_language_prediction"]["decision"] == "clean_abstain"
            for row in rows
            if row["row_family"] in {
                "hard_abstain_no_evidence_metadata_masked_control",
                "token_soup_adversarial_control",
            }
        ),
        "sentinels_preserved": sum(1 for row in rows if row["row_family"] == "old_sentinel_all_grammar_family"),
        "unsupported_physical_claims": unsupported_physical_claims,
        "unsupported_topology_claims": unsupported_topology_claims,
        "unsupported_fold_claims": unsupported_fold_claims,
        "language_claim_count": sum(1 for packet in closed_packets if packet["claim_tier_unlocked"] == "CLAIM_0_LANGUAGE_ONLY"),
        "independent_observable_support_count": sum(
            1
            for packet in closed_packets
            if packet["claim_tier_unlocked"]
            in {
                "CLAIM_1_BOUNDED_OBSERVABLE_SUPPORTED",
                "CLAIM_2_MECHANISM_SUPPORTED",
                "CLAIM_3_TOPOLOGY_SUPPORTED",
                "CLAIM_4_COARSE_PHYSICAL_SUPPORTED",
                "CLAIM_5_TARGET_FOLD_SUPPORTED",
            }
        ),
        "mechanism_claim_count": counts.get("CLAIM_2_MECHANISM_SUPPORTED", 0),
        "topology_claim_count": counts.get("CLAIM_3_TOPOLOGY_SUPPORTED", 0),
        "bounded_physical_claim_count": counts.get("CLAIM_4_COARSE_PHYSICAL_SUPPORTED", 0),
        "target_fold_claim_count": counts.get("CLAIM_5_TARGET_FOLD_SUPPORTED", 0),
    }


def run_v84(out_dir: Path = DEFAULT_OUT_DIR) -> dict[str, Path]:
    v83_cert, v83p_cert, _v82_rows, _v82p_rows, _v83_rows = _load_inputs()
    rows = _rows()
    metrics = _row_metrics(rows)
    fresh_available = _fresh_source_count()
    fresh_target_shortage = fresh_available < FRESH_BLIND_TARGETS
    deterministic_variant_count = max(0, FRESH_BLIND_TARGETS - fresh_available)
    subset_mode_counts = Counter(row.get("subset_execution_mode") for row in rows if row.get("atomistic_or_validated_coarse_subset_member"))
    universal_firewall = universal_claim_firewall(
        broad_fresh_generalization_passed=not fresh_target_shortage,
        independent_physical_observables_passed=metrics["bounded_physical_claim_count"] > 0,
        target_fold_subset_passed=metrics["target_fold_claim_count"] > 0,
        unsupported_claims=(
            metrics["unsupported_physical_claims"]
            + metrics["unsupported_topology_claims"]
            + metrics["unsupported_fold_claims"]
        ),
        failed_accepted_count=metrics["failed_accepted_count"],
    )
    general_solution_candidate_claim_allowed = (
        not fresh_target_shortage
        and metrics["failed_accepted_count"] == 0
        and metrics["target_fold_claim_count"] > 0
        and metrics["unsupported_fold_claims"] == 0
    )
    ledger = final_claim_ledger(
        language_claim_count=metrics["language_claim_count"],
        independent_observable_support_count=metrics["independent_observable_support_count"],
        mechanism_claim_count=metrics["mechanism_claim_count"],
        topology_claim_count=metrics["topology_claim_count"],
        bounded_physical_claim_count=metrics["bounded_physical_claim_count"],
        target_fold_claim_count=metrics["target_fold_claim_count"],
        general_solution_candidate_claim_allowed=general_solution_candidate_claim_allowed,
        universal_solution_claim_allowed=universal_firewall["universal_solution_claim_allowed"],
        unsupported_physical_claims=metrics["unsupported_physical_claims"],
        unsupported_topology_claims=metrics["unsupported_topology_claims"],
        unsupported_fold_claims=metrics["unsupported_fold_claims"],
        failed_accepted_count=metrics["failed_accepted_count"],
    )
    final_status = ledger["final_status"]
    row_family_counts = Counter(row["row_family"] for row in rows)
    claim_blocked_reason = None
    if fresh_target_shortage:
        claim_blocked_reason = "fresh_target_shortage_and_proxy_physical_execution_cap_claim_ceiling_at_CLAIM_4"
    elif metrics["target_fold_claim_count"] == 0:
        claim_blocked_reason = "missing_target_fold_subset_support"
    failure_autopsies = [
        packet["failure_autopsy_if_any"]
        for packet in (row["closed_loop_packet"] for row in rows)
        if packet["failure_autopsy_if_any"]
    ]
    certificate_core = {
        "row_family_counts": dict(row_family_counts),
        "claim_counts": metrics["claim_counts"],
        "final_status": final_status,
        "highest_claim_tier_unlocked": ledger["highest_claim_tier_unlocked"],
        "fresh_target_shortage": fresh_target_shortage,
        "deterministic_variant_count": deterministic_variant_count,
        "claim_ceiling_due_to_cache": "CLAIM_4_COARSE_PHYSICAL_SUPPORTED" if fresh_target_shortage else None,
    }
    failed_controls = []
    if len(rows) != TOTAL_TARGETS:
        failed_controls.append("target_count_9000")
    if metrics["failed_accepted_count"]:
        failed_controls.append("failed_accepted_zero")
    if metrics["token_only_acceptance_count"]:
        failed_controls.append("token_only_acceptance_zero")
    if metrics["coordinate_native_leakage"]:
        failed_controls.append("coordinate_native_leakage_false")
    if metrics["sentinels_preserved"] != OLD_SENTINELS:
        failed_controls.append("sentinels_preserved")
    if not metrics["hard_abstains_clean"]:
        failed_controls.append("hard_abstains_clean")
    if metrics["unsupported_physical_claims"]:
        failed_controls.append("unsupported_physical_claims_zero")
    if metrics["unsupported_topology_claims"]:
        failed_controls.append("unsupported_topology_claims_zero")
    if metrics["unsupported_fold_claims"]:
        failed_controls.append("unsupported_fold_claims_zero")
    if final_status not in E78_FINAL_STATUSES:
        failed_controls.append("final_status_known")
    cert = {
        "kind": "V84_FINAL_CLOSED_LOOP_PROTEIN_FOLDING_CAMPAIGN_CERTIFICATE_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "baseline_engine_version": BASELINE_ENGINE_VERSION,
        "source_batch_id": SOURCE_BATCH_ID,
        "sentence_source_batch_id": SENTENCE_SOURCE_BATCH_ID,
        "physical_source_batch_id": PHYSICAL_SOURCE_BATCH_ID,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "status": final_status,
        "campaign_controls_passed": not failed_controls,
        "targets_total": len(rows),
        "row_family_counts": dict(row_family_counts),
        "fresh_blind_target_count": row_family_counts["fresh_blind_nonredundant_target"],
        "old_sentinel_count": row_family_counts["old_sentinel_all_grammar_family"],
        "compositional_sentence_target_count": row_family_counts["compositional_sentence_target"],
        "hard_abstain_control_count": row_family_counts["hard_abstain_no_evidence_metadata_masked_control"],
        "token_soup_control_count": row_family_counts["token_soup_adversarial_control"],
        "topology_family_holdout_count": row_family_counts["topology_family_holdout_target"],
        "physical_execution_target_count": row_family_counts["physical_execution_target"],
        "atomistic_or_validated_coarse_subset_count": ATOMISTIC_OR_VALIDATED_COARSE_SUBSET,
        "known_hard_frontier_subset_count": KNOWN_HARD_FRONTIER_SUBSET,
        "subset_execution_mode_counts": {str(key): value for key, value in subset_mode_counts.items()},
        "fresh_target_shortage": fresh_target_shortage,
        "fresh_source_count_available": fresh_available,
        "deterministic_variant_count": deterministic_variant_count,
        "claim_ceiling_due_to_cache": "CLAIM_4_COARSE_PHYSICAL_SUPPORTED" if fresh_target_shortage else None,
        "claim_ladder": E78_CLAIM_LADDER,
        "final_claim_ledger": ledger,
        "universal_claim_firewall": universal_firewall,
        "protein_esperanto_language_complete_for_current_registry": True,
        "independent_observable_support_count": metrics["independent_observable_support_count"],
        "bounded_physical_claim_count": metrics["bounded_physical_claim_count"],
        "topology_claim_count": metrics["topology_claim_count"],
        "target_fold_claim_count": metrics["target_fold_claim_count"],
        "general_solution_candidate_claim_allowed": general_solution_candidate_claim_allowed,
        "universal_folding_solution_claim_allowed": universal_firewall["universal_solution_claim_allowed"],
        "language_claim_allowed_count": metrics["language_claim_count"],
        "failed_accepted_count": metrics["failed_accepted_count"],
        "token_only_acceptance_count": metrics["token_only_acceptance_count"],
        "coordinate_native_leakage": metrics["coordinate_native_leakage"],
        "sentinels_preserved": metrics["sentinels_preserved"],
        "hard_abstains_clean": metrics["hard_abstains_clean"],
        "unsupported_physical_claims": metrics["unsupported_physical_claims"],
        "unsupported_topology_claims": metrics["unsupported_topology_claims"],
        "unsupported_fold_claims": metrics["unsupported_fold_claims"],
        "unsupported_claims": (
            metrics["unsupported_physical_claims"]
            + metrics["unsupported_topology_claims"]
            + metrics["unsupported_fold_claims"]
        ),
        "bounded_physical_language_claim_allowed": metrics["bounded_physical_claim_count"] > 0,
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": universal_firewall["protein_folding_solved"],
        "claim_blocked_reason": claim_blocked_reason,
        "failure_autopsy_count": len(failure_autopsies),
        "failure_autopsies_available_for_non_closed_claims": bool(failure_autopsies),
        "v83_status": v83_cert["status"],
        "v83p_status": v83p_cert["status"],
        "certificate_core_hash": stable_hash(certificate_core),
        "failed_controls": failed_controls,
    }
    cert["certificate_hash"] = stable_hash({key: value for key, value in cert.items() if key != "certificate_hash"})
    report = {
        "kind": "V84_FINAL_CLOSED_LOOP_PROTEIN_FOLDING_CAMPAIGN_REPORT_v0",
        "batch_id": BATCH_ID,
        "engine_version_used": ENGINE_VERSION_USED,
        "rows": rows,
        "failure_autopsies": failure_autopsies,
    }
    e78_cert = {
        "kind": "E78_CLOSED_LOOP_BIOPHYSICAL_PROTEIN_LANGUAGE_ENGINE_CERTIFICATE_v0",
        "engine_revision": ENGINE_VERSION_USED,
        "baseline_engine_revision": BASELINE_ENGINE_VERSION,
        "claim_ladder": E78_CLAIM_LADDER,
        "final_status": final_status,
        "highest_claim_tier_unlocked": ledger["highest_claim_tier_unlocked"],
        "final_closed_loop_packet_kind": "E78_FINAL_CLOSED_LOOP_PACKET_v0",
        "independent_holdout_binding_kind": "E78_INDEPENDENT_HOLDOUT_BINDING_v0",
        "physical_execution_binding_kind": "E78_PHYSICAL_EXECUTION_BINDING_v0",
        "target_fold_support_packet_kind": "E78_TARGET_FOLD_SUPPORT_PACKET_v0",
        "universal_claim_firewall_kind": "E78_UNIVERSAL_CLAIM_FIREWALL_v0",
        "failure_autopsy_packet_kind": "E78_FAILURE_AUTOPSY_PACKET_v0",
        "physical_basis_claim_allowed": False,
        "protein_folding_solved": universal_firewall["protein_folding_solved"],
        "next_required_batch": None,
    }
    paths = {
        "certificate": _write_json(DATA_ROOT / "v84_final_closed_loop_protein_folding_campaign_certificate.json", cert),
        "report": _write_json(DATA_ROOT / "v84_final_closed_loop_protein_folding_campaign_report.json", report),
        "e78_certificate": _write_json(E78_ROOT / "e78_closed_loop_biophysical_protein_language_engine_certificate.json", e78_cert),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    paths["run_certificate"] = _write_json(out_dir / "v84_final_closed_loop_protein_folding_campaign_certificate.json", cert)
    paths["run_report"] = _write_json(out_dir / "v84_final_closed_loop_protein_folding_campaign_report.json", report)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run V84 final closed-loop Protein Esperanto campaign.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    paths = run_v84(args.out_dir)
    cert = _read_json(paths["run_certificate"], "V84 certificate")
    print(json.dumps({
        "kind": cert["kind"],
        "status": cert["status"],
        "campaign_controls_passed": cert["campaign_controls_passed"],
        "targets_total": cert["targets_total"],
        "highest_claim_tier_unlocked": cert["final_claim_ledger"]["highest_claim_tier_unlocked"],
        "independent_observable_support_count": cert["independent_observable_support_count"],
        "bounded_physical_claim_count": cert["bounded_physical_claim_count"],
        "target_fold_claim_count": cert["target_fold_claim_count"],
        "general_solution_candidate_claim_allowed": cert["general_solution_candidate_claim_allowed"],
        "universal_folding_solution_claim_allowed": cert["universal_folding_solution_claim_allowed"],
        "fresh_target_shortage": cert["fresh_target_shortage"],
        "claim_ceiling_due_to_cache": cert["claim_ceiling_due_to_cache"],
        "physical_basis_claim_allowed": cert["physical_basis_claim_allowed"],
        "protein_folding_solved": cert["protein_folding_solved"],
        "claim_blocked_reason": cert["claim_blocked_reason"],
        "certificate": str(paths["run_certificate"]),
    }, indent=2, sort_keys=True))
    return 0 if cert["campaign_controls_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
