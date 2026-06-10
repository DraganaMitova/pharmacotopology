#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V15_4AKE_BALANCED_CANDIDATE_READOUT" / "v15_4ake_balanced_candidate_readout_certificate.json"


def main() -> None:
    if not CERT.exists():
        raise SystemExit(f"missing 4AKE balanced-candidate certificate: {CERT}\nrun scripts/run_v15_4ake_balanced_candidate_readout.sh first")
    payload = json.loads(CERT.read_text(encoding="utf-8"))
    band = payload.get("selected_frequency_band") or {}
    print("=== V15 4AKE BALANCED CANDIDATE READOUT ===")
    print("run_mode:", payload.get("run_mode"))
    print("artifact_status:", payload.get("artifact_status"))
    print("target_role:", payload.get("target_role"))
    print("grammar_policy:", payload.get("grammar_policy"))
    print("trajectory_count:", payload.get("trajectory_count"))
    print("effective_balanced_count:", payload.get("effective_balanced_count"))
    print("selected_threshold:", band.get("threshold"))
    print("selected_pairs:", payload.get("selected_pairs"))
    print("selected_domain_core:", payload.get("selected_domain_core"))
    print("selected_hinge_or_interdomain:", payload.get("selected_hinge_or_interdomain"))
    print("support_by_selected_pair:", payload.get("support_by_selected_pair"))
    print("mean_frequency_by_selected_pair:", payload.get("mean_frequency_by_selected_pair"))
    print("chemical_score_by_selected_pair:", payload.get("chemical_score_by_selected_pair"))
    print("dca_score_by_selected_pair:", payload.get("dca_score_by_selected_pair"))
    print("noise_added:", payload.get("noise_added"))
    print("long_range_evidence_polluted:", payload.get("long_range_evidence_polluted"))
    print("classification_coverage_ratio:", payload.get("classification_coverage_ratio"))
    print("positive_evidence_found:", payload.get("positive_evidence_found"))
    print("claim_lock_status:", payload.get("claim_lock_status"))
    print("claim_lock_failed_checks:", payload.get("claim_lock_failed_checks"))
    print("final_status:", payload.get("final_status"))
    print("claim_allowed:", payload.get("claim_allowed"))
    print("\nCandidate pairs:")
    for key, role in (payload.get("candidate_pair_roles") or {}).items():
        print(f"{key}:")
        print("  role_decision:", role.get("role_decision"))
        print("  domain_relation:", role.get("domain_relation"))
        print("  tail_frequency_mean:", role.get("tail_frequency_mean"))
        print("  tail_frequency_min:", role.get("tail_frequency_min"))
        print("  tail_frequency_max:", role.get("tail_frequency_max"))
        print("  tail_presence_count_at_0_50:", role.get("tail_presence_count_at_0_50"))
        print("  chemical_score:", role.get("chemical_score"))
        print("  dca_score:", role.get("dca_score"))
        print("  fixed_residue_cutoff_used:", role.get("fixed_residue_cutoff_used"))


if __name__ == "__main__":
    main()
