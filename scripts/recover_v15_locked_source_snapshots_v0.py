#!/usr/bin/env python3
from __future__ import annotations

"""Recover the small V15 source certificates from the locked conversation outputs.

This is NOT a raw-trajectory restore.  It writes small, claim-disabled summary
certificates that were already observed in the locked V15/V20+ run logs so the
postprocess readouts can be rebuilt after the generated runtime folder was
accidentally deleted/overwritten.

Heavy MD trajectories are not reconstructed here.  If a future workflow needs
raw trajectories, rerun/restore those explicitly.
"""

import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
SNAPSHOT_DIR = REPO_ROOT / "data" / "locked_runtime_certificates" / "V15_RECOVERED_SOURCE_SNAPSHOTS"

RECOVERY_NOTE = {
    "recovery_source": "recovered_from_prior_locked_outputs_shared_in_chat_transcript",
    "recovery_scope": "small_claim_disabled_summary_certificate_not_raw_trajectory_restore",
    "raw_md_trajectory_restored": False,
    "claim_allowed": False,
}


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def pair_role(i: int, j: int, role: str, relation: str, mean: float, tmin: float, tmax: float, chem: float, dca: float, selected: bool=False) -> dict[str, Any]:
    return {
        "sequence_separation": j - i,
        "normalized_sequence_separation": round((j - i) / 213.0, 6),
        "domain_relation": relation,
        "role_decision": role,
        "evidence_class": (
            "D4_domain_core_evidence" if "D4" in role else
            "D2_domain_core_evidence" if "D2" in role else
            "domain_hinge_or_interdomain_closure_monitor_evidence"
        ),
        "selected": selected,
        "tail_frequency_mean": mean,
        "tail_frequency_min": tmin,
        "tail_frequency_max": tmax,
        "tail_presence_count_at_0_50": 10,
        "chemical_score": chem,
        "dca_score": dca,
        "separation_filter_applied": False,
        "fixed_residue_cutoff_used": False,
    }


def build_1ubq() -> dict[str, Any]:
    cert = {
        "kind": "V13a_1UBQ_ADAPTIVE_CHEMICAL_POLICY_READOUT_v0_RECOVERED_SNAPSHOT",
        "target_id": "1UBQ",
        "target_role": "single_domain_compact",
        "sequence_length": 76,
        "chemical_policy": "adaptive_soft_guard",
        "claim_lock_status": "claim_locked_pending_cross_target_validation",
        "claim_lock_failed_checks": ["dca_background_enrichment_pass"],
        "claim_allowed": False,
        "new_md_executed": False,
        "selected_frequency_band": {
            "selected_balanced_core": [[23, 48]],
            "support_by_selected_pair": {"23-48": 9},
            "mean_frequency_by_selected_pair": {"23-48": 0.784444},
            "chemical_score_by_selected_pair": {"23-48": 0.1},
            "dca_score_by_selected_pair": {"23-48": 0.947233},
            "dca_mean_selected": 0.947233,
            "dca_background_enrichment_ratio": 0.969389,
            "dca_background_enrichment_pass": False,
            "dca_absolute_support_pass": True,
            "noise_added": 0,
            "long_range_evidence_polluted": False,
            "classification_coverage_ratio": 1.0,
        },
        "positive_evidence_found": True,
        "final_status": "single_domain_compact_signal_found_with_dynamic_separation_context;claim_allowed=false",
        "recovery_metadata": RECOVERY_NOTE,
    }
    return cert


def build_1cll() -> dict[str, Any]:
    role_by_candidate_pair = {
        "22-29": {
            "domain_relation": "intradomain_N",
            "topology_role": "local_support",
            "evidence_class": "local_support",
            "mean_frequency": 0.511,
            "tail_frequency_min": 0.44,
            "tail_frequency_max": 0.66,
            "chemical_score": 0.1,
            "dca_score": 0.72114,
        },
        "24-60": {
            "domain_relation": "intradomain_N",
            "topology_role": "N_domain_compact_core_candidate",
            "evidence_class": "N_domain_core_evidence",
            "mean_frequency": 0.735,
            "tail_frequency_min": 0.65,
            "tail_frequency_max": 0.82,
            "chemical_score": 0.1,
            "dca_score": 1.0,
        },
        "95-102": {
            "domain_relation": "intradomain_C",
            "topology_role": "local_support",
            "evidence_class": "local_support",
            "mean_frequency": 0.492,
            "tail_frequency_min": 0.4,
            "tail_frequency_max": 0.56,
            "chemical_score": 0.1,
            "dca_score": 0.72114,
        },
        "97-133": {
            "domain_relation": "intradomain_C",
            "topology_role": "C_domain_compact_core_candidate",
            "evidence_class": "C_domain_core_evidence",
            "mean_frequency": 0.742,
            "tail_frequency_min": 0.66,
            "tail_frequency_max": 0.82,
            "chemical_score": 0.1,
            "dca_score": 1.0,
        },
    }
    return {
        "kind": "V13b_1CLL_HIERARCHICAL_PURPOSE_TOPOLOGY_READOUT_v0_RECOVERED_SNAPSHOT",
        "target_id": "1CLL",
        "target_role": "multi_domain_composite",
        "sequence_length": 148,
        "topology_policy": "hierarchical_domain_core_plus_interdomain",
        "chemical_policy": "adaptive_soft_guard",
        "selected_N_domain_core": [],
        "selected_C_domain_core": [[97, 133]],
        "selected_interdomain_hinge": [],
        "selected_local_support": [],
        "selected_medium_support": [],
        "role_by_candidate_pair": role_by_candidate_pair,
        "selected_frequency_band": {
            "selected_C_domain_core": [[97, 133]],
            "selected_N_domain_core": [],
            "selected_interdomain_hinge": [],
            "support_by_selected_pair": {"97-133": 7},
            "mean_frequency_by_selected_pair": {"97-133": 0.742},
            "chemical_score_by_selected_pair": {"97-133": 0.1},
            "dca_score_by_selected_pair": {"97-133": 1.0},
            "dca_background_enrichment_ratio": 1.162021,
            "dca_background_enrichment_pass": True,
            "noise_added": 0,
            "long_range_evidence_polluted": False,
            "classification_coverage_ratio": 1.0,
        },
        "claim_lock_status": "claim_lock_passed_but_claim_still_disabled",
        "claim_lock_failed_checks": [],
        "claim_allowed": False,
        "new_md_executed": False,
        "positive_evidence_found": True,
        "final_status": "multi_domain_composite_domain_core_signal_found_with_dynamic_separation_context;interdomain_hinge_not_yet_proven;claim_allowed=false",
        "recovery_metadata": RECOVERY_NOTE,
    }


def build_4ake_balanced() -> dict[str, Any]:
    candidate_roles = {
        "113-136": pair_role(113, 136, "domain_hinge_or_interdomain_closure_candidate", "interdomain_D3_D4", 0.746, 0.68, 0.82, 0.9, 0.738199),
        "124-135": pair_role(124, 135, "D4_domain_compact_core_candidate", "intradomain_D4", 0.768, 0.74, 0.80, 0.1, 0.866565, selected=True),
        "126-146": pair_role(126, 146, "D4_domain_compact_core_candidate", "intradomain_D4", 0.721, 0.67, 0.80, 0.9, 1.0),
        "126-149": pair_role(126, 149, "D4_domain_compact_core_candidate", "intradomain_D4", 0.732, 0.65, 0.78, 0.1, 0.719678),
        "127-154": pair_role(127, 154, "D4_domain_compact_core_candidate", "intradomain_D4", 0.717, 0.67, 0.78, 0.3, 0.882066),
        "129-146": pair_role(129, 146, "D4_domain_compact_core_candidate", "intradomain_D4", 0.747, 0.70, 0.80, 0.1, 0.750517),
        "135-155": pair_role(135, 155, "D4_domain_compact_core_candidate", "intradomain_D4", 0.719, 0.67, 0.78, 0.3, 0.757759),
        "145-152": pair_role(145, 152, "D4_domain_compact_core_candidate", "intradomain_D4", 0.764, 0.68, 0.87, 0.9, 0.577392),
        "35-64": pair_role(35, 64, "D2_domain_compact_core_candidate", "intradomain_D2", 0.723, 0.68, 0.77, 1.0, 0.824903),
        "52-60": pair_role(52, 60, "D2_domain_compact_core_candidate", "intradomain_D2", 0.682, 0.62, 0.77, 0.3, 0.756135),
        "65-86": pair_role(65, 86, "domain_hinge_or_interdomain_closure_candidate", "interdomain_D2_D3", 0.657, 0.59, 0.72, 1.0, 0.727478),
    }
    return {
        "kind": "V15_4AKE_BALANCED_CANDIDATE_READOUT_v0_RECOVERED_SNAPSHOT",
        "run_mode": "postprocess_only_no_new_simulation_recovered_snapshot",
        "artifact_status": "present_machine_readable_4ake_balanced_candidate_readout_positive",
        "target_id": "4AKE",
        "target_role": "domain_hinge_closure_object",
        "grammar_policy": "domain_hinge_dynamic_balanced_candidate_frequency_readout",
        "trajectory_count": 10,
        "effective_balanced_count": 11,
        "selected_threshold": 0.75,
        "selected_pairs": ["124-135"],
        "selected_domain_core": ["124-135"],
        "selected_hinge_or_interdomain": [],
        "selected_balanced_core": ["124-135"],
        "support_by_selected_pair": {"124-135": 7},
        "mean_frequency_by_selected_pair": {"124-135": 0.768},
        "chemical_score_by_selected_pair": {"124-135": 0.1},
        "dca_score_by_selected_pair": {"124-135": 0.866565},
        "candidate_pair_roles": candidate_roles,
        "noise_added": 0,
        "long_range_evidence_polluted": False,
        "classification_coverage_ratio": 1.0,
        "positive_evidence_found": True,
        "claim_lock_status": "claim_locked_4ake_candidate_readout_claim_disabled",
        "claim_lock_failed_checks": [],
        "claim_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "final_status": "domain_hinge_object_balanced_candidate_signal_found_under_dynamic_frequency_readout;claim_allowed=false",
        "recovery_metadata": RECOVERY_NOTE,
    }


def main() -> None:
    outputs = {
        RUN_ROOT / "V13a_1UBQ_ADAPTIVE_CHEMICAL_POLICY_READOUT" / "v13a_1ubq_purpose_gate_readout_certificate.json": build_1ubq(),
        RUN_ROOT / "V13b_1CLL_HIERARCHICAL_PURPOSE_TOPOLOGY_READOUT" / "v13b_1cll_hierarchical_purpose_topology_readout_certificate.json": build_1cll(),
        RUN_ROOT / "V15_4AKE_BALANCED_CANDIDATE_READOUT" / "v15_4ake_balanced_candidate_readout_certificate.json": build_4ake_balanced(),
    }
    for path, payload in outputs.items():
        write_json(path, payload)

    # Export a copy to data/locked_runtime_certificates so runtime output can be restored later.
    snapshot_outputs = {}
    for path, payload in outputs.items():
        rel = path.relative_to(RUN_ROOT)
        snapshot_path = SNAPSHOT_DIR / rel
        write_json(snapshot_path, payload)
        snapshot_outputs[str(rel)] = str(snapshot_path.relative_to(REPO_ROOT))

    manifest = {
        "kind": "V15_RECOVERED_SOURCE_SNAPSHOTS_MANIFEST_v0",
        "claim_allowed": False,
        "recovery_scope": "small_locked_summary_certificates_not_raw_trajectories",
        "raw_md_trajectory_restored": False,
        "snapshot_outputs": snapshot_outputs,
    }
    write_json(SNAPSHOT_DIR / "v15_recovered_source_snapshots_manifest.json", manifest)

    print(json.dumps({
        "status": "V15_RECOVERED_SOURCE_SNAPSHOTS_WRITTEN",
        "written_runtime_certificates": [str(p) for p in outputs],
        "snapshot_manifest": str(SNAPSHOT_DIR / "v15_recovered_source_snapshots_manifest.json"),
        "claim_allowed": False,
        "raw_md_trajectory_restored": False,
    }, indent=2))


if __name__ == "__main__":
    main()
