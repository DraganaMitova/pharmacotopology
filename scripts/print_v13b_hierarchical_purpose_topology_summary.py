#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

root = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1]))
cert_path = (
    root
    / "first_contact_clean_pharmacotopology_layer_run"
    / "V13b_1CLL_HIERARCHICAL_PURPOSE_TOPOLOGY_READOUT"
    / "v13b_1cll_hierarchical_purpose_topology_readout_certificate.json"
)
if not cert_path.exists():
    raise SystemExit(f"missing certificate: {cert_path}")

o = json.loads(cert_path.read_text(encoding="utf-8"))
band = o.get("selected_frequency_band") or {}
lock = o.get("claim_lock_check") or {}
print("=== V13b 1CLL HIERARCHICAL PURPOSE TOPOLOGY READOUT ===")
print("run_mode:", o.get("run_mode"))
print("trajectory_count:", o.get("trajectory_count"))
print("target_role:", o.get("target_role"))
print("topology_policy:", o.get("topology_policy"))
print("domain_boundaries:", o.get("domain_boundaries"))
print("domain_roles:", o.get("domain_roles"))
print("input_preflight_status:", o.get("input_preflight_status"))
legacy = o.get("legacy_v5_preflight") or {}
legacy_req = legacy.get("v5_requirements") if isinstance(legacy, dict) else {}
print("legacy_v5_failed_checks:", legacy_req.get("failed_checks") if isinstance(legacy_req, dict) else None)
print("legacy_runtime_v10_selected_pair_count:", o.get("legacy_runtime_v10_selected_pair_count"))
print("hierarchical_topology_decision:", o.get("hierarchical_topology_decision"))
print("selected_threshold:", o.get("selected_threshold"))
print("selected_pair_count:", o.get("selected_pair_count"))
print("selected_N_domain_core:", o.get("selected_N_domain_core"))
print("selected_C_domain_core:", o.get("selected_C_domain_core"))
print("selected_interdomain_hinge:", o.get("selected_interdomain_hinge"))
print("selected_local_support:", o.get("selected_local_support"))
print("selected_medium_support:", o.get("selected_medium_support"))
print("support_by_selected_pair:", band.get("support_by_selected_pair", {}))
print("mean_frequency_by_selected_pair:", band.get("mean_frequency_by_selected_pair", {}))
print("chemical_score_by_selected_pair:", band.get("chemical_score_by_selected_pair", {}))
print("dca_score_by_selected_pair:", band.get("dca_score_by_selected_pair", {}))
print("noise_added:", band.get("noise_added"))
print("long_range_evidence_polluted:", band.get("long_range_evidence_polluted"))
print("classification_coverage_ratio:", band.get("classification_coverage_ratio"))
print("dca_absolute_support_pass:", band.get("dca_absolute_support_pass"))
print("dca_background_enrichment_ratio:", band.get("dca_background_enrichment_ratio"))
print("dca_background_enrichment_pass:", band.get("dca_background_enrichment_pass"))
print("dca_pass_semantics:", band.get("dca_pass_semantics"))
print("claim_lock_status:", o.get("claim_lock_status"))
print("claim_lock_failed_checks:", o.get("claim_lock_failed_checks"))
print("claim_allowed:", o.get("claim_allowed"))

print("\nCandidate pair roles:")
for key, payload in sorted((o.get("role_by_candidate_pair") or {}).items()):
    print(f"{key}:")
    for field in [
        "sequence_separation", "domain_relation", "topology_role", "evidence_class",
        "dca_score", "chemical_score", "tail_frequency_mean", "tail_frequency_min",
        "tail_frequency_max", "tail_presence_count_at_0_50", "inside_effective_balanced",
    ]:
        print(f"  {field}: {payload.get(field)}")
