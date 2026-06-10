#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

root = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1]))
cert_path = (
    root
    / "first_contact_clean_pharmacotopology_layer_run"
    / "V14_UNIFIED_PROTEIN_ESPERANTO_GRAMMAR_READOUT"
    / "v14_unified_protein_esperanto_grammar_readout_certificate.json"
)
if not cert_path.exists():
    raise SystemExit(f"missing certificate: {cert_path}")

o = json.loads(cert_path.read_text(encoding="utf-8"))
checks = o.get("coherence_checks") or {}
print("=== V14 UNIFIED PROTEIN ESPERANTO GRAMMAR READOUT ===")
print("run_mode:", o.get("run_mode"))
print("global_status:", o.get("global_status"))
print("claim_allowed:", o.get("claim_allowed"))
print("positive_evidence_proteins:", checks.get("positive_evidence_proteins"))
print("missing_artifacts:", checks.get("missing_artifacts"))
print("no_claim_allowed_anywhere:", checks.get("no_claim_allowed_anywhere"))
print("all_rows_have_core_grammar_axes:", checks.get("all_rows_have_core_grammar_axes"))
print("grammar_axes:")
for axis in o.get("unified_grammar_axes", []):
    print("  -", axis)

print("\nProtein rows:")
for row in o.get("protein_rows", []):
    print(f"\n{row.get('protein')}:")
    print("  artifact_status:", row.get("artifact_status"))
    print("  target_role:", row.get("target_role"))
    print("  grammar_policy:", row.get("grammar_policy"))
    print("  selected_domain_core:", row.get("selected_domain_core"))
    print("  selected_balanced_core:", row.get("selected_balanced_core"))
    print("  selected_hinge_or_interdomain:", row.get("selected_hinge_or_interdomain"))
    print("  selected_local_support:", row.get("selected_local_support"))
    print("  replica_support:", row.get("replica_support"))
    print("  chemical_policy:", row.get("chemical_policy"))
    print("  topology_policy:", row.get("topology_policy"))
    print("  noise_added:", row.get("noise_added"))
    print("  long_range_evidence_polluted:", row.get("long_range_evidence_polluted"))
    print("  classification_coverage_ratio:", row.get("classification_coverage_ratio"))
    print("  claim_lock_status:", row.get("claim_lock_status"))
    print("  claim_lock_failed_checks:", row.get("claim_lock_failed_checks"))
    print("  final_status:", row.get("final_status"))
    print("  claim_allowed:", row.get("claim_allowed"))
