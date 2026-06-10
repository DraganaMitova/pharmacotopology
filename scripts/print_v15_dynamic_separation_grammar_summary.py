#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = (
    REPO_ROOT
    / "first_contact_clean_pharmacotopology_layer_run"
    / "V15_DYNAMIC_SEPARATION_GRAMMAR_READOUT"
    / "v15_dynamic_separation_grammar_readout_certificate.json"
)


def main() -> None:
    if not CERT.exists():
        raise SystemExit(f"missing V15 certificate: {CERT}\nrun scripts/run_v15_dynamic_separation_grammar_readout.sh first")
    payload = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V15 DYNAMIC SEPARATION GRAMMAR READOUT ===")
    print("run_mode:", payload.get("run_mode"))
    print("global_status:", payload.get("global_status"))
    print("separation_policy:", payload.get("separation_policy"))
    print("fixed_residue_cutoff_used:", payload.get("fixed_residue_cutoff_used"))
    print("claim_allowed:", payload.get("claim_allowed"))
    print("positive_evidence_proteins:", payload.get("positive_evidence_proteins"))
    print("missing_artifacts:", payload.get("missing_artifacts"))
    checks = payload.get("coherence_checks", {})
    print("no_fixed_residue_cutoff_used_anywhere:", checks.get("no_fixed_residue_cutoff_used_anywhere"))
    print("no_claim_allowed_anywhere:", checks.get("no_claim_allowed_anywhere"))
    print("\nGrammar axes:")
    for axis in payload.get("grammar_axes", []):
        print("  -", axis)
    print("\nProtein rows:")
    for row in payload.get("protein_rows", []):
        print(f"\n{row.get('protein')}:")
        for key in [
            "artifact_status",
            "target_role",
            "grammar_policy",
            "separation_policy",
            "fixed_residue_cutoff_used",
            "selected_pairs",
            "replica_support",
            "chemical_policy",
            "topology_policy",
            "noise_added",
            "long_range_evidence_polluted",
            "classification_coverage_ratio",
            "claim_lock_status",
            "claim_lock_failed_checks",
            "final_status",
            "claim_allowed",
        ]:
            print(f"  {key}:", row.get(key))
        roles = row.get("dynamic_pair_roles") or {}
        if roles:
            print("  dynamic_pair_roles:")
            for pair, role in sorted(roles.items()):
                print(f"    {pair}:")
                for role_key in [
                    "sequence_separation",
                    "normalized_sequence_separation",
                    "domain_relation",
                    "role_decision",
                    "evidence_class",
                    "selected",
                    "support",
                    "mean_frequency",
                    "chemical_score",
                    "dca_score",
                    "separation_filter_applied",
                    "fixed_residue_cutoff_used",
                ]:
                    print(f"      {role_key}:", role.get(role_key))


if __name__ == "__main__":
    main()
