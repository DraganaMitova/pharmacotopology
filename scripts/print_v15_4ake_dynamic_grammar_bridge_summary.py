#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V15_4AKE_DYNAMIC_GRAMMAR_BRIDGE" / "v15_4ake_dynamic_grammar_bridge_certificate.json"


def main() -> None:
    if not CERT.exists():
        raise SystemExit(f"missing 4AKE bridge certificate: {CERT}\nrun scripts/run_v15_4ake_dynamic_grammar_bridge.sh first")
    payload = json.loads(CERT.read_text(encoding="utf-8"))
    row = payload.get("protein_row", {})
    print("=== V15 4AKE DYNAMIC GRAMMAR BRIDGE ===")
    print("run_mode:", payload.get("run_mode"))
    print("bridge_policy:", payload.get("bridge_policy"))
    print("artifact_status:", row.get("artifact_status"))
    print("target_role:", row.get("target_role"))
    print("grammar_policy:", row.get("grammar_policy"))
    print("selected_pairs:", row.get("selected_pairs"))
    print("positive_evidence_found:", row.get("positive_evidence_found"))
    print("claim_lock_status:", row.get("claim_lock_status"))
    print("claim_lock_failed_checks:", row.get("claim_lock_failed_checks"))
    print("final_status:", row.get("final_status"))
    print("claim_allowed:", row.get("claim_allowed"))
    print("legacy_visual_files_count:", len(row.get("legacy_visual_files") or []))
    print("legacy_global_certificates_present:", row.get("legacy_global_certificates_present"))
    print("input_material_present:", row.get("input_material_present"))


if __name__ == "__main__":
    main()
