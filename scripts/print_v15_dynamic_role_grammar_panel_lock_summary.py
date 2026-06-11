#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCKED" / "v15_dynamic_role_grammar_panel_locked_certificate.json"


def main() -> None:
    if not CERT.exists():
        raise SystemExit(f"missing V15 lock certificate: {CERT}\nrun scripts/run_v15_dynamic_role_grammar_panel_lock.sh first")
    p = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V15 DYNAMIC ROLE GRAMMAR PANEL LOCK ===")
    print("lock_status:", p.get("lock_status"))
    print("lock_failed_checks:", p.get("lock_failed_checks"))
    print("source_global_status:", p.get("source_global_status"))
    print("positive_evidence_proteins:", p.get("positive_evidence_proteins"))
    print("missing_artifacts:", p.get("missing_artifacts"))
    print("claim_allowed:", p.get("claim_allowed"))
    print("locked_claim:", p.get("locked_claim"))
    print("locked_interpretation:", p.get("locked_interpretation"))
    print("lock_checks:", p.get("lock_checks"))


if __name__ == "__main__":
    main()
