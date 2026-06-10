#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path

root = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1]))
cert = root / "first_contact_clean_pharmacotopology_layer_run" / "V13a_1UBQ_PURPOSE_GATE_READOUT" / "v13a_1ubq_purpose_gate_readout_certificate.json"
if not cert.exists():
    raise SystemExit(f"missing certificate: {cert}\nrun: bash $REPO_ROOT/scripts/run_v13a_1ubq_purpose_gate_readout.sh")
o = json.loads(cert.read_text(encoding="utf-8"))
sel = o.get("selected_frequency_band") or {}
pre = o.get("v13_transfer_preflight", {})
print("=== V13a 1UBQ PURPOSE GATE READOUT ===")
print("run_mode:", o.get("run_mode"))
print("trajectory_count:", o.get("trajectory_count"))
print("target_role:", o.get("target_purpose", {}).get("target_role"))
print("v13_preflight_status:", pre.get("status"))
print("v13_failed_checks:", pre.get("failed_checks"))
print("legacy_v5_failed_checks:", o.get("legacy_v5_preflight", {}).get("v5_requirements", {}).get("failed_checks"))
print("legacy_runtime_v10_selected_pair_count:", o.get("legacy_runtime_v10_selected_pair_count"))
print("purpose_gate_decision:", o.get("purpose_gate_decision"))
print("selected_threshold:", sel.get("threshold"))
print("selected_pair_count:", sel.get("selected_pair_count", 0))
print("selected_balanced_core:", sel.get("selected_balanced_core", []))
print("noise_added:", sel.get("noise_added"))
print("long_range_evidence_polluted:", sel.get("long_range_evidence_polluted"))
print("classification_coverage_ratio:", sel.get("classification_coverage_ratio"))
print("dca_mean_selected:", sel.get("dca_mean_selected"))
print("claim_allowed:", o.get("claim_allowed"))
for pair, payload in (sel.get("audit_pairs") or {}).items():
    print(f"audit_pair {pair}:")
    print("  selected_at_threshold:", payload.get("selected_at_threshold"))
    print("  support_at_threshold:", payload.get("support_at_threshold"))
    print("  tail_frequency_mean:", payload.get("tail_frequency_mean"))
    print("  dca_score:", payload.get("dca_score"))
    print("  blocked_reason_at_threshold:", payload.get("blocked_reason_at_threshold"))
