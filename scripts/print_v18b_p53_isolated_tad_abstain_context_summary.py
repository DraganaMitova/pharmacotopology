#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT" / "v18b_p53_isolated_tad_abstain_context_certificate.json"

def main() -> None:
    data = json.loads(CERT.read_text(encoding="utf-8"))
    print("=== V18b P53 ISOLATED-TAD ABSTAIN CONTEXT ===")
    for key in [
        "run_mode",
        "test_status",
        "target_id",
        "role_class",
        "positive_pressure_evidence_found",
        "positive_folding_evidence_found",
        "bound_complex_partner_induced_evidence_preserved",
        "partner_induced_role_evidence_found",
        "isolated_TAD_material_available",
        "isolated_TAD_status",
        "isolated_TAD_autonomous_core_selected",
        "isolated_TAD_fold_claim_made",
        "isolated_TAD_clean_abstain_valid",
        "selected_core_or_clean_abstain",
        "contrast_status",
        "claim_allowed",
        "new_md_executed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
    ]:
        print(f"{key}: {data.get(key)}")
    print("forbidden_misclassification_violations:", data.get("forbidden_misclassification_violations"))
    print("available_evidence:", data.get("available_evidence"))
    print("missing_evidence:", data.get("missing_evidence"))
    print("failed_checks:", data.get("failed_checks"))
    print("next_recommended_panel:", data.get("next_recommended_panel"))

if __name__ == "__main__":
    main()
