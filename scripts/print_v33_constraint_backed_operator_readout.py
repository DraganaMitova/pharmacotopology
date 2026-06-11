#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
CERT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "V33_CONSTRAINT_BACKED_OPERATOR_READOUT" / "v33_constraint_backed_operator_readout_certificate.json"


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def main() -> int:
    if not CERT.exists():
        raise SystemExit(f"missing V33 certificate: {CERT}\nRun: python3 scripts/run_v33_constraint_backed_operator_readout_v0.py")
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    readout = cert.get("constraint_backed_operator_readout") or {}
    print("=== V33 CONSTRAINT-BACKED OPERATOR READOUT ===")
    for key in [
        "run_mode",
        "source_v32_status",
        "readout_status",
        "selected_V33_target",
        "selected_V33_panel",
        "claim_allowed",
        "new_md_executed",
        "membrane_md_executed",
        "new_MD_allowed",
        "new_MD_recommended",
        "fixed_threshold_policy",
        "target_specific_threshold_tuning_allowed",
        "fixed_residue_cutoff_used",
        "native_metrics_used_for_selection",
        "positive_folding_evidence_found",
        "folding_problem_solved",
        "provenance_clean",
    ]:
        print(f"{key}: {cert.get(key)}")
    print(f"failed_checks: {cert.get('failed_checks')}")
    print(f"missing_operator_evidence: {cert.get('missing_operator_evidence')}")
    print("")
    print("KcsA readout:")
    for key in [
        "source_boundary",
        "pore_filter_source_file_count",
        "assembly_interface_source_file_count",
        "pore_filter_row_count",
        "assembly_interface_row_count",
        "pore_filter_chains",
        "pore_filter_residue_numbers_or_indices",
        "pore_filter_motif_labels",
        "pore_filter_ion_names",
        "pore_filter_constraint_classes",
        "assembly_interface_chain_pair_count",
        "assembly_interface_chain_pairs_sample",
        "assembly_interface_constraint_classes",
        "operator_buckets_assigned",
        "constraint_backed_operator_readout_found",
    ]:
        print(f"  {key}: {readout.get(key)}")
    print("")
    print("Next decision:")
    print(json.dumps(cert.get("next_decision", {}), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
