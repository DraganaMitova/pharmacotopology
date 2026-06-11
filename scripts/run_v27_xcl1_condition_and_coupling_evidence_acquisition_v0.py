#!/usr/bin/env python3
from __future__ import annotations

"""V27 XCL1 condition/coupling evidence acquisition.

Fast zero-MD acquisition sprint after V26 selected XCL1 as the next evidence
acquisition target. This is not a fold-switch claim and not a simulation. It
locks the state-condition label context that is already safe, scans for local
state-specific coupling/constraint files without inventing them, preserves the
mixed-state leakage guard, and decides the next V28 mechanism-evidence step.
"""

import argparse
import glob
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V26_CERT = RUN_ROOT / "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST" / "v26_xcl1_state_separation_operator_test_certificate.json"
DEFAULT_V20_CERT = RUN_ROOT / "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT" / "v20_xcl1_state_specific_evidence_readout_certificate.json"
DEFAULT_PREFLIGHT_CERT = RUN_ROOT / "V16_TRANSFER_DATA_PREFLIGHT" / "v16_transfer_data_preflight_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION"

XCL1_COUPLING_PATTERNS = [
    "external_msa/**/*xcl1*",
    "external_msa/**/*XCL1*",
    "external_msa/**/*lymphotactin*",
    "external_msa/**/*Lymphotactin*",
    "external_msa_free_predictors/**/*xcl1*",
    "external_msa_free_predictors/**/*XCL1*",
    "data/**/*xcl1*coupl*",
    "data/**/*XCL1*coupl*",
    "data/**/*lymphotactin*coupl*",
    "data/**/*state*specific*coupl*",
    "data/**/*xcl1*constraint*",
    "data/**/*lymphotactin*constraint*",
]

XCL1_CONDITION_LABEL_PATTERNS = [
    "data/**/*xcl1*condition*.json",
    "data/**/*XCL1*condition*.json",
    "data/**/*lymphotactin*condition*.json",
    "data/**/*xcl1*state*label*.json",
    "data/**/*lymphotactin*state*label*.json",
    "data/**/*xcl1*monomer*dimer*.json",
    "data/**/*lymphotactin*monomer*dimer*.json",
]


def _read_json(path: Path, label: str, required: bool = True) -> dict[str, Any]:
    if not path.exists():
        if required:
            raise SystemExit(f"missing {label}: {path}")
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _bool(value: Any) -> bool:
    return value is True


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _glob_repo(patterns: list[str], limit: int = 120) -> list[str]:
    hits: list[str] = []
    for pattern in patterns:
        for hit in glob.glob(str(REPO_ROOT / pattern), recursive=True):
            p = Path(hit)
            if p.is_file():
                try:
                    hits.append(str(p.relative_to(REPO_ROOT)))
                except ValueError:
                    hits.append(str(p))
    return sorted(set(hits))[:limit]


def _target_row_by_id(cert: dict[str, Any], target_id: str) -> dict[str, Any]:
    for key in ("target_results", "target_rows", "targets"):
        rows = cert.get(key) if isinstance(cert.get(key), list) else []
        for row in rows:
            if isinstance(row, dict) and row.get("target_id") == target_id:
                return row
    return {}


def _state_a_b_from_v20(v20: dict[str, Any]) -> tuple[bool, bool]:
    if _bool(v20.get("state_A_role_evidence_found")) or _bool(v20.get("state_A_detected")):
        a = True
    else:
        a = False
    if _bool(v20.get("state_B_role_evidence_found")) or _bool(v20.get("state_B_detected")):
        b = True
    else:
        b = False
    return a, b


def _state_a_b_from_v26(v26: dict[str, Any]) -> tuple[bool, bool]:
    op = v26.get("operator_readout") if isinstance(v26.get("operator_readout"), dict) else {}
    return _bool(op.get("state_A_detected")), _bool(op.get("state_B_detected"))


def build_v27(v26: dict[str, Any], v20: dict[str, Any], preflight: dict[str, Any]) -> dict[str, Any]:
    target_id = "XCL1_lymphotactin"
    preflight_row = _target_row_by_id(preflight, target_id)

    v26_passed = v26.get("test_status") == "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST_PASSED_CLAIM_DISABLED" and _bool(v26.get("state_separation_operator_passed"))
    op = v26.get("operator_readout") if isinstance(v26.get("operator_readout"), dict) else {}
    v20_passed = v20.get("test_status") == "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED"

    a26, b26 = _state_a_b_from_v26(v26)
    a20, b20 = _state_a_b_from_v20(v20)
    state_a_locked = bool(a26 or a20)
    state_b_locked = bool(b26 or b20)
    state_labels_locked = bool(state_a_locked and state_b_locked)

    monomer_dimer_context = bool(
        v20.get("monomer_dimer_context_present") is True
        or "monomer_dimer_context" in _as_list(v20.get("available_evidence"))
        or "monomer_dimer_context" in _as_list(op.get("available_evidence"))
    )
    # The preflight has both state structures; for this stage, that is enough
    # to lock the context labels even if raw condition files are not local yet.
    preflight_ready = preflight.get("data_preflight_status") == "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT"
    state_structures_present = preflight_row.get("target_material_status") == "ready_for_zero_md_role_readout" or preflight_ready

    mixed_state_leakage_guard = bool(
        _bool(op.get("state_separation_guard_passed"))
        and not _bool(op.get("mixed_state_pollution"))
        and not _bool(op.get("mixed_state_contact_pooling_used"))
        and not _bool(op.get("single_fold_forcing"))
        and not _bool(op.get("fold_switch_claim_made"))
    )
    if not op and v20_passed:
        mixed_state_leakage_guard = bool(
            not _bool(v20.get("mixed_state_pollution"))
            and not _bool(v20.get("single_fold_claim_made"))
            and not _bool(v20.get("fold_switch_claim_made"))
        )

    condition_files = _glob_repo(XCL1_CONDITION_LABEL_PATTERNS)
    coupling_files = _glob_repo(XCL1_COUPLING_PATTERNS)
    condition_file_present = bool(condition_files)
    coupling_present = bool(coupling_files)

    condition_label_context_locked = bool(
        v26_passed
        and state_labels_locked
        and monomer_dimer_context
        and state_structures_present
        and mixed_state_leakage_guard
    )

    failed_checks: list[str] = []
    if not v26_passed:
        failed_checks.append("V26_state_separation_operator_passed")
    if not state_labels_locked:
        failed_checks.append("state_A_and_state_B_context_labels_locked")
    if not monomer_dimer_context:
        failed_checks.append("monomer_dimer_context_locked")
    if not state_structures_present:
        failed_checks.append("state_context_structures_present")
    if not mixed_state_leakage_guard:
        failed_checks.append("mixed_state_leakage_guard_preserved")

    positive_pressure = bool(condition_label_context_locked and not failed_checks)
    status = (
        "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION_PASSED_CLAIM_DISABLED"
        if positive_pressure
        else "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION_CLEAN_ABSTAIN_OR_BLOCKED_CLAIM_DISABLED"
    )

    manifest = {
        "kind": "V27_XCL1_CONDITION_LABEL_MANIFEST_v0",
        "target_id": target_id,
        "role_class": "metamorphic_switch_object",
        "mechanism_operator": "state_separation",
        "state_A_label": "state_A_chemokine_like_or_monomer_context",
        "state_B_label": "state_B_alternative_beta_sandwich_or_dimer_context",
        "state_label_source_policy": "local_state_context_structures_plus_prior_V20_V26_state_buckets_not_fold_switch_claim",
        "state_A_context_label_locked": state_a_locked,
        "state_B_context_label_locked": state_b_locked,
        "state_labels_locked": state_labels_locked,
        "monomer_dimer_context_locked": monomer_dimer_context,
        "condition_label_context_locked": condition_label_context_locked,
        "raw_local_condition_annotation_files_present": condition_file_present,
        "raw_local_condition_annotation_files": condition_files,
        "claim_allowed": False,
        "positive_folding_evidence_found": False,
    }

    preflight_obj = {
        "kind": "V27_XCL1_CONDITION_PREFLIGHT_v0",
        "target_id": target_id,
        "source_v26_status": v26.get("test_status"),
        "source_v20_status": v20.get("test_status"),
        "source_v16_preflight_status": preflight.get("data_preflight_status"),
        "condition_preflight_status": "ready_for_zero_md_condition_contrast" if condition_label_context_locked else "blocked_or_clean_abstain_condition_context_incomplete",
        "state_context_structures_present": state_structures_present,
        "state_A_context_label_locked": state_a_locked,
        "state_B_context_label_locked": state_b_locked,
        "monomer_dimer_context_locked": monomer_dimer_context,
        "mixed_state_leakage_guard_preserved": mixed_state_leakage_guard,
        "mixed_state_contact_pooling_used": False,
        "single_native_state_forced": False,
        "fold_switch_claim_made": False,
        "failed_checks": failed_checks,
        "claim_allowed": False,
        "new_md_executed": False,
    }

    coupling_scan = {
        "kind": "V27_XCL1_STATE_SPECIFIC_COUPLING_AVAILABILITY_v0",
        "target_id": target_id,
        "state_specific_couplings_present": coupling_present,
        "state_specific_external_constraints_present": coupling_present,
        "external_coupling_or_constraint_files": coupling_files,
        "external_couplings_required_for_V27": False,
        "external_couplings_required_for_future_MD_or_contact_test": True,
        "coupling_policy": "do_not_synthesize_state_specific_couplings_missing_is_reported_not_backfilled",
        "claim_allowed": False,
        "new_md_executed": False,
    }

    next_decision = {
        "kind": "V27_XCL1_NEXT_DECISION_v0",
        "decision_status": "V28_TARGET_SELECTED" if condition_label_context_locked else "V27_BLOCKED_OR_CLEAN_ABSTAIN",
        "selected_next_panel": "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST" if condition_label_context_locked else None,
        "selected_V28_target": target_id if condition_label_context_locked else None,
        "selected_V28_test": "XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST" if condition_label_context_locked else None,
        "reason": (
            "XCL1 state-condition labels and mixed-state leakage guard are locked; state-specific couplings remain optional for V28 and required before any MD/contact claim"
            if condition_label_context_locked
            else "XCL1 condition labels or mixed-state leakage guard are incomplete"
        ),
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "new_MD_allowed_policy": "only_after_state_specific_external_constraints_are_locked_and_false_win_risk_is_low",
        "parallel_MD_paths_allowed": False,
        "claim_allowed": False,
        "positive_folding_evidence_targets": [],
        "required_before_any_MD": [
            "state_specific_external_couplings_or_constraints_if_available",
            "condition_labels_or_state_context_annotations_preserved",
            "mixed_state_leakage_guard_preserved",
        ],
    }

    available = [
        "state_A_and_state_B_context_structures",
        "state_A_and_state_B_labels" if state_labels_locked else None,
        "monomer_dimer_context" if monomer_dimer_context else None,
        "mixed_state_leakage_guard" if mixed_state_leakage_guard else None,
        "state_specific_external_couplings_or_constraints" if coupling_present else None,
        "raw_local_condition_annotation_files" if condition_file_present else None,
    ]
    missing = []
    if not condition_file_present:
        missing.append("raw_local_condition_annotation_files_if_needed")
    if not coupling_present:
        missing.append("state_specific_external_couplings_or_constraints_if_available")
    if not coupling_present:
        missing.append("state_specific_coupling_support_for_future_MD_or_contact_test")

    return {
        "kind": "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION_v0",
        "run_mode": "zero_md_condition_coupling_acquisition_no_simulation_no_threshold_tuning",
        "test_status": status,
        "target_id": target_id,
        "role_class": "metamorphic_switch_object",
        "mechanism_operator": "state_separation",
        "source_v26_status": v26.get("test_status"),
        "condition_label_context_locked": condition_label_context_locked,
        "state_A_context_label_locked": state_a_locked,
        "state_B_context_label_locked": state_b_locked,
        "state_labels_locked": state_labels_locked,
        "monomer_dimer_context_locked": monomer_dimer_context,
        "mixed_state_leakage_guard_preserved": mixed_state_leakage_guard,
        "state_specific_couplings_present": coupling_present,
        "state_specific_external_constraints_present": coupling_present,
        "external_couplings_required_for_V27": False,
        "external_couplings_required_for_future_MD_or_contact_test": True,
        "positive_pressure_evidence_found": positive_pressure,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "mixed_state_contact_pooling_used": False,
        "mixed_state_pollution": False,
        "single_fold_forcing": False,
        "single_fold_claim_made": False,
        "fold_switch_claim_made": False,
        "forbidden_misclassification_violations": [],
        "failed_checks": failed_checks,
        "available_evidence": [x for x in available if x],
        "missing_evidence": missing,
        "condition_label_manifest": manifest,
        "condition_preflight": preflight_obj,
        "state_specific_coupling_availability": coupling_scan,
        "next_decision": next_decision,
        "locked_interpretation": (
            "V27 acquires the XCL1 condition/state-label layer and scans for state-specific couplings without inventing them. "
            "State A and state B remain separate; mixed-state fake core, single-fold forcing, fold-switch claim, MD, threshold tuning, and claim upgrade remain forbidden."
        ),
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V27 XCL1 Condition and Coupling Evidence Acquisition",
        "",
        f"Status: `{cert.get('test_status')}`",
        f"Condition label context locked: `{cert.get('condition_label_context_locked')}`",
        f"State-specific couplings present: `{cert.get('state_specific_couplings_present')}`",
        f"Positive pressure evidence: `{cert.get('positive_pressure_evidence_found')}`",
        f"Positive folding evidence: `{cert.get('positive_folding_evidence_found')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        "",
        "## Locked interpretation",
        cert.get("locked_interpretation", ""),
        "",
        "## Next decision",
        json.dumps(cert.get("next_decision", {}), indent=2, sort_keys=True),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v27_xcl1_condition_and_coupling_evidence_acquisition_certificate.json",
        "condition_label_manifest": out_dir / "v27_xcl1_condition_label_manifest.json",
        "condition_preflight": out_dir / "v27_xcl1_condition_preflight.json",
        "state_specific_coupling_availability": out_dir / "v27_xcl1_state_specific_coupling_availability.json",
        "next_decision": out_dir / "v27_xcl1_next_decision.json",
        "report": out_dir / "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {k: str(v) for k, v in paths.items()}
    paths["condition_label_manifest"].write_text(json.dumps(cert["condition_label_manifest"], indent=2, sort_keys=True), encoding="utf-8")
    paths["condition_preflight"].write_text(json.dumps(cert["condition_preflight"], indent=2, sort_keys=True), encoding="utf-8")
    paths["state_specific_coupling_availability"].write_text(json.dumps(cert["state_specific_coupling_availability"], indent=2, sort_keys=True), encoding="utf-8")
    paths["next_decision"].write_text(json.dumps(cert["next_decision"], indent=2, sort_keys=True), encoding="utf-8")
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v26-cert", type=Path, default=DEFAULT_V26_CERT)
    parser.add_argument("--v20-cert", type=Path, default=DEFAULT_V20_CERT)
    parser.add_argument("--preflight-cert", type=Path, default=DEFAULT_PREFLIGHT_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    v26 = _read_json(args.v26_cert, "V26 XCL1 state-separation operator certificate")
    v20 = _read_json(args.v20_cert, "V20 XCL1 state-specific evidence certificate")
    preflight = _read_json(args.preflight_cert, "V16 transfer data preflight certificate")
    cert = build_v27(v26, v20, preflight)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "test_status": cert["test_status"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "positive_pressure_evidence_found": cert["positive_pressure_evidence_found"],
        "positive_folding_evidence_found": cert["positive_folding_evidence_found"],
        "condition_label_context_locked": cert["condition_label_context_locked"],
        "state_specific_couplings_present": cert["state_specific_couplings_present"],
        "selected_next_panel": cert["next_decision"].get("selected_next_panel"),
        "certificate": str(paths["certificate"]),
        "condition_label_manifest": str(paths["condition_label_manifest"]),
        "condition_preflight": str(paths["condition_preflight"]),
        "state_specific_coupling_availability": str(paths["state_specific_coupling_availability"]),
        "next_decision": str(paths["next_decision"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
