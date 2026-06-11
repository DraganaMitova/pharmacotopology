#!/usr/bin/env python3
from __future__ import annotations

"""V26 XCL1 state-separation operator test.

This is the first locked mechanism-operator test after V25 selected XCL1 as
V26 target. It is not a new MD run, not a fold-switch claim, and not a
single-native-state selection. The test asks whether the repeated mechanism
operator `state_separation` remains valid for XCL1: state A and state B must
stay in separate evidence buckets, mixed-state fake core selection must remain
forbidden, and claim_allowed must remain false.
"""

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V25_CERT = RUN_ROOT / "V25_FAST_MECHANISM_EVIDENCE_SPRINT" / "v25_fast_mechanism_evidence_sprint_certificate.json"
DEFAULT_V20_CERT = RUN_ROOT / "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT" / "v20_xcl1_state_specific_evidence_readout_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST"


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _bool(value: Any) -> bool:
    return bool(value is True)


def _list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _xcl1_from_v25(v25: dict[str, Any]) -> dict[str, Any]:
    x = v25.get("xcl1_state_specific_readout")
    return x if isinstance(x, dict) else {}


def build_v26(v25: dict[str, Any], v20: dict[str, Any]) -> dict[str, Any]:
    decision = v25.get("next_mechanism_test_decision") if isinstance(v25.get("next_mechanism_test_decision"), dict) else {}
    x25 = _xcl1_from_v25(v25)

    selected_by_v25 = (
        decision.get("selected_V26_target") == "XCL1_lymphotactin"
        and decision.get("selected_V26_test") == "XCL1_STATE_SEPARATION_OPERATOR_TEST"
    )
    state_a = _bool(x25.get("state_A_detected")) or _bool(v20.get("state_A_role_evidence_found"))
    state_b = _bool(x25.get("state_B_detected")) or _bool(v20.get("state_B_role_evidence_found"))
    state_specific = _bool(x25.get("state_specific_role_evidence_found")) or _bool(v20.get("state_specific_role_evidence_found"))
    mixed_pollution = _bool(x25.get("mixed_state_pollution")) or _bool(v20.get("mixed_state_pollution"))
    single_forcing = _bool(x25.get("single_fold_forcing")) or _bool(v20.get("single_fold_forcing")) or _bool(v20.get("single_fold_claim_made"))
    mixed_pooling = _bool(x25.get("mixed_state_contact_pooling_used")) or _bool(v20.get("mixed_state_contact_pooling_used"))
    fold_switch_claim = _bool(x25.get("fold_switch_claim_made")) or _bool(v20.get("fold_switch_claim_made"))

    no_claim = (
        v25.get("claim_allowed") is False
        and v25.get("new_md_executed") is False
        and v20.get("claim_allowed") is False
        and v20.get("new_md_executed") is False
        and v20.get("positive_folding_evidence_found") is False
    )
    state_separation_operator_passed = bool(
        selected_by_v25
        and state_a
        and state_b
        and state_specific
        and not mixed_pollution
        and not single_forcing
        and not mixed_pooling
        and not fold_switch_claim
        and no_claim
    )

    failed_checks: list[str] = []
    if not selected_by_v25:
        failed_checks.append("V25_selected_XCL1_state_separation_operator")
    if not state_a:
        failed_checks.append("state_A_detected")
    if not state_b:
        failed_checks.append("state_B_detected")
    if not state_specific:
        failed_checks.append("state_specific_role_evidence_found")
    if mixed_pollution:
        failed_checks.append("mixed_state_pollution_absent")
    if single_forcing:
        failed_checks.append("single_fold_forcing_absent")
    if mixed_pooling:
        failed_checks.append("mixed_state_contact_pooling_absent")
    if fold_switch_claim:
        failed_checks.append("fold_switch_claim_absent")
    if not no_claim:
        failed_checks.append("claim_and_MD_boundaries_clean")

    status = (
        "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST_PASSED_CLAIM_DISABLED"
        if state_separation_operator_passed
        else "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST_CLEAN_ABSTAIN_OR_BLOCKED_CLAIM_DISABLED"
    )
    missing = sorted(set(str(x) for x in (_list(x25.get("missing_evidence")) + _list(v20.get("missing_evidence")))))
    if "condition_labels_if_available" not in missing:
        missing.append("condition_labels_if_available")
    if "state_specific_external_couplings_or_constraints_if_available" not in missing:
        missing.append("state_specific_external_couplings_or_constraints_if_available")

    operator_readout = {
        "kind": "V26_XCL1_STATE_SEPARATION_OPERATOR_READOUT_v0",
        "target_id": "XCL1_lymphotactin",
        "role_class": "metamorphic_switch_object",
        "mechanism_operator": "state_separation",
        "operator_status": "passed" if state_separation_operator_passed else "clean_abstain_or_blocked",
        "state_A_detected": state_a,
        "state_B_detected": state_b,
        "state_specific_role_evidence_found": state_specific,
        "state_specific_buckets_preserved": bool(state_a and state_b and state_specific),
        "mixed_state_fake_core_selected": False,
        "mixed_state_pollution": mixed_pollution,
        "mixed_state_contact_pooling_used": mixed_pooling,
        "single_fold_forcing": single_forcing,
        "single_fold_claim_made": single_forcing,
        "fold_switch_claim_made": fold_switch_claim,
        "state_separation_guard_passed": bool(not mixed_pollution and not single_forcing and not mixed_pooling),
        "forbidden_misclassification_violations": [],
        "selection_threshold_used": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "positive_pressure_evidence_found": state_separation_operator_passed,
        "positive_folding_evidence_found": False,
        "available_evidence": sorted(set(str(x) for x in (_list(x25.get("available_evidence")) + _list(v20.get("available_evidence"))))),
        "missing_evidence": missing,
        "operator_policy": "state_specific_buckets_no_cross_state_pooling_no_single_native_assumption_no_fold_switch_claim",
    }
    next_decision = {
        "kind": "V26_NEXT_MECHANISM_DECISION_v0",
        "decision_status": "V26_OPERATOR_TEST_LOCKED_NEXT_EVIDENCE_ACQUISITION" if state_separation_operator_passed else "V26_OPERATOR_TEST_BLOCKED",
        "selected_next_panel": "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION" if state_separation_operator_passed else None,
        "reason": (
            "XCL1 state-separation operator is clean and false-win risk remains controlled; next step is condition/coupling acquisition, not MD"
            if state_separation_operator_passed
            else "XCL1 state-separation operator is not clean enough for the next acquisition step"
        ),
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "new_MD_allowed_policy": "only_after_condition_labels_and_state_specific_external_constraints_are_locked_and_false_win_risk_is_low",
        "parallel_MD_paths_allowed": False,
        "claim_allowed": False,
        "positive_folding_evidence_targets": [],
        "required_before_any_MD": [
            "condition_labels_or_state_context_annotations",
            "state_specific_external_couplings_or_constraints_if_available",
            "explicit_mixed_state_leakage_guard_preserved",
        ],
    }
    return {
        "kind": "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST_v0",
        "run_mode": "zero_md_state_separation_operator_test_no_simulation_no_threshold_tuning",
        "test_status": status,
        "target_id": "XCL1_lymphotactin",
        "role_class": "metamorphic_switch_object",
        "mechanism_operator_tested": "state_separation",
        "source_v25_status": v25.get("sprint_status"),
        "source_v20_status": v20.get("test_status"),
        "selected_by_v25": selected_by_v25,
        "state_separation_operator_passed": state_separation_operator_passed,
        "mechanism_operator_evidence_found": state_separation_operator_passed,
        "positive_pressure_evidence_found": state_separation_operator_passed,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "failed_checks": failed_checks,
        "operator_readout": operator_readout,
        "next_mechanism_decision": next_decision,
        "locked_interpretation": (
            "V26 tests the repeated mechanism operator `state_separation` on XCL1. State A and state B remain separate; "
            "mixed-state fake core, single-fold forcing, fold-switch claim, native-metric selection, and MD remain forbidden. "
            "This is mechanism-operator evidence, not positive folding evidence."
        ),
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V26 XCL1 State-Separation Operator Test",
        "",
        f"Status: `{cert.get('test_status')}`",
        f"Mechanism operator: `{cert.get('mechanism_operator_tested')}`",
        f"State-separation operator passed: `{cert.get('state_separation_operator_passed')}`",
        f"Positive pressure evidence: `{cert.get('positive_pressure_evidence_found')}`",
        f"Positive folding evidence: `{cert.get('positive_folding_evidence_found')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        "",
        "## Locked interpretation",
        cert.get("locked_interpretation", ""),
        "",
        "## Next decision",
        json.dumps(cert.get("next_mechanism_decision", {}), indent=2, sort_keys=True),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v26_xcl1_state_separation_operator_test_certificate.json",
        "operator_readout": out_dir / "v26_xcl1_state_separation_operator_readout.json",
        "decision": out_dir / "v26_next_mechanism_decision.json",
        "report": out_dir / "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["operator_readout"].write_text(json.dumps(cert["operator_readout"], indent=2, sort_keys=True), encoding="utf-8")
    paths["decision"].write_text(json.dumps(cert["next_mechanism_decision"], indent=2, sort_keys=True), encoding="utf-8")
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v25-cert", type=Path, default=DEFAULT_V25_CERT)
    parser.add_argument("--v20-cert", type=Path, default=DEFAULT_V20_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    v25 = _read_json(args.v25_cert, "V25 fast mechanism sprint certificate")
    v20 = _read_json(args.v20_cert, "V20 XCL1 state-specific evidence certificate")
    cert = build_v26(v25, v20)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "test_status": cert["test_status"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "positive_pressure_evidence_found": cert["positive_pressure_evidence_found"],
        "positive_folding_evidence_found": cert["positive_folding_evidence_found"],
        "state_separation_operator_passed": cert["state_separation_operator_passed"],
        "selected_next_panel": cert["next_mechanism_decision"].get("selected_next_panel"),
        "certificate": str(paths["certificate"]),
        "operator_readout": str(paths["operator_readout"]),
        "decision": str(paths["decision"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
