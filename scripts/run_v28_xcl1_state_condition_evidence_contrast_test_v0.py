#!/usr/bin/env python3
from __future__ import annotations

"""V28 XCL1 state-condition evidence contrast test.

Zero-MD mechanism/evidence contrast after V26/V27. The test verifies that
XCL1's two valid state contexts remain separated under the locked condition
labels:

  state A -> chemokine-like / monomer context
  state B -> alternative beta-sandwich / dimer context

It is not a fold-switch claim, does not pool contacts across states, does not
force a single native state, and does not require state-specific couplings for
this contrast-only stage. Couplings remain required before any MD/contact claim.
"""

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V27_CERT = RUN_ROOT / "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION" / "v27_xcl1_condition_and_coupling_evidence_acquisition_certificate.json"
DEFAULT_V26_CERT = RUN_ROOT / "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST" / "v26_xcl1_state_separation_operator_test_certificate.json"
DEFAULT_V20_CERT = RUN_ROOT / "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT" / "v20_xcl1_state_specific_evidence_readout_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST"


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _bool(value: Any) -> bool:
    return value is True


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _state_readout(v20: dict[str, Any], state_key: str) -> dict[str, Any]:
    readouts = v20.get("state_specific_readouts") if isinstance(v20.get("state_specific_readouts"), dict) else {}
    item = readouts.get(state_key) if isinstance(readouts.get(state_key), dict) else {}
    return dict(item)


def build_v28(v27: dict[str, Any], v26: dict[str, Any], v20: dict[str, Any]) -> dict[str, Any]:
    target_id = "XCL1_lymphotactin"
    failed: list[str] = []
    violations: list[str] = []

    v27_passed = (
        v27.get("test_status") == "V27_XCL1_CONDITION_AND_COUPLING_EVIDENCE_ACQUISITION_PASSED_CLAIM_DISABLED"
        and _bool(v27.get("condition_label_context_locked"))
        and _bool(v27.get("state_A_context_label_locked"))
        and _bool(v27.get("state_B_context_label_locked"))
        and _bool(v27.get("monomer_dimer_context_locked"))
        and _bool(v27.get("mixed_state_leakage_guard_preserved"))
    )
    if not v27_passed:
        failed.append("V27_condition_state_label_context_locked")

    next_decision = v27.get("next_decision") if isinstance(v27.get("next_decision"), dict) else {}
    if next_decision.get("selected_next_panel") != "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST":
        failed.append("V27_selected_V28_state_condition_contrast")

    v26_passed = (
        v26.get("test_status") == "V26_XCL1_STATE_SEPARATION_OPERATOR_TEST_PASSED_CLAIM_DISABLED"
        and _bool(v26.get("state_separation_operator_passed"))
    )
    if not v26_passed:
        failed.append("V26_state_separation_operator_passed")

    v20_passed = (
        v20.get("test_status") == "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED"
        and _bool(v20.get("state_A_role_evidence_found"))
        and _bool(v20.get("state_B_role_evidence_found"))
        and _bool(v20.get("state_specific_role_evidence_found"))
    )
    if not v20_passed:
        failed.append("V20_state_specific_role_evidence_found")

    state_a = _state_readout(v20, "state_A")
    state_b = _state_readout(v20, "state_B")
    state_a_condition_evidence = bool(
        _bool(v27.get("state_A_context_label_locked"))
        and _bool(v20.get("state_A_role_evidence_found"))
        and state_a.get("state_role_bucket") in {"state_A_chemokine_monomer_support_context", "state_A_chemokine_monomer_support"}
    )
    state_b_condition_evidence = bool(
        _bool(v27.get("state_B_context_label_locked"))
        and _bool(v20.get("state_B_role_evidence_found"))
        and state_b.get("state_role_bucket") in {"state_B_beta_sandwich_dimer_or_alternative_state_support_context", "state_B_beta_sandwich_dimer_support"}
    )
    if not state_a_condition_evidence:
        failed.append("state_A_condition_evidence_found")
    if not state_b_condition_evidence:
        failed.append("state_B_condition_evidence_found")

    # These are deliberately false by policy; V28 must never mix states or force
    # a single native fold to manufacture a win.
    mixed_state_contact_pooling_used = False
    mixed_state_fake_core_selected = False
    single_native_state_forced = False
    single_fold_claim_made = False
    fold_switch_claim_made = False
    if mixed_state_contact_pooling_used or mixed_state_fake_core_selected:
        violations.append("mixing_two_states_into_one_false_core")
    if single_native_state_forced or single_fold_claim_made:
        violations.append("forcing_single_canonical_fold")
    if fold_switch_claim_made:
        violations.append("claiming_fold_switch_without_state_specific_support")

    # Preserve upstream violations if any appeared.
    for upstream in [v27, v26, v20]:
        for violation in _as_list(upstream.get("forbidden_misclassification_violations")):
            violations.append(str(violation))
    violations = sorted(set(violations))
    if violations:
        failed.append("forbidden_misclassification_violations_empty")

    for label, cert in [("V27", v27), ("V26", v26), ("V20", v20)]:
        if cert.get("new_md_executed") is not False:
            failed.append(f"{label}_no_new_md")
        if cert.get("fixed_residue_cutoff_used") is not False:
            failed.append(f"{label}_fixed_residue_cutoff_retired")
        if cert.get("native_metrics_used_for_selection") is not False:
            failed.append(f"{label}_native_metric_selection_forbidden")
        if cert.get("claim_allowed") is not False:
            failed.append(f"{label}_claim_disabled")
        if cert.get("positive_folding_evidence_found") is not False:
            failed.append(f"{label}_positive_folding_evidence_not_allowed")

    state_specific_buckets_preserved = bool(state_a_condition_evidence and state_b_condition_evidence and not violations)
    state_condition_contrast_preserved = bool(v27_passed and v26_passed and v20_passed and state_specific_buckets_preserved)
    pass_status = bool(state_condition_contrast_preserved and not failed)

    status = (
        "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST_PASSED_CLAIM_DISABLED"
        if pass_status
        else "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST_BLOCKED_OR_CLEAN_ABSTAIN_CLAIM_DISABLED"
    )

    coupling = v27.get("state_specific_coupling_availability") if isinstance(v27.get("state_specific_coupling_availability"), dict) else {}
    coupling_files = _as_list(coupling.get("external_coupling_or_constraint_files"))
    state_specific_couplings_present = _bool(v27.get("state_specific_couplings_present")) or _bool(coupling.get("state_specific_couplings_present")) or bool(coupling_files)

    available_evidence = [
        item for item, present in {
            "state_A_condition_label_context": _bool(v27.get("state_A_context_label_locked")),
            "state_B_condition_label_context": _bool(v27.get("state_B_context_label_locked")),
            "monomer_dimer_context": _bool(v27.get("monomer_dimer_context_locked")),
            "state_A_state_B_context_structures": True,
            "state_A_condition_role_evidence": state_a_condition_evidence,
            "state_B_condition_role_evidence": state_b_condition_evidence,
            "mixed_state_leakage_guard": not violations,
            "state_specific_external_couplings_or_constraints": state_specific_couplings_present,
        }.items() if present
    ]
    missing_evidence = []
    if not state_specific_couplings_present:
        missing_evidence.append("state_specific_external_couplings_or_constraints_if_available")
    missing_evidence.append("independent_condition_labels_or_experimental_context_if_available")
    missing_evidence.append("state_specific_coupling_support_for_future_MD_or_contact_test")

    state_condition_contrast = {
        "kind": "V28_XCL1_STATE_CONDITION_CONTRAST_READOUT_v0",
        "target_id": target_id,
        "state_A": {
            "label": v27.get("condition_label_manifest", {}).get("state_A_label", "state_A_chemokine_like_or_monomer_context") if isinstance(v27.get("condition_label_manifest"), dict) else "state_A_chemokine_like_or_monomer_context",
            "source_pdb": state_a.get("pdb_id", "2HDM"),
            "condition_role_evidence_found": state_a_condition_evidence,
            "role_bucket": state_a.get("state_role_bucket"),
            "expected_context": state_a.get("expected_context"),
            "chain_ca_counts": state_a.get("chain_ca_counts", {}),
            "model_count": state_a.get("model_count"),
        },
        "state_B": {
            "label": v27.get("condition_label_manifest", {}).get("state_B_label", "state_B_alternative_beta_sandwich_or_dimer_context") if isinstance(v27.get("condition_label_manifest"), dict) else "state_B_alternative_beta_sandwich_or_dimer_context",
            "source_pdb": state_b.get("pdb_id", "2JP1"),
            "condition_role_evidence_found": state_b_condition_evidence,
            "role_bucket": state_b.get("state_role_bucket"),
            "expected_context": state_b.get("expected_context"),
            "chain_ca_counts": state_b.get("chain_ca_counts", {}),
            "model_count": state_b.get("model_count"),
            "interface_pairs_present": (state_b.get("chain_interface_readout") or {}).get("interface_pairs_present", []),
        },
        "guard": {
            "guard_policy": "state_condition_buckets_no_cross_state_pooling_no_single_native_assumption_no_fold_switch_claim",
            "state_A_state_B_contact_pooling_used": mixed_state_contact_pooling_used,
            "mixed_state_fake_core_selected": mixed_state_fake_core_selected,
            "single_native_state_forced": single_native_state_forced,
            "single_fold_claim_made": single_fold_claim_made,
            "fold_switch_claim_made": fold_switch_claim_made,
            "selection_threshold_used": False,
            "state_specific_buckets_preserved": state_specific_buckets_preserved,
        },
    }

    next_decision_out = {
        "kind": "V28_XCL1_NEXT_DECISION_v0",
        "decision_status": "V28_CONTRAST_LOCKED_NEXT_EXTERNAL_CONSTRAINT_ACQUISITION" if pass_status else "V28_BLOCKED_OR_CLEAN_ABSTAIN",
        "selected_next_panel": "V29_MECHANISM_OPERATOR_PANEL_SUMMARY_AND_MD_READINESS_DECISION" if pass_status else None,
        "reason": (
            "XCL1 condition/state contrast is preserved and false-win risk remains controlled; summarize mechanism operators and MD readiness before any simulation"
            if pass_status else "XCL1 state-condition contrast incomplete or guard failed"
        ),
        "new_MD_allowed": False,
        "new_MD_recommended": False,
        "new_MD_allowed_policy": "only_after_state_specific_external_constraints_are_locked_and_false_win_risk_is_low",
        "parallel_MD_paths_allowed": False,
        "claim_allowed": False,
        "required_before_any_MD": [
            "state_specific_external_couplings_or_constraints_if_available",
            "independent_condition_or_oligomerization_context_if_available",
            "mixed_state_leakage_guard_preserved",
        ],
    }

    return {
        "kind": "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST_v0",
        "run_mode": "zero_md_state_condition_evidence_contrast_no_simulation_no_threshold_tuning",
        "test_status": status,
        "target_id": target_id,
        "role_class": "metamorphic_switch_object",
        "mechanism_operator_tested": "state_condition_separation",
        "source_v27_status": v27.get("test_status"),
        "source_v26_status": v26.get("test_status"),
        "source_v20_status": v20.get("test_status"),
        "condition_label_context_locked": _bool(v27.get("condition_label_context_locked")),
        "state_A_condition_evidence_found": state_a_condition_evidence,
        "state_B_condition_evidence_found": state_b_condition_evidence,
        "state_condition_contrast_preserved": state_condition_contrast_preserved,
        "state_specific_buckets_preserved": state_specific_buckets_preserved,
        "mixed_state_pollution": bool(mixed_state_contact_pooling_used or mixed_state_fake_core_selected),
        "mixed_state_contact_pooling_used": mixed_state_contact_pooling_used,
        "mixed_state_fake_core_selected": mixed_state_fake_core_selected,
        "single_fold_forcing": single_native_state_forced,
        "single_fold_claim_made": single_fold_claim_made,
        "fold_switch_claim_made": fold_switch_claim_made,
        "state_specific_couplings_present": state_specific_couplings_present,
        "external_couplings_required_for_V28": False,
        "external_couplings_required_for_future_MD_or_contact_test": True,
        "positive_pressure_evidence_found": pass_status,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "forbidden_misclassification_violations": violations,
        "failed_checks": sorted(set(failed)),
        "available_evidence": sorted(set(available_evidence)),
        "missing_evidence": sorted(set(missing_evidence)),
        "state_condition_contrast_readout": state_condition_contrast,
        "next_decision": next_decision_out,
        "locked_interpretation": (
            "V28 locks the XCL1 state-condition contrast: state A and state B keep separate condition/evidence buckets, no mixed-state fake core is selected, no single canonical fold is forced, and no fold-switch claim is made. This is pressure-context evidence, not positive folding evidence."
        ),
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V28 XCL1 State-Condition Evidence Contrast Test",
        "",
        f"Status: `{cert.get('test_status')}`",
        f"Mechanism operator tested: `{cert.get('mechanism_operator_tested')}`",
        f"Positive pressure evidence: `{cert.get('positive_pressure_evidence_found')}`",
        f"Positive folding evidence: `{cert.get('positive_folding_evidence_found')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        "",
        "## State contrast",
        f"State A condition evidence found: `{cert.get('state_A_condition_evidence_found')}`",
        f"State B condition evidence found: `{cert.get('state_B_condition_evidence_found')}`",
        f"State condition contrast preserved: `{cert.get('state_condition_contrast_preserved')}`",
        f"State-specific buckets preserved: `{cert.get('state_specific_buckets_preserved')}`",
        "",
        "## Guards",
        f"Mixed-state pollution: `{cert.get('mixed_state_pollution')}`",
        f"Mixed-state contact pooling used: `{cert.get('mixed_state_contact_pooling_used')}`",
        f"Single-fold claim made: `{cert.get('single_fold_claim_made')}`",
        f"Fold-switch claim made: `{cert.get('fold_switch_claim_made')}`",
        f"Forbidden violations: `{cert.get('forbidden_misclassification_violations')}`",
        "",
        "## Missing evidence",
        f"`{cert.get('missing_evidence')}`",
        "",
        "## Locked interpretation",
        cert.get("locked_interpretation", ""),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v28_xcl1_state_condition_evidence_contrast_certificate.json",
        "state_condition_contrast_readout": out_dir / "v28_xcl1_state_condition_contrast_readout.json",
        "next_decision": out_dir / "v28_xcl1_next_decision.json",
        "report": out_dir / "V28_XCL1_STATE_CONDITION_EVIDENCE_CONTRAST_TEST_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {k: str(v) for k, v in paths.items()}
    paths["state_condition_contrast_readout"].write_text(json.dumps(cert["state_condition_contrast_readout"], indent=2, sort_keys=True), encoding="utf-8")
    paths["next_decision"].write_text(json.dumps(cert["next_decision"], indent=2, sort_keys=True), encoding="utf-8")
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v27-cert", type=Path, default=DEFAULT_V27_CERT)
    parser.add_argument("--v26-cert", type=Path, default=DEFAULT_V26_CERT)
    parser.add_argument("--v20-cert", type=Path, default=DEFAULT_V20_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    v27 = _read_json(args.v27_cert, "V27 XCL1 condition/coupling acquisition certificate")
    v26 = _read_json(args.v26_cert, "V26 XCL1 state-separation operator certificate")
    v20 = _read_json(args.v20_cert, "V20 XCL1 state-specific evidence certificate")
    cert = build_v28(v27, v26, v20)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "test_status": cert["test_status"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "positive_pressure_evidence_found": cert["positive_pressure_evidence_found"],
        "positive_folding_evidence_found": cert["positive_folding_evidence_found"],
        "state_condition_contrast_preserved": cert["state_condition_contrast_preserved"],
        "state_specific_buckets_preserved": cert["state_specific_buckets_preserved"],
        "selected_next_panel": cert["next_decision"].get("selected_next_panel"),
        "certificate": str(paths["certificate"]),
        "state_condition_contrast_readout": str(paths["state_condition_contrast_readout"]),
        "next_decision": str(paths["next_decision"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
