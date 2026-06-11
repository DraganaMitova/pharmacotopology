#!/usr/bin/env python3
from __future__ import annotations

"""V16 zero-MD role-transfer readout.

This stage reads the locked V16 manifest and the V16 data preflight certificate,
then performs a conservative role-context readout for the new pressure classes.
It does not run MD, does not tune thresholds, does not use native metrics for
selection, and does not permit folding claims.
"""

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_ROLE_MANIFEST = REPO_ROOT / "data" / "v16_locked_grammar_transfer_target_manifest.json"
DEFAULT_PREFLIGHT_CERT = RUN_ROOT / "V16_TRANSFER_DATA_PREFLIGHT" / "v16_transfer_data_preflight_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V16_ZERO_MD_ROLE_TRANSFER_READOUT"
REQUIRED_TARGET_IDS = {"p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # pragma: no cover
        raise SystemExit(f"invalid JSON in {label}: {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _role_manifest_by_id(role_manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    targets = role_manifest.get("targets") if isinstance(role_manifest.get("targets"), list) else []
    return {str(t.get("target_id")): t for t in targets if isinstance(t, dict)}


def _material_by_id(target_result: dict[str, Any]) -> dict[str, dict[str, Any]]:
    mats = target_result.get("material_results") if isinstance(target_result.get("material_results"), list) else []
    return {str(m.get("material_id")): m for m in mats if isinstance(m, dict)}


def _chain_counts(mat: dict[str, Any]) -> dict[str, int]:
    raw = mat.get("pdb_summary", {}).get("chain_ca_counts", {}) if isinstance(mat.get("pdb_summary"), dict) else {}
    if not isinstance(raw, dict):
        return {}
    out: dict[str, int] = {}
    for key, value in raw.items():
        try:
            out[str(key)] = int(value)
        except Exception:
            continue
    return out


def _base_target(target_result: dict[str, Any], role_target: dict[str, Any]) -> dict[str, Any]:
    return {
        "target_id": target_result.get("target_id"),
        "pressure_class": target_result.get("pressure_class"),
        "expected_role_class": target_result.get("expected_role_class"),
        "target_material_status": target_result.get("target_material_status"),
        "target_role_preflight_status": target_result.get("target_role_preflight_status"),
        "clean_abstain_allowed": target_result.get("clean_abstain_allowed") is True,
        "allowed_evidence_roles": role_target.get("allowed_evidence_roles", []),
        "forbidden_misclassification": role_target.get("forbidden_misclassification", []),
        "claim_allowed": False,
        "new_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "selected_core_or_clean_abstain": None,
        "role_buckets_assigned": [],
        "monitor_only_roles": [],
        "clean_abstain_roles": [],
        "forbidden_misclassification_violations": [],
        "positive_role_context_found": False,
        "target_transfer_status": "not_evaluated",
        "limitations": [],
        "material_observations": {},
    }


def _p53_readout(target_result: dict[str, Any], role_target: dict[str, Any]) -> dict[str, Any]:
    row = _base_target(target_result, role_target)
    if target_result.get("target_material_status") != "ready_for_zero_md_role_readout":
        row["target_transfer_status"] = "blocked_preflight_not_ready"
        row["clean_abstain_roles"] = ["clean_abstain_missing_partner_complex_context"]
        return row

    mats = _material_by_id(target_result)
    mat = mats.get("p53_TAD_MDM2_complex_structure", {})
    chains = _chain_counts(mat)
    values = sorted(chains.values())
    small_chain = values[0] if values else 0
    large_chain = values[-1] if values else 0
    partner_complex_present = len(chains) >= 2 and small_chain <= 30 and large_chain >= 50

    row["material_observations"] = {
        "chain_ca_counts": chains,
        "small_chain_ca_count": small_chain,
        "large_chain_ca_count": large_chain,
        "partner_complex_context_present": partner_complex_present,
    }
    if partner_complex_present:
        row["selected_core_or_clean_abstain"] = "partner_induced_interface_or_helix_signal_found"
        row["role_buckets_assigned"] = [
            "partner_induced_interface_or_helix_signal_found",
            "binding_stabilized_interface_core_context",
        ]
        row["clean_abstain_roles"] = ["clean_abstain_no_autonomous_isolated_TAD_core_claim"]
        row["positive_role_context_found"] = True
        row["target_transfer_status"] = "partner_induced_complex_role_context_detected_no_autonomous_fold_claim"
    else:
        row["selected_core_or_clean_abstain"] = "clean_abstain_no_partner_induced_context"
        row["clean_abstain_roles"] = ["clean_abstain_no_partner_induced_context"]
        row["target_transfer_status"] = "clean_abstain_partner_context_not_readable"
    row["limitations"] = [
        "isolated_TAD_material_not_required_at_this_stage",
        "disorder_prior_and_external_couplings_later",
        "no_universal_p53_fold_claim",
    ]
    return row


def _kcsa_readout(target_result: dict[str, Any], role_target: dict[str, Any]) -> dict[str, Any]:
    row = _base_target(target_result, role_target)
    if target_result.get("target_material_status") != "ready_for_zero_md_role_readout":
        row["target_transfer_status"] = "blocked_preflight_not_ready"
        row["clean_abstain_roles"] = ["clean_abstain_missing_membrane_or_oligomer_context"]
        return row

    mats = _material_by_id(target_result)
    mat = mats.get("KcsA_Fab_complex_structure", {})
    chains = _chain_counts(mat)
    ca_total = int(mat.get("pdb_summary", {}).get("ca_atom_count") or 0) if isinstance(mat.get("pdb_summary"), dict) else 0
    complex_context_present = ca_total >= 100 and len(chains) >= 3

    row["material_observations"] = {
        "chain_ca_counts": chains,
        "ca_atom_count": ca_total,
        "complex_chain_context_present": complex_context_present,
        "biological_tetramer_annotation_present": False,
        "membrane_topology_annotation_present": False,
    }
    if complex_context_present:
        row["selected_core_or_clean_abstain"] = "membrane_pore_roles_detected_without_soluble_core_misclassification"
        row["role_buckets_assigned"] = [
            "membrane_pore_oligomer_context",
            "transmembrane_helix_scaffold_context",
        ]
        row["monitor_only_roles"] = [
            "pore_selectivity_filter_core_later_annotation_required",
            "tetramer_interface_support_later_biological_assembly_required",
            "ion_filter_diagnostic_shell_later_annotation_required",
        ]
        row["positive_role_context_found"] = True
        row["target_transfer_status"] = "membrane_pore_context_detected_soluble_core_misclassification_avoided_no_tetramer_claim"
    else:
        row["selected_core_or_clean_abstain"] = "clean_abstain_missing_membrane_or_oligomer_context"
        row["clean_abstain_roles"] = ["clean_abstain_missing_membrane_or_oligomer_context"]
        row["target_transfer_status"] = "clean_abstain_kcsa_context_not_readable"
    row["limitations"] = [
        "biological_assembly_or_oligomer_annotation_later",
        "membrane_topology_annotation_later",
        "pore_filter_residue_context_later",
        "no_heavy_membrane_MD_executed",
        "no_whole_fold_claim",
    ]
    return row


def _xcl1_readout(target_result: dict[str, Any], role_target: dict[str, Any]) -> dict[str, Any]:
    row = _base_target(target_result, role_target)
    if target_result.get("target_material_status") != "ready_for_zero_md_role_readout":
        row["target_transfer_status"] = "blocked_preflight_not_ready"
        row["clean_abstain_roles"] = ["clean_state_specific_abstain_missing_state_context"]
        return row

    mats = _material_by_id(target_result)
    state_a = mats.get("XCL1_state_A_chemokine_like_structure", {})
    state_b = mats.get("XCL1_state_B_alternative_conformation_structure", {})
    a_chains = _chain_counts(state_a)
    b_chains = _chain_counts(state_b)
    two_state_context = bool(a_chains) and bool(b_chains)

    row["material_observations"] = {
        "state_A_chain_ca_counts": a_chains,
        "state_B_chain_ca_counts": b_chains,
        "two_state_context_present": two_state_context,
        "state_specific_external_couplings_present": False,
    }
    if two_state_context:
        row["selected_core_or_clean_abstain"] = "multiple_state_roles_supported_or_clean_state_specific_abstain"
        row["role_buckets_assigned"] = [
            "state_A_chemokine_monomer_support_context",
            "state_B_alternative_or_beta_sandwich_state_support_context",
            "pressure_condition_switch_context",
        ]
        row["monitor_only_roles"] = ["shared_local_support_later_coupling_readout_required"]
        row["positive_role_context_found"] = True
        row["target_transfer_status"] = "metamorphic_two_state_context_detected_single_fold_forcing_avoided_no_switch_claim"
    else:
        row["selected_core_or_clean_abstain"] = "clean_state_specific_abstain"
        row["clean_abstain_roles"] = ["clean_state_specific_abstain"]
        row["target_transfer_status"] = "clean_state_specific_abstain_missing_state_context"
    row["limitations"] = [
        "condition_or_oligomerization_context_later",
        "state_specific_external_couplings_later",
        "no_fold_switch_claim_without_state_specific_support",
    ]
    return row


def _target_readout(target_result: dict[str, Any], role_target: dict[str, Any]) -> dict[str, Any]:
    tid = str(target_result.get("target_id"))
    if tid == "p53_TAD_MDM2":
        return _p53_readout(target_result, role_target)
    if tid == "KcsA":
        return _kcsa_readout(target_result, role_target)
    if tid == "XCL1_lymphotactin":
        return _xcl1_readout(target_result, role_target)
    row = _base_target(target_result, role_target)
    row["target_transfer_status"] = "unknown_target_clean_abstain"
    row["selected_core_or_clean_abstain"] = "clean_abstain_unknown_target"
    return row


def build_readout(role_manifest: dict[str, Any], preflight: dict[str, Any]) -> dict[str, Any]:
    role_by_id = _role_manifest_by_id(role_manifest)
    target_results = preflight.get("target_results") if isinstance(preflight.get("target_results"), list) else []
    rows: list[dict[str, Any]] = []
    for target in target_results:
        if not isinstance(target, dict):
            continue
        rows.append(_target_readout(target, role_by_id.get(str(target.get("target_id")), {})))

    target_ids = {str(r.get("target_id")) for r in rows}
    blocked = [str(r.get("target_id")) for r in rows if str(r.get("target_transfer_status", "")).startswith("blocked")]
    role_passed = [str(r.get("target_id")) for r in rows if r.get("positive_role_context_found") is True]
    abstain = [str(r.get("target_id")) for r in rows if r.get("positive_role_context_found") is not True]
    violations = {
        str(r.get("target_id")): r.get("forbidden_misclassification_violations", [])
        for r in rows
        if r.get("forbidden_misclassification_violations")
    }

    checks = {
        "source_preflight_ready": preflight.get("data_preflight_status") == "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT",
        "role_manifest_kind_valid": role_manifest.get("kind") == "V16_LOCKED_GRAMMAR_TRANSFER_TARGET_MANIFEST_v0",
        "required_targets_present": REQUIRED_TARGET_IDS.issubset(target_ids),
        "claim_disabled": role_manifest.get("claim_allowed") is False and preflight.get("claim_allowed") is False,
        "no_new_md": preflight.get("new_md_executed") is False,
        "fixed_residue_cutoff_retired": preflight.get("fixed_residue_cutoff_used") is False,
        "native_metrics_not_used_for_selection": preflight.get("native_metrics_not_used_for_selection") is True,
        "forbidden_misclassification_avoided": not violations,
    }
    failed = [k for k, ok in checks.items() if not ok]
    if failed:
        status = "V16_ZERO_MD_ROLE_TRANSFER_BLOCKED_PREFLIGHT_OR_GUARD_FAILURE"
    elif blocked:
        status = "V16_ZERO_MD_ROLE_TRANSFER_BLOCKED_TARGETS_PRESENT"
    else:
        status = "V16_ZERO_MD_ROLE_TRANSFER_READOUT_COMPLETED_CLAIM_DISABLED"

    return {
        "kind": "V16_ZERO_MD_ROLE_TRANSFER_READOUT_v0",
        "run_mode": "zero_md_role_transfer_readout_no_simulation_no_threshold_tuning",
        "role_transfer_status": status,
        "readout_failed_checks": failed,
        "source_preflight_status": preflight.get("data_preflight_status"),
        "source_v16_lock_status": preflight.get("source_v16_lock_status"),
        "claim_allowed": False,
        "new_md_executed": False,
        "data_preflight_executed": True,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "native_metrics_not_used_for_selection": True,
        "selection_policy": "role_context_or_clean_abstain_only_no_native_metric_selection",
        # Naming is intentionally conservative: this zero-MD stage reports role classification / pressure-transfer only, not folding evidence.
        "role_classification_passed_targets": role_passed,
        "pressure_role_transfer_passed_targets": role_passed,
        "positive_folding_evidence_targets": [],
        "clean_abstain_targets": abstain,
        "blocked_targets": blocked,
        "forbidden_misclassification_violations": violations,
        "target_rows": rows,
        "interpretation": (
            "This is a zero-MD role-transfer readout. It checks whether the locked V15/V16 grammar can assign or abstain from pressure-class roles "
            "without soluble-core misclassification, target-specific threshold tuning, native-metric selection, new MD, or folding-solved claims."
        ),
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V16 Zero-MD Role Transfer Readout",
        "",
        "This readout does not run MD, does not tune thresholds, and does not claim solved folding.",
        "",
        f"Status: `{cert.get('role_transfer_status')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        f"Role-classification passed targets: `{cert.get('role_classification_passed_targets')}`",
        f"Positive folding-evidence targets: `{cert.get('positive_folding_evidence_targets')}`",
        f"Clean abstain targets: `{cert.get('clean_abstain_targets')}`",
        "",
        "## Target rows",
    ]
    for row in cert.get("target_rows", []):
        lines.extend([
            f"### {row.get('target_id')}",
            f"- Status: `{row.get('target_transfer_status')}`",
            f"- Selected/abstain: `{row.get('selected_core_or_clean_abstain')}`",
            f"- Role buckets: `{row.get('role_buckets_assigned')}`",
            f"- Monitor-only roles: `{row.get('monitor_only_roles')}`",
            f"- Clean abstain roles: `{row.get('clean_abstain_roles')}`",
            f"- Forbidden violations: `{row.get('forbidden_misclassification_violations')}`",
            f"- Limitations: `{row.get('limitations')}`",
            "",
        ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run V16 zero-MD role-transfer readout from the V16 preflight certificate.")
    parser.add_argument("--role-manifest", default=str(DEFAULT_ROLE_MANIFEST))
    parser.add_argument("--preflight-cert", default=str(DEFAULT_PREFLIGHT_CERT))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    args = parser.parse_args()

    role_manifest = _read_json(Path(args.role_manifest), "V16 role manifest")
    preflight = _read_json(Path(args.preflight_cert), "V16 transfer data preflight certificate")
    cert = build_readout(role_manifest, preflight)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v16_zero_md_role_transfer_readout_certificate.json"
    report_path = out_dir / "V16_ZERO_MD_ROLE_TRANSFER_READOUT_REPORT.md"
    cert_path.write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(report_path, cert)

    print(json.dumps({
        "kind": cert.get("kind"),
        "certificate": str(cert_path),
        "report": str(report_path),
        "role_transfer_status": cert.get("role_transfer_status"),
        "role_classification_passed_targets": cert.get("role_classification_passed_targets"),
        "pressure_role_transfer_passed_targets": cert.get("pressure_role_transfer_passed_targets"),
        "positive_folding_evidence_targets": cert.get("positive_folding_evidence_targets"),
        "clean_abstain_targets": cert.get("clean_abstain_targets"),
        "blocked_targets": cert.get("blocked_targets"),
        "claim_allowed": False,
        "new_md_executed": False,
    }, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
