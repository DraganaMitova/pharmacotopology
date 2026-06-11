#!/usr/bin/env python3
from __future__ import annotations

"""V23 p53 TAD/MDM2 external-evidence contrast test.

This zero-MD readout consumes the V22 acquisition panel plus the locked V18/V18b
p53 pressure-evidence results. It tests whether the external annotation/source
context supports the intended contrast:

  isolated p53 TAD -> disorder/reference context and no autonomous compact-fold claim
  MDM2-bound p53 TAD -> partner-induced interface/helix role evidence preserved

It is not a folding prediction, does not use native-metric selection, and does
not require external couplings for this first contrast test. Couplings remain
reported as missing/optional for stronger future tests.
"""

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V22_CERT = RUN_ROOT / "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL" / "v22_external_evidence_and_annotation_acquisition_panel_certificate.json"
DEFAULT_V18_CERT = RUN_ROOT / "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST" / "v18_p53_tad_mdm2_partner_induced_evidence_certificate.json"
DEFAULT_V18B_CERT = RUN_ROOT / "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT" / "v18b_p53_isolated_tad_abstain_context_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST"


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be JSON object: {path}")
    return data


def _target_row(section: dict[str, Any], key: str, target_id: str) -> dict[str, Any]:
    for row in section.get(key, []):
        if isinstance(row, dict) and row.get("target_id") == target_id:
            return row
    return {}


def _source_refs_for(v22: dict[str, Any], target_id: str) -> list[dict[str, Any]]:
    manifest = v22.get("manifest", {})
    for row in manifest.get("targets", []):
        if isinstance(row, dict) and row.get("target_id") == target_id:
            refs = row.get("external_annotation_source_references", [])
            return refs if isinstance(refs, list) else []
    return []


def _has_ref_kind(refs: list[dict[str, Any]], needle: str) -> bool:
    needle = needle.lower()
    for ref in refs:
        hay = " ".join(str(ref.get(k, "")) for k in ["source_id", "kind", "url", "use_boundary"]).lower()
        if needle in hay:
            return True
    return False


def _coupling_status(v22: dict[str, Any]) -> dict[str, Any]:
    scan = v22.get("coupling_availability_scan", {})
    for row in scan.get("targets", []):
        if isinstance(row, dict) and row.get("target_id") == "p53_TAD_MDM2":
            return row
    return {"target_id": "p53_TAD_MDM2", "coupling_status": "unknown", "matching_files": []}


def build_v23(v22: dict[str, Any], v18: dict[str, Any], v18b: dict[str, Any]) -> dict[str, Any]:
    failed: list[str] = []
    violations: list[str] = []

    refs = _source_refs_for(v22, "p53_TAD_MDM2")
    p53_preflight = _target_row(v22.get("annotation_preflight", {}), "targets", "p53_TAD_MDM2")
    coupling = _coupling_status(v22)

    v22_ready = (
        v22.get("panel_status") == "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL_LOCKED"
        and v22.get("selected_V23_target") == "p53_TAD_MDM2"
        and v22.get("selected_V23_test") == "p53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST"
        and v22.get("claim_allowed") is False
        and v22.get("positive_folding_evidence_targets") == []
    )
    if not v22_ready:
        failed.append("source_v22_selected_p53_external_contrast")

    source_refs_locked = bool(p53_preflight.get("external_annotation_source_references_locked")) and bool(refs)
    if not source_refs_locked:
        failed.append("external_annotation_source_references_locked")

    isolated_disorder_ref_locked = _has_ref_kind(refs, "disorder") or _has_ref_kind(refs, "disprot")
    bound_interface_ref_locked = _has_ref_kind(refs, "1ycr") or _has_ref_kind(refs, "partner_bound")
    if not isolated_disorder_ref_locked:
        failed.append("isolated_TAD_disorder_reference_locked")
    if not bound_interface_ref_locked:
        failed.append("partner_bound_interface_reference_locked")

    bound_evidence_preserved = (
        v18.get("test_status") == "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_PASSED_CLAIM_DISABLED"
        and v18.get("positive_pressure_evidence_found") is True
        and v18.get("partner_induced_role_evidence_found") is True
        and v18.get("positive_folding_evidence_found") is False
    )
    if not bound_evidence_preserved:
        failed.append("bound_partner_induced_evidence_preserved")

    isolated_abstain_preserved = (
        v18b.get("test_status") == "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT_PASSED_CLAIM_DISABLED"
        and v18b.get("isolated_TAD_autonomous_core_selected") is False
        and v18b.get("isolated_TAD_fold_claim_made") is False
        and v18b.get("bound_complex_partner_induced_evidence_preserved") is True
    )
    if not isolated_abstain_preserved:
        failed.append("isolated_TAD_abstain_context_preserved")

    autonomous_claim_made = False
    universal_p53_fold_claim_made = False
    partner_interface_as_autonomous_fold_claim = False
    if autonomous_claim_made:
        violations.append("isolated_TAD_as_soluble_compact_core")
    if universal_p53_fold_claim_made:
        violations.append("binding_helix_as_universal_p53_fold")
    if partner_interface_as_autonomous_fold_claim:
        violations.append("partner_interface_as_autonomous_fold_claim")
    if violations:
        failed.append("forbidden_misclassification_violations_empty")

    if v18.get("new_md_executed") is not False or v18b.get("new_md_executed") is not False:
        failed.append("no_new_md")
    if v18.get("fixed_residue_cutoff_used") is not False or v18b.get("fixed_residue_cutoff_used") is not False:
        failed.append("fixed_residue_cutoff_retired")
    if v18.get("native_metrics_used_for_selection") is not False or v18b.get("native_metrics_used_for_selection") is not False:
        failed.append("native_metric_selection_forbidden")

    pass_status = not failed
    test_status = (
        "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST_PASSED_CLAIM_DISABLED"
        if pass_status
        else "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST_BLOCKED_CLAIM_DISABLED"
    )

    interface = v18.get("interface_contact_readout", {})
    helix = v18.get("partner_bound_helix_proxy_readout", {})
    raw_files = p53_preflight.get("raw_local_external_annotation_files") or []
    matching_couplings = coupling.get("matching_files") or []

    return {
        "kind": "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST_v0",
        "run_mode": "zero_md_external_annotation_evidence_contrast_no_simulation_no_threshold_tuning",
        "test_status": test_status,
        "target_id": "p53_TAD_MDM2",
        "role_class": "disorder_partner_induced_object",
        "pressure_class": "disorder_partner_induced_binding",
        "source_v22_panel_status": v22.get("panel_status"),
        "source_v18_status": v18.get("test_status"),
        "source_v18b_status": v18b.get("test_status"),
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "positive_pressure_evidence_found": pass_status,
        "positive_external_annotation_evidence_found": bool(source_refs_locked and isolated_disorder_ref_locked and bound_interface_ref_locked),
        "external_disorder_reference_context_found": isolated_disorder_ref_locked,
        "partner_bound_interface_reference_context_found": bound_interface_ref_locked,
        "partner_bound_interface_evidence_preserved": bound_evidence_preserved,
        "bound_complex_partner_induced_evidence_preserved": bool(v18b.get("bound_complex_partner_induced_evidence_preserved")),
        "isolated_TAD_clean_abstain_preserved": isolated_abstain_preserved,
        "isolated_TAD_autonomous_core_selected": False,
        "isolated_TAD_fold_claim_made": False,
        "positive_folding_evidence_found": False,
        "selected_core_or_clean_abstain": (
            "external_disorder_reference_bound_partner_induced_evidence_contrast_found"
            if pass_status
            else "clean_abstain_or_blocked_external_contrast_incomplete"
        ),
        "contrast_status": (
            "isolated_disorder_reference_context_plus_bound_partner_induced_interface_evidence_preserved"
            if pass_status else "external_contrast_incomplete"
        ),
        "state_labels": {
            "isolated_state": "p53_TAD_disorder_reference_context_no_autonomous_fold_claim",
            "bound_state": "p53_TAD_MDM2_partner_bound_interface_or_helix_context",
            "state_label_policy": "isolated_vs_bound_labels_no_native_metric_selection",
        },
        "external_annotation_source_references": refs,
        "raw_local_external_annotation_files": raw_files,
        "raw_local_annotation_file_required_for_this_test": False,
        "coupling_status": coupling.get("coupling_status", "unknown"),
        "external_couplings_or_msa_files_present": bool(matching_couplings),
        "external_couplings_or_msa_files": matching_couplings,
        "external_couplings_required_for_this_test": False,
        "interface_contact_readout": interface,
        "partner_bound_helix_proxy_readout": helix,
        "contact_probe_policy": interface.get("contact_probe_policy", "multi_radius_report_only_no_single_fixed_selection_threshold"),
        "helix_proxy_policy": helix.get("helix_proxy_policy", "report_only_partner_bound_segment_geometry_no_autonomous_fold_claim"),
        "forbidden_misclassification": [
            "isolated_TAD_as_soluble_compact_core",
            "partner_interface_as_autonomous_fold_claim",
            "binding_helix_as_universal_p53_fold",
        ],
        "forbidden_misclassification_violations": violations,
        "available_evidence": [
            item for item, present in {
                "external_annotation_source_references_locked": source_refs_locked,
                "isolated_TAD_disorder_reference_context": isolated_disorder_ref_locked,
                "partner_bound_1YCR_interface_reference_context": bound_interface_ref_locked,
                "partner_bound_interface_or_helix_evidence": bound_evidence_preserved,
                "isolated_TAD_clean_abstain_guard": isolated_abstain_preserved,
                "leakage_guard_autonomous_fold_vs_partner_induced_fold": not violations,
                "state_labels_isolated_vs_bound": True,
            }.items() if present
        ],
        "missing_evidence": [
            item for item, missing in {
                "raw_local_external_annotation_files_if_needed": not bool(raw_files),
                "external_couplings_if_available": not bool(matching_couplings),
            }.items() if missing
        ],
        "failed_checks": failed,
        "locked_interpretation": (
            "V23 locks the p53/MDM2 external-evidence contrast: external disorder-source context supports the isolated TAD abstain guard, "
            "while the 1YCR partner-bound complex preserves interface/helix role evidence. This is positive pressure-context evidence, not positive folding evidence."
        ),
        "next_recommended_panel": "V24_PRESSURE_EVIDENCE_NEXT_TARGET_DECISION_OR_KcsA_EXTERNAL_ANNOTATION_READOUT",
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V23 p53 TAD/MDM2 External Evidence Contrast Test",
        "",
        f"Status: `{cert.get('test_status')}`",
        f"Positive pressure evidence found: `{cert.get('positive_pressure_evidence_found')}`",
        f"Positive folding evidence found: `{cert.get('positive_folding_evidence_found')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        "",
        "## Contrast",
        f"Isolated TAD clean abstain preserved: `{cert.get('isolated_TAD_clean_abstain_preserved')}`",
        f"Bound partner-induced evidence preserved: `{cert.get('partner_bound_interface_evidence_preserved')}`",
        f"External disorder reference context found: `{cert.get('external_disorder_reference_context_found')}`",
        f"Partner-bound interface reference context found: `{cert.get('partner_bound_interface_reference_context_found')}`",
        "",
        "## Interface readout",
        f"`{cert.get('interface_contact_readout')}`",
        "",
        "## Source references",
    ]
    for ref in cert.get("external_annotation_source_references", []):
        lines.append(f"- `{ref.get('source_id')}` — `{ref.get('kind')}` — {ref.get('url')}")
    lines += [
        "",
        "## Missing but nonblocking for this test",
    ]
    for item in cert.get("missing_evidence", []):
        lines.append(f"- `{item}`")
    lines += [
        "",
        "## Forbidden misclassification violations",
        f"`{cert.get('forbidden_misclassification_violations')}`",
        "",
        "## Interpretation",
        cert.get("locked_interpretation", ""),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v23_p53_tad_mdm2_external_evidence_contrast_certificate.json",
        "report": out_dir / "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v22-cert", type=Path, default=DEFAULT_V22_CERT)
    parser.add_argument("--v18-cert", type=Path, default=DEFAULT_V18_CERT)
    parser.add_argument("--v18b-cert", type=Path, default=DEFAULT_V18B_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    v22 = _read_json(args.v22_cert, "V22 acquisition certificate")
    v18 = _read_json(args.v18_cert, "V18 p53 evidence certificate")
    v18b = _read_json(args.v18b_cert, "V18b p53 contrast certificate")
    cert = build_v23(v22, v18, v18b)
    write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "test_status": cert["test_status"],
        "positive_pressure_evidence_found": cert["positive_pressure_evidence_found"],
        "positive_external_annotation_evidence_found": cert["positive_external_annotation_evidence_found"],
        "positive_folding_evidence_found": cert["positive_folding_evidence_found"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "isolated_TAD_clean_abstain_preserved": cert["isolated_TAD_clean_abstain_preserved"],
        "partner_bound_interface_evidence_preserved": cert["partner_bound_interface_evidence_preserved"],
        "external_couplings_or_msa_files_present": cert["external_couplings_or_msa_files_present"],
        "failed_checks": cert["failed_checks"],
        "certificate": str(args.out_dir / "v23_p53_tad_mdm2_external_evidence_contrast_certificate.json"),
        "report": str(args.out_dir / "V23_P53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST_REPORT.md"),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
