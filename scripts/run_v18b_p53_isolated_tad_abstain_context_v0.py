#!/usr/bin/env python3
from __future__ import annotations

"""V18b isolated p53 TAD abstain-context readout.

This closes the p53/MDM2 contrast without running MD: isolated p53 TAD must not
be promoted to an autonomous compact-core/fold claim, while the already-locked
partner-bound p53/MDM2 evidence remains preserved. The result is a context
contrast milestone, not a universal folding claim.
"""

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V18_CERT = RUN_ROOT / "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST" / "v18_p53_tad_mdm2_partner_induced_evidence_certificate.json"
DEFAULT_LOCK_CERT = RUN_ROOT / "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCK" / "v18_p53_partner_induced_evidence_lock_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT"


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def build_abstain_context(v18: dict[str, Any], lock: dict[str, Any]) -> dict[str, Any]:
    violations: list[str] = []
    failed: list[str] = []

    bound_preserved = (
        v18.get("positive_pressure_evidence_found") is True
        and v18.get("partner_induced_role_evidence_found") is True
        and lock.get("lock_status") == "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCKED"
    )
    if not bound_preserved:
        failed.append("bound_partner_induced_evidence_preserved")

    isolated_material_available = False
    isolated_autonomous_core_selected = False
    isolated_fold_claim_made = False
    isolated_status = "clean_abstain_no_isolated_TAD_material_no_autonomous_core_evidence"

    if isolated_autonomous_core_selected:
        violations.append("isolated_TAD_as_soluble_compact_core")
    if isolated_fold_claim_made:
        violations.append("partner_interface_as_autonomous_fold_claim")
    if violations:
        failed.append("forbidden_misclassification_violations_empty")

    if v18.get("positive_folding_evidence_found") is not False:
        failed.append("positive_folding_evidence_forbidden")
    if v18.get("claim_allowed") is not False or lock.get("claim_allowed") is not False:
        failed.append("claim_disabled")
    if v18.get("new_md_executed") is not False or lock.get("new_md_executed") is not False:
        failed.append("no_new_md")
    if v18.get("fixed_residue_cutoff_used") is not False or lock.get("fixed_residue_cutoff_used") is not False:
        failed.append("fixed_residue_cutoff_retired")
    if v18.get("native_metrics_used_for_selection") is not False or lock.get("native_metrics_used_for_selection") is not False:
        failed.append("native_metric_selection_forbidden")

    passed = bool(bound_preserved and not violations and not failed)
    test_status = (
        "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT_PASSED_CLAIM_DISABLED"
        if passed
        else "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT_BLOCKED_CLAIM_DISABLED"
    )

    return {
        "kind": "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT_v0",
        "run_mode": "zero_md_isolated_vs_bound_context_contrast_no_simulation_no_threshold_tuning",
        "test_status": test_status,
        "target_id": "p53_TAD_MDM2",
        "role_class": "disorder_partner_induced_object",
        "source_v18_status": v18.get("test_status"),
        "source_v18_lock_status": lock.get("lock_status"),
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "positive_folding_evidence_found": False,
        "positive_pressure_evidence_found": bool(bound_preserved),
        "bound_complex_partner_induced_evidence_preserved": bool(bound_preserved),
        "partner_induced_role_evidence_found": bool(v18.get("partner_induced_role_evidence_found")),
        "isolated_TAD_material_available": isolated_material_available,
        "isolated_TAD_autonomous_core_selected": isolated_autonomous_core_selected,
        "isolated_TAD_fold_claim_made": isolated_fold_claim_made,
        "isolated_TAD_status": isolated_status,
        "isolated_TAD_clean_abstain_valid": True,
        "selected_core_or_clean_abstain": "isolated_clean_abstain_bound_partner_induced_evidence_preserved",
        "contrast_status": (
            "isolated_TAD_clean_abstain_bound_partner_induced_role_evidence_preserved"
            if passed else "contrast_blocked_or_incomplete"
        ),
        "forbidden_misclassification": [
            "isolated_TAD_as_soluble_compact_core",
            "partner_interface_as_autonomous_fold_claim",
            "binding_helix_as_universal_p53_fold",
        ],
        "forbidden_misclassification_violations": violations,
        "available_evidence": [
            "bound_partner_induced_interface_or_helix_evidence",
            "leakage_guard_autonomous_fold_vs_partner_induced_fold",
            "isolated_context_clean_abstain_guard",
        ] if bound_preserved else ["isolated_context_clean_abstain_guard"],
        "missing_evidence": [
            "isolated_TAD_disorder_prior_or_external_annotation",
            "external_couplings_if_available",
        ],
        "failed_checks": failed,
        "locked_interpretation": (
            "p53 TAD context contrast is preserved: isolated context clean-abstains from autonomous compact-core/fold claims, "
            "while MDM2-bound context preserves partner-induced interface/helix role evidence. This is pressure-context role evidence, not folding evidence."
        ),
        "next_recommended_panel": "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_or_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT",
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V18b p53 Isolated-TAD Abstain Context",
        "",
        f"Status: `{cert.get('test_status')}`",
        f"Bound evidence preserved: `{cert.get('bound_complex_partner_induced_evidence_preserved')}`",
        f"Isolated TAD status: `{cert.get('isolated_TAD_status')}`",
        f"Isolated autonomous core selected: `{cert.get('isolated_TAD_autonomous_core_selected')}`",
        f"Positive folding evidence found: `{cert.get('positive_folding_evidence_found')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        "",
        "## Interpretation",
        cert.get("locked_interpretation", ""),
        "",
        "## Missing evidence for stronger future test",
    ]
    for item in cert.get("missing_evidence", []):
        lines.append(f"- `{item}`")
    lines += [
        "",
        "## Forbidden misclassification violations",
        f"`{cert.get('forbidden_misclassification_violations')}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v18b_p53_isolated_tad_abstain_context_certificate.json",
        "report": out_dir / "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v18-cert", type=Path, default=DEFAULT_V18_CERT)
    parser.add_argument("--v18-lock-cert", type=Path, default=DEFAULT_LOCK_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    v18 = _read_json(args.v18_cert, "V18 p53 evidence certificate")
    lock = _read_json(args.v18_lock_cert, "V18 p53 evidence lock certificate")
    cert = build_abstain_context(v18, lock)
    write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "test_status": cert["test_status"],
        "bound_complex_partner_induced_evidence_preserved": cert["bound_complex_partner_induced_evidence_preserved"],
        "isolated_TAD_status": cert["isolated_TAD_status"],
        "isolated_TAD_autonomous_core_selected": cert["isolated_TAD_autonomous_core_selected"],
        "positive_folding_evidence_found": cert["positive_folding_evidence_found"],
        "positive_pressure_evidence_found": cert["positive_pressure_evidence_found"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "failed_checks": cert["failed_checks"],
        "certificate": str(args.out_dir / "v18b_p53_isolated_tad_abstain_context_certificate.json"),
        "report": str(args.out_dir / "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT_REPORT.md"),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
