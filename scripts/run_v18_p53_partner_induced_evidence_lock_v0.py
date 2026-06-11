#!/usr/bin/env python3
from __future__ import annotations

"""Lock the V18 p53/MDM2 partner-induced evidence milestone.

This is a reporting/locking layer only. It consumes the V18 zero-MD
partner-induced evidence certificate and freezes the interpretation: positive
pressure evidence exists in the partner-bound context, but no folding claim or
isolated p53 TAD autonomous-fold claim is allowed.
"""

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V18_CERT = RUN_ROOT / "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST" / "v18_p53_tad_mdm2_partner_induced_evidence_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCK"


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def build_lock(v18: dict[str, Any]) -> dict[str, Any]:
    failed: list[str] = []
    if v18.get("test_status") != "V18_P53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST_PASSED_CLAIM_DISABLED":
        failed.append("v18_partner_induced_test_passed")
    if v18.get("positive_pressure_evidence_found") is not True:
        failed.append("positive_pressure_evidence_found")
    if v18.get("partner_induced_role_evidence_found") is not True:
        failed.append("partner_induced_role_evidence_found")
    if v18.get("positive_folding_evidence_found") is not False:
        failed.append("positive_folding_evidence_forbidden")
    if v18.get("claim_allowed") is not False:
        failed.append("claim_disabled")
    if v18.get("new_md_executed") is not False:
        failed.append("no_new_md")
    if v18.get("fixed_residue_cutoff_used") is not False:
        failed.append("fixed_residue_cutoff_retired")
    if v18.get("native_metrics_used_for_selection") is not False:
        failed.append("native_metric_selection_forbidden")
    if v18.get("forbidden_misclassification_violations") not in ([], None):
        failed.append("no_forbidden_misclassification")
    if v18.get("isolated_TAD_autonomous_fold_status") != "clean_abstain_no_isolated_TAD_material_no_autonomous_fold_claim":
        failed.append("isolated_TAD_autonomous_fold_claim_not_made")

    iface = v18.get("interface_contact_readout") if isinstance(v18.get("interface_contact_readout"), dict) else {}
    helix = v18.get("partner_bound_helix_proxy_readout") if isinstance(v18.get("partner_bound_helix_proxy_readout"), dict) else {}
    if iface.get("interface_contact_evidence_present") is not True:
        failed.append("interface_contact_evidence_present")
    if iface.get("selection_threshold_used") is not False:
        failed.append("interface_probe_not_fixed_selection_threshold")
    if helix.get("selection_threshold_used") is not False:
        failed.append("helix_proxy_not_fixed_selection_threshold")

    lock_status = (
        "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCKED"
        if not failed
        else "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCK_BLOCKED"
    )
    return {
        "kind": "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCK_v0",
        "run_mode": "lock_only_no_new_simulation_no_threshold_tuning",
        "lock_status": lock_status,
        "source_v18_status": v18.get("test_status"),
        "target_id": "p53_TAD_MDM2",
        "role_class": "disorder_partner_induced_object",
        "locked_claim": "positive_pressure_context_evidence_not_folding_claim",
        "positive_pressure_evidence_found": v18.get("positive_pressure_evidence_found") is True,
        "partner_induced_role_evidence_found": v18.get("partner_induced_role_evidence_found") is True,
        "positive_folding_evidence_found": False,
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "autonomous_TAD_fold_claim": "forbidden_not_made",
        "isolated_TAD_autonomous_fold_status": v18.get("isolated_TAD_autonomous_fold_status"),
        "selected_core_or_clean_abstain": v18.get("selected_core_or_clean_abstain"),
        "available_evidence": v18.get("available_evidence", []),
        "missing_evidence": v18.get("missing_evidence", []),
        "interface_contact_readout": iface,
        "partner_bound_helix_proxy_readout": helix,
        "forbidden_misclassification_violations": v18.get("forbidden_misclassification_violations", []),
        "lock_failed_checks": failed,
        "locked_interpretation": (
            "p53/MDM2 partner-bound complex contains zero-MD partner-induced interface/helix role evidence, "
            "while isolated p53 TAD autonomous compact-fold claims remain forbidden. This is pressure evidence, not folding evidence."
        ),
        "next_step": "V18b_P53_ISOLATED_TAD_ABSTAIN_CONTEXT",
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V18 p53 Partner-Induced Evidence Lock",
        "",
        f"Lock status: `{cert.get('lock_status')}`",
        f"Source V18 status: `{cert.get('source_v18_status')}`",
        f"Positive pressure evidence found: `{cert.get('positive_pressure_evidence_found')}`",
        f"Positive folding evidence found: `{cert.get('positive_folding_evidence_found')}`",
        f"Autonomous TAD fold claim: `{cert.get('autonomous_TAD_fold_claim')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        "",
        "## Interpretation",
        cert.get("locked_interpretation", ""),
        "",
        "## Next",
        f"`{cert.get('next_step')}`",
        "",
        "## Failed checks",
        f"`{cert.get('lock_failed_checks')}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v18_p53_partner_induced_evidence_lock_certificate.json",
        "report": out_dir / "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCK_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v18-cert", type=Path, default=DEFAULT_V18_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    v18 = _read_json(args.v18_cert, "V18 p53 evidence certificate")
    cert = build_lock(v18)
    write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "lock_status": cert["lock_status"],
        "positive_pressure_evidence_found": cert["positive_pressure_evidence_found"],
        "positive_folding_evidence_found": cert["positive_folding_evidence_found"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "lock_failed_checks": cert["lock_failed_checks"],
        "certificate": str(args.out_dir / "v18_p53_partner_induced_evidence_lock_certificate.json"),
        "report": str(args.out_dir / "V18_P53_PARTNER_INDUCED_EVIDENCE_LOCK_REPORT.md"),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
