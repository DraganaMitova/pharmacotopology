#!/usr/bin/env python3
from __future__ import annotations

"""Lock the V15 dynamic role grammar panel as a milestone.

This is postprocess-only. It can lock either the earlier partial panel or the
three-protein panel after 4AKE has a machine-readable role artifact.  In all
cases the lock remains claim-disabled.
"""

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V15_CERT = RUN_ROOT / "V15_DYNAMIC_SEPARATION_GRAMMAR_READOUT" / "v15_dynamic_separation_grammar_readout_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCKED"


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing V15 source certificate: {path}\nrun scripts/run_v15_dynamic_separation_grammar_readout.sh first")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"invalid V15 certificate: {path}")
    return data


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V15 Dynamic Role Grammar Panel Locked",
        "",
        "This milestone is claim-disabled. It freezes the current dynamic grammar panel without claiming universal protein folding.",
        "",
        f"Lock status: `{cert['lock_status']}`",
        f"Source global status: `{cert['source_global_status']}`",
        f"Claim allowed: `{cert['claim_allowed']}`",
        "",
        "## Locked interpretation",
        "",
        str(cert.get("locked_interpretation")),
        "",
        "## Locked rows",
    ]
    for row in cert.get("locked_protein_rows", []):
        lines.extend([
            f"### {row.get('protein')}",
            f"- Artifact status: `{row.get('artifact_status')}`",
            f"- Target role: `{row.get('target_role')}`",
            f"- Final status: `{row.get('final_status')}`",
            f"- Claim allowed: `{row.get('claim_allowed')}`",
            "",
        ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Lock V15 dynamic role grammar panel milestone.")
    parser.add_argument("--v15-cert", default=str(DEFAULT_V15_CERT))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    args = parser.parse_args()

    source = _read_json(Path(args.v15_cert))
    rows = source.get("protein_rows") if isinstance(source.get("protein_rows"), list) else []
    checks = source.get("coherence_checks") if isinstance(source.get("coherence_checks"), dict) else {}
    positive = source.get("positive_evidence_proteins") or checks.get("positive_evidence_proteins") or []
    missing = source.get("missing_artifacts") or checks.get("missing_artifacts") or []

    lock_checks = {
        "source_certificate_present": True,
        "fixed_residue_cutoff_retired": source.get("fixed_residue_cutoff_used") is False and checks.get("no_fixed_residue_cutoff_used_anywhere") is True,
        "claim_disabled": source.get("claim_allowed") is False and checks.get("no_claim_allowed_anywhere") is True,
        "1ubq_positive": "1UBQ" in positive,
        "1cll_positive": "1CLL" in positive,
        "4ake_not_required_for_lock": True,
        "4ake_positive_when_present": ("4AKE" in positive) or ("4AKE" in missing) or bool(checks.get("bridge_pending_artifacts")),
    }
    failed = [name for name, ok in lock_checks.items() if not ok]
    all_three_positive = all(name in positive for name in ["4AKE", "1UBQ", "1CLL"]) and not missing
    lock_status = (
        "V15_THREE_PROTEIN_DYNAMIC_GRAMMAR_PANEL_LOCKED"
        if all_three_positive and not failed
        else "V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCKED"
        if not failed
        else "V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCK_BLOCKED"
    )
    locked_claim = (
        "unified_role_aware_evidence_grammar_three_protein_panel_locked_not_folding_solved"
        if all_three_positive
        else "unified_role_aware_evidence_grammar_partial_panel_locked_not_folding_solved"
    )
    locked_interpretation = (
        "V15 demonstrates a unified, dynamic, role-aware evidence grammar across three protein object types: "
        "single-domain compact, multi-domain composite, and domain-hinge closure. It does not claim universal protein folding, "
        "does not claim full hinge/closure recovery for 4AKE or 1CLL, does not use a fixed residue-distance cutoff, and keeps claim_allowed=false."
        if all_three_positive
        else "V15 successfully removed fixed sequence-separation cutoff as a decision rule. Protein evidence is assigned by dynamic role context, not residue-distance cutoff. 4AKE machine-readable bridge remains pending or claim-disabled."
    )

    cert = {
        "kind": "V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCKED_v0",
        "run_mode": "postprocess_only_no_new_simulation",
        "lock_status": lock_status,
        "lock_failed_checks": failed,
        "source_certificate": str(Path(args.v15_cert)),
        "source_global_status": source.get("global_status"),
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "locked_claim": locked_claim,
        "locked_interpretation": locked_interpretation,
        "lock_checks": lock_checks,
        "positive_evidence_proteins": positive,
        "missing_artifacts": missing,
        "grammar_axes": source.get("grammar_axes", []),
        "locked_protein_rows": rows,
    }

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v15_dynamic_role_grammar_panel_locked_certificate.json"
    report_path = out_dir / "V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCKED_REPORT.md"
    cert_path.write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(report_path, cert)

    print(json.dumps({
        "kind": cert["kind"],
        "certificate": str(cert_path),
        "report": str(report_path),
        "lock_status": lock_status,
        "lock_failed_checks": failed,
        "positive_evidence_proteins": positive,
        "missing_artifacts": missing,
        "claim_allowed": False,
    }, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
