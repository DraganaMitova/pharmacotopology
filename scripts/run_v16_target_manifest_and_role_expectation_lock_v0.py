#!/usr/bin/env python3
from __future__ import annotations

"""V16 target manifest and role-expectation lock.

This is a manifest/preflight lock only. It does not download data, does not run
MD, does not change the V15 grammar, and does not tune thresholds.  It freezes
what V16 is allowed to test before any target-specific work happens.
"""

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_MANIFEST = REPO_ROOT / "data" / "v16_locked_grammar_transfer_target_manifest.json"
DEFAULT_V15_LOCK_CERT = (
    RUN_ROOT
    / "V15_DYNAMIC_ROLE_GRAMMAR_PANEL_LOCKED"
    / "v15_dynamic_role_grammar_panel_locked_certificate.json"
)
DEFAULT_OUT_DIR = RUN_ROOT / "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCK"

FORBIDDEN_NUMERIC_THRESHOLD_KEYS = {
    "selection_threshold",
    "frequency_threshold",
    "chemical_threshold",
    "dca_threshold",
    "fixed_sequence_separation_cutoff",
    "min_separation",
    "fixed_residue_cutoff",
}

REQUIRED_TARGET_IDS = {"p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"}


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # pragma: no cover - defensive CLI path
        raise SystemExit(f"invalid JSON in {label}: {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _iter_dicts(value: Any):
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from _iter_dicts(child)
    elif isinstance(value, list):
        for child in value:
            yield from _iter_dicts(child)


def _find_forbidden_thresholds(manifest: dict[str, Any]) -> list[str]:
    found: list[str] = []
    for node in _iter_dicts(manifest.get("targets", [])):
        for key, value in node.items():
            if key in FORBIDDEN_NUMERIC_THRESHOLD_KEYS and isinstance(value, (int, float)):
                found.append(f"{key}={value}")
    return found


def _target_check(target: dict[str, Any]) -> dict[str, Any]:
    required_inputs = target.get("required_inputs")
    allowed_roles = target.get("allowed_evidence_roles")
    forbidden = target.get("forbidden_misclassification")
    states = target.get("states", [])
    return {
        "target_id_present": bool(target.get("target_id")),
        "pressure_class_present": bool(target.get("pressure_class")),
        "expected_role_class_present": bool(target.get("expected_role_class")),
        "allowed_evidence_roles_present": isinstance(allowed_roles, list) and bool(allowed_roles),
        "forbidden_misclassification_present": isinstance(forbidden, list) and bool(forbidden),
        "required_inputs_present": isinstance(required_inputs, list) and bool(required_inputs),
        "clean_abstain_allowed": target.get("clean_abstain_allowed") is True or any(
            isinstance(state, dict) and state.get("clean_abstain_allowed") is True for state in states
        ),
        "claim_disabled": target.get("claim_allowed") is False,
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V16 Target Manifest and Role-Expectation Lock",
        "",
        "This is not a new tuning panel and not a folding-solved claim.",
        "",
        f"Lock status: `{cert.get('lock_status')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        f"Data preflight executed: `{cert.get('data_preflight_executed')}`",
        "",
        "## Locked interpretation",
        "",
        str(cert.get("locked_interpretation")),
        "",
        "## Targets",
    ]
    for target in cert.get("locked_targets", []):
        lines.extend([
            f"### {target.get('target_id')}",
            f"- Pressure class: `{target.get('pressure_class')}`",
            f"- Expected role class: `{target.get('expected_role_class')}`",
            f"- Clean abstain allowed: `{target.get('clean_abstain_allowed')}`",
            f"- Claim allowed: `{target.get('claim_allowed')}`",
            "",
            "Allowed evidence roles:",
        ])
        for role in target.get("allowed_evidence_roles", []):
            lines.append(f"- `{role}`")
        lines.extend(["", "Forbidden misclassification:"])
        for role in target.get("forbidden_misclassification", []):
            lines.append(f"- `{role}`")
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def build_lock(manifest: dict[str, Any], v15_lock: dict[str, Any]) -> dict[str, Any]:
    policy = manifest.get("panel_policy") if isinstance(manifest.get("panel_policy"), dict) else {}
    targets = manifest.get("targets") if isinstance(manifest.get("targets"), list) else []
    target_ids = {str(t.get("target_id")) for t in targets if isinstance(t, dict)}
    target_checks = {str(t.get("target_id")): _target_check(t) for t in targets if isinstance(t, dict)}
    forbidden_thresholds = _find_forbidden_thresholds(manifest)

    lock_checks = {
        "v15_three_protein_panel_locked": v15_lock.get("lock_status") == "V15_THREE_PROTEIN_DYNAMIC_GRAMMAR_PANEL_LOCKED",
        "v15_claim_disabled": v15_lock.get("claim_allowed") is False,
        "manifest_kind_valid": manifest.get("kind") == "V16_LOCKED_GRAMMAR_TRANSFER_TARGET_MANIFEST_v0",
        "manifest_claim_disabled": manifest.get("claim_allowed") is False and policy.get("claim_allowed") is False,
        "not_new_tuning_panel": policy.get("not_new_tuning_panel") is True and policy.get("target_specific_threshold_tuning_allowed") is False,
        "not_folding_solved_claim": policy.get("not_proof_of_solved_folding") is True,
        "no_new_md_in_manifest_lock": policy.get("new_md_allowed_in_manifest_lock") is False,
        "no_download_in_manifest_lock": policy.get("data_download_allowed_in_manifest_lock") is False,
        "grammar_changes_forbidden": policy.get("grammar_changes_allowed") is False,
        "native_metrics_not_used_for_selection": policy.get("native_metrics_not_used_for_selection") is True,
        "fixed_threshold_policy_forbidden": policy.get("fixed_threshold_policy") == "forbidden",
        "fixed_residue_cutoff_retired": policy.get("fixed_residue_cutoff_used") is False,
        "required_targets_present": REQUIRED_TARGET_IDS.issubset(target_ids),
        "no_numeric_target_thresholds": not forbidden_thresholds,
        "all_targets_have_role_expectations": all(all(check.values()) for check in target_checks.values()) and bool(target_checks),
    }
    failed = [key for key, ok in lock_checks.items() if not ok]
    lock_status = "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCKED" if not failed else "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCK_BLOCKED"

    return {
        "kind": "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCK_v0",
        "run_mode": "manifest_lock_only_no_data_preflight_no_new_simulation",
        "lock_status": lock_status,
        "lock_failed_checks": failed,
        "lock_checks": lock_checks,
        "target_checks": target_checks,
        "forbidden_thresholds_found": forbidden_thresholds,
        "source_manifest_kind": manifest.get("kind"),
        "source_manifest_panel_name": manifest.get("panel_name"),
        "source_v15_lock_status": v15_lock.get("lock_status"),
        "source_v15_locked_claim": v15_lock.get("locked_claim"),
        "transfer_panel_claim": "locked_transfer_manifest_only_not_new_tuning_not_folding_solved",
        "locked_interpretation": (
            "V16 tests whether the locked V15 dynamic role grammar can transfer from protein object types "
            "to new pressure regimes: disorder/binding, membrane/pore/oligomer, and fold-switching. "
            "It is not a new tuning panel, does not run MD at lock time, forbids fixed target-specific thresholds, "
            "does not use native metrics for selection, and keeps claim_allowed=false."
        ),
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "data_preflight_executed": False,
        "new_md_executed": False,
        "fixed_threshold_policy": policy.get("fixed_threshold_policy"),
        "threshold_policy": policy.get("threshold_policy"),
        "native_metrics_not_used_for_selection": policy.get("native_metrics_not_used_for_selection"),
        "fixed_residue_cutoff_used": False,
        "locked_targets": targets,
        "shared_transfer_axes": manifest.get("shared_transfer_axes", []),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Lock V16 target manifest and role expectations without running data preflight or MD.")
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--v15-lock-cert", default=str(DEFAULT_V15_LOCK_CERT))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    args = parser.parse_args()

    manifest = _read_json(Path(args.manifest), "V16 target manifest")
    v15_lock = _read_json(Path(args.v15_lock_cert), "V15 lock certificate")
    cert = build_lock(manifest, v15_lock)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v16_target_manifest_and_role_expectation_lock_certificate.json"
    report_path = out_dir / "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCK_REPORT.md"
    cert_path.write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(report_path, cert)

    print(json.dumps({
        "kind": cert.get("kind"),
        "certificate": str(cert_path),
        "report": str(report_path),
        "lock_status": cert.get("lock_status"),
        "lock_failed_checks": cert.get("lock_failed_checks"),
        "locked_targets": [target.get("target_id") for target in cert.get("locked_targets", [])],
        "claim_allowed": False,
        "data_preflight_executed": False,
        "new_md_executed": False,
    }, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
