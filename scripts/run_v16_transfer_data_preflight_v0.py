#!/usr/bin/env python3
from __future__ import annotations

"""V16 transfer data preflight.

This checks whether the locked V16 transfer targets have enough public target/context
material to proceed to a zero-MD role-transfer readout. It does not download data,
does not run MD, does not tune thresholds, and does not permit claims.
"""

import argparse
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_ROLE_MANIFEST = REPO_ROOT / "data" / "v16_locked_grammar_transfer_target_manifest.json"
DEFAULT_MATERIAL_MANIFEST = REPO_ROOT / "data" / "v16_transfer_data_material_manifest.json"
DEFAULT_V16_LOCK_CERT = (
    RUN_ROOT
    / "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCK"
    / "v16_target_manifest_and_role_expectation_lock_certificate.json"
)
DEFAULT_OUT_DIR = RUN_ROOT / "V16_TRANSFER_DATA_PREFLIGHT"
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


def _repo_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path


def _summarize_pdb(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"exists": False, "is_file": False, "ca_atom_count": 0, "chain_ca_counts": {}}
    chain_counts: dict[str, int] = {}
    ca_count = 0
    atom_lines = 0
    hetero_lines = 0
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            record = line[0:6].strip()
            if record == "ATOM":
                atom_lines += 1
            elif record == "HETATM":
                hetero_lines += 1
            else:
                continue
            if len(line) >= 22 and line[12:16].strip() == "CA":
                chain = (line[21].strip() or "_") if len(line) > 21 else "_"
                chain_counts[chain] = chain_counts.get(chain, 0) + 1
                ca_count += 1
    return {
        "exists": True,
        "is_file": path.is_file(),
        "ca_atom_count": ca_count,
        "chain_ca_counts": dict(sorted(chain_counts.items())),
        "chain_count_with_ca": len(chain_counts),
        "atom_line_count": atom_lines,
        "hetero_atom_line_count": hetero_lines,
    }


def _provenance_ok(path: Path, expected_pdb_id: str) -> dict[str, Any]:
    if not path.exists():
        return {
            "exists": False,
            "valid_json": False,
            "pdb_id_matches": False,
            "claim_disabled": False,
            "usage_boundary_safe": False,
        }
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {
            "exists": True,
            "valid_json": False,
            "pdb_id_matches": False,
            "claim_disabled": False,
            "usage_boundary_safe": False,
        }
    usage = str(data.get("usage_boundary", ""))
    return {
        "exists": True,
        "valid_json": isinstance(data, dict),
        "pdb_id_matches": str(data.get("pdb_id", "")).upper() == expected_pdb_id.upper(),
        "claim_disabled": data.get("claim_allowed") is False,
        "usage_boundary_safe": "not_native_selection" in usage and "not_folding_solved_claim" in usage,
        "source": data.get("source"),
        "role_use": data.get("role_use"),
    }


def _material_check(material: dict[str, Any]) -> dict[str, Any]:
    path = _repo_path(str(material.get("path", "")))
    prov_path = _repo_path(str(material.get("provenance_path", "")))
    pdb_id = str(material.get("pdb_id", ""))
    pdb_summary = _summarize_pdb(path)
    prov = _provenance_ok(prov_path, pdb_id)
    min_ca = int(material.get("min_ca_atoms", 0) or 0)
    min_chains = int(material.get("min_chains", 0) or 0)

    checks = {
        "pdb_file_present": pdb_summary.get("exists") is True and pdb_summary.get("is_file") is True,
        "pdb_has_ca_atoms": int(pdb_summary.get("ca_atom_count") or 0) >= min_ca,
        "pdb_has_expected_chain_context": int(pdb_summary.get("chain_count_with_ca") or 0) >= min_chains,
        "provenance_present": prov.get("exists") is True,
        "provenance_valid": prov.get("valid_json") is True,
        "provenance_pdb_id_matches": prov.get("pdb_id_matches") is True,
        "provenance_claim_disabled": prov.get("claim_disabled") is True,
        "provenance_usage_boundary_safe": prov.get("usage_boundary_safe") is True,
    }
    failed = [key for key, ok in checks.items() if not ok]
    return {
        "material_id": material.get("material_id"),
        "kind": material.get("kind"),
        "pdb_id": pdb_id,
        "role_use": material.get("role_use"),
        "path": str(path),
        "provenance_path": str(prov_path),
        "material_status": "present_and_provenanced" if not failed else "missing_or_incomplete",
        "failed_checks": failed,
        "checks": checks,
        "pdb_summary": pdb_summary,
        "provenance_summary": prov,
        "material_sanity_floor_not_selection_threshold": {
            "min_ca_atoms": min_ca,
            "min_chains": min_chains,
            "selection_threshold": False,
        },
    }


def _target_role_status(target_id: str, all_ok: bool, material_results: list[dict[str, Any]]) -> str:
    if not all_ok:
        return "blocked_missing_required_material"
    if target_id == "p53_TAD_MDM2":
        return "partner_complex_context_present_isolated_disorder_material_later_no_claim"
    if target_id == "KcsA":
        return "membrane_pore_oligomer_context_present_membrane_annotation_later_no_claim"
    if target_id == "XCL1_lymphotactin":
        return "two_state_metamorphic_context_present_state_specific_readout_later_no_claim"
    return "material_context_present_no_claim"


def build_preflight(role_manifest: dict[str, Any], material_manifest: dict[str, Any], v16_lock: dict[str, Any]) -> dict[str, Any]:
    role_targets = role_manifest.get("targets") if isinstance(role_manifest.get("targets"), list) else []
    role_by_id = {str(t.get("target_id")): t for t in role_targets if isinstance(t, dict)}
    material_targets = material_manifest.get("targets") if isinstance(material_manifest.get("targets"), list) else []
    target_results: list[dict[str, Any]] = []

    for target in material_targets:
        if not isinstance(target, dict):
            continue
        target_id = str(target.get("target_id"))
        role_target = role_by_id.get(target_id, {})
        mats = target.get("required_material") if isinstance(target.get("required_material"), list) else []
        material_results = [_material_check(mat) for mat in mats if isinstance(mat, dict)]
        missing_ids = [str(item.get("material_id")) for item in material_results if item.get("material_status") != "present_and_provenanced"]
        all_ok = bool(material_results) and not missing_ids
        target_results.append({
            "target_id": target_id,
            "pressure_class": target.get("pressure_class"),
            "expected_role_class": target.get("expected_role_class"),
            "clean_abstain_allowed": target.get("clean_abstain_allowed") is True,
            "claim_allowed": False,
            "target_material_status": "ready_for_zero_md_role_readout" if all_ok else "blocked_missing_required_material",
            "target_role_preflight_status": _target_role_status(target_id, all_ok, material_results),
            "missing_required_material_ids": missing_ids,
            "optional_or_later_material": target.get("optional_or_later_material", []),
            "allowed_evidence_roles": role_target.get("allowed_evidence_roles", []),
            "forbidden_misclassification": role_target.get("forbidden_misclassification", []),
            "material_results": material_results,
        })

    target_ids = {t.get("target_id") for t in target_results}
    ready_targets = [t.get("target_id") for t in target_results if t.get("target_material_status") == "ready_for_zero_md_role_readout"]
    blocked_targets = [t.get("target_id") for t in target_results if t.get("target_material_status") != "ready_for_zero_md_role_readout"]

    checks = {
        "v16_manifest_lock_present": v16_lock.get("lock_status") == "V16_TARGET_MANIFEST_AND_ROLE_EXPECTATION_LOCKED",
        "source_manifest_kind_valid": role_manifest.get("kind") == "V16_LOCKED_GRAMMAR_TRANSFER_TARGET_MANIFEST_v0",
        "material_manifest_kind_valid": material_manifest.get("kind") == "V16_TRANSFER_DATA_MATERIAL_MANIFEST_v0",
        "claim_disabled": role_manifest.get("claim_allowed") is False and material_manifest.get("claim_allowed") is False,
        "no_new_md": material_manifest.get("new_md_allowed") is False,
        "grammar_changes_forbidden": material_manifest.get("grammar_changes_allowed") is False,
        "native_metrics_not_used_for_selection": material_manifest.get("native_metrics_not_used_for_selection") is True,
        "fixed_residue_cutoff_retired": material_manifest.get("fixed_residue_cutoff_used") is False,
        "required_targets_present": REQUIRED_TARGET_IDS.issubset({str(x) for x in target_ids}),
    }
    failed = [key for key, ok in checks.items() if not ok]
    if failed:
        status = "V16_TRANSFER_DATA_PREFLIGHT_BLOCKED_MANIFEST_OR_LOCK_FAILURE"
    elif blocked_targets:
        status = "V16_TRANSFER_DATA_PREFLIGHT_BLOCKED_MISSING_REQUIRED_MATERIAL"
    else:
        status = "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT"

    return {
        "kind": "V16_TRANSFER_DATA_PREFLIGHT_v0",
        "run_mode": "data_preflight_only_no_new_simulation_no_threshold_tuning",
        "data_preflight_status": status,
        "preflight_failed_checks": failed,
        "ready_targets": ready_targets,
        "blocked_targets": blocked_targets,
        "missing_material_by_target": {
            str(t.get("target_id")): t.get("missing_required_material_ids", [])
            for t in target_results
            if t.get("missing_required_material_ids")
        },
        "source_v16_lock_status": v16_lock.get("lock_status"),
        "source_v15_lock_status": v16_lock.get("source_v15_lock_status"),
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "data_preflight_executed": True,
        "new_md_executed": False,
        "download_attempted_by_this_script": False,
        "fixed_threshold_policy": material_manifest.get("fixed_threshold_policy"),
        "native_metrics_not_used_for_selection": material_manifest.get("native_metrics_not_used_for_selection"),
        "fixed_residue_cutoff_used": False,
        "target_specific_threshold_tuning_allowed": False,
        "target_results": target_results,
        "interpretation": (
            "This preflight only verifies whether V16 pressure targets have clean target/context material. "
            "It does not select folding evidence, does not run MD, does not tune grammar, and does not permit claims."
        ),
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V16 Transfer Data Preflight",
        "",
        "This is a data/material preflight only. It does not run MD, tune thresholds, or claim solved folding.",
        "",
        f"Status: `{cert.get('data_preflight_status')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        f"Ready targets: `{cert.get('ready_targets')}`",
        f"Blocked targets: `{cert.get('blocked_targets')}`",
        "",
        "## Targets",
    ]
    for target in cert.get("target_results", []):
        lines.extend([
            f"### {target.get('target_id')}",
            f"- Pressure class: `{target.get('pressure_class')}`",
            f"- Expected role class: `{target.get('expected_role_class')}`",
            f"- Target material status: `{target.get('target_material_status')}`",
            f"- Role preflight status: `{target.get('target_role_preflight_status')}`",
            f"- Missing required material: `{target.get('missing_required_material_ids')}`",
            "",
            "Required material:",
        ])
        for mat in target.get("material_results", []):
            lines.extend([
                f"- `{mat.get('material_id')}` / `{mat.get('pdb_id')}`: `{mat.get('material_status')}`",
                f"  - Failed checks: `{mat.get('failed_checks')}`",
                f"  - CA atoms: `{mat.get('pdb_summary', {}).get('ca_atom_count')}`",
                f"  - Chains with CA: `{mat.get('pdb_summary', {}).get('chain_count_with_ca')}`",
            ])
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run V16 transfer data preflight without downloading, MD, or threshold tuning.")
    parser.add_argument("--role-manifest", default=str(DEFAULT_ROLE_MANIFEST))
    parser.add_argument("--material-manifest", default=str(DEFAULT_MATERIAL_MANIFEST))
    parser.add_argument("--v16-lock-cert", default=str(DEFAULT_V16_LOCK_CERT))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    args = parser.parse_args()

    role_manifest = _read_json(Path(args.role_manifest), "V16 role manifest")
    material_manifest = _read_json(Path(args.material_manifest), "V16 material manifest")
    v16_lock = _read_json(Path(args.v16_lock_cert), "V16 manifest lock certificate")
    cert = build_preflight(role_manifest, material_manifest, v16_lock)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v16_transfer_data_preflight_certificate.json"
    report_path = out_dir / "V16_TRANSFER_DATA_PREFLIGHT_REPORT.md"
    cert_path.write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(report_path, cert)

    print(json.dumps({
        "kind": cert.get("kind"),
        "certificate": str(cert_path),
        "report": str(report_path),
        "data_preflight_status": cert.get("data_preflight_status"),
        "ready_targets": cert.get("ready_targets"),
        "blocked_targets": cert.get("blocked_targets"),
        "missing_material_by_target": cert.get("missing_material_by_target"),
        "claim_allowed": False,
        "new_md_executed": False,
    }, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
