#!/usr/bin/env python3
from __future__ import annotations

"""V15 4AKE dynamic grammar bridge.

Postprocess-only.  This does not run MD, does not synthesize steering evidence,
and does not promote visual/input-only artifacts into biological claims.

Purpose: make the 4AKE state explicit for the unified V15 grammar panel.  If a
machine-readable role-aware 4AKE steering certificate exists, the bridge can
summarize it.  If only legacy visual/input material exists, the bridge records a
claim-disabled bridge-pending row instead of pretending that 4AKE is evaluated.
"""

import argparse
import json
from pathlib import Path
from typing import Any, Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_OUT_DIR = RUN_ROOT / "V15_4AKE_DYNAMIC_GRAMMAR_BRIDGE"

DEFAULT_CANDIDATE_ROLE_CERTS = [
    RUN_ROOT / "V13c_4AKE_ADENYLATE_KINASE" / "openmm_tmd_replicas_v0_certificate.json",
    RUN_ROOT / "V13c_4AKE_ADENYLATE_KINASE" / "v15_4ake_dynamic_grammar_bridge_source.json",
    RUN_ROOT / "comparison_v9" / "4AKE" / "garage_role_aware_rescue_selector_v9_certificate.json",
    RUN_ROOT / "comparison_v9" / "4AKE" / "role_aware_rescue_selector_v9_certificate.json",
]

DEFAULT_LEGACY_VISUAL_DIR = RUN_ROOT / "real_coordinate_visuals" / "coord_007_pdb_4AKE_A_adenylate_kinase"
DEFAULT_LEGACY_CERTS = [
    RUN_ROOT / "real_coordinate_visual_8_certificate.json",
    RUN_ROOT / "visual_mechanism_12_certificate.json",
    RUN_ROOT / "visual_mechanism_audit_certificate.json",
]
DEFAULT_INPUT_MATERIAL = [
    REPO_ROOT / "data" / "audit_4ake_v13c.json",
    REPO_ROOT / "data" / "independent_contact_sources" / "AF-P69441-F1-model_v4.pdb",
    REPO_ROOT / "data" / "independent_contact_sources" / "AF-P69441-F1-model_v4.provenance.json",
    REPO_ROOT / "external_msa" / "4ake_pfam00406" / "4ake_pfam00406_focus.fasta",
    REPO_ROOT / "external_msa_free_predictors" / "4ake_esmfold_api.esmfold_api_report.json",
    REPO_ROOT / "external_msa_free_predictors" / "tryhard_runs_live" / "4ake_esmfold_api.pdb",
]


def _read_json(path: Path) -> Optional[dict[str, Any]]:
    if not path.exists() or not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if isinstance(data, dict):
        return data
    return None


def _existing(paths: list[Path]) -> list[str]:
    return [str(p) for p in paths if p.exists()]


def _pair_key(pair: Any) -> str:
    if isinstance(pair, str):
        return pair
    if isinstance(pair, (list, tuple)) and len(pair) == 2:
        return f"{int(pair[0])}-{int(pair[1])}"
    return str(pair)


def _pairs_to_keys(value: Any) -> list[str]:
    if not isinstance(value, list):
        return []
    return [_pair_key(pair) for pair in value]


def _normalize_role_cert(path: Path, payload: dict[str, Any]) -> dict[str, Any]:
    """Best-effort normalizer for a future/restored 4AKE machine-readable artifact."""
    selected_strict = _pairs_to_keys(payload.get("selected_strict_scaffold") or [])
    selected_balanced = _pairs_to_keys(payload.get("selected_balanced_core") or [])
    selected_rescue = _pairs_to_keys(payload.get("selected_border_rescue") or [])
    selected_hinge = _pairs_to_keys(payload.get("selected_hinge_or_interdomain") or payload.get("selected_interdomain_hinge") or [])
    selected_local = _pairs_to_keys(payload.get("selected_local_support") or [])
    runtime_pairs = _pairs_to_keys(payload.get("runtime_v10_selected_pairs") or [])

    selected_pairs = sorted(set(selected_strict + selected_balanced + selected_rescue + selected_hinge + selected_local + runtime_pairs))
    positive = bool(selected_pairs)

    dynamic_pair_roles: dict[str, dict[str, Any]] = {}
    for pair in selected_pairs:
        if pair in selected_strict:
            role = "strict_scaffold"
            evidence = "domain_hinge_strict_scaffold_evidence"
        elif pair in selected_hinge:
            role = "domain_hinge_or_interdomain_closure_candidate"
            evidence = "domain_hinge_closure_evidence"
        elif pair in selected_rescue:
            role = "border_rescue_or_closure_support"
            evidence = "domain_hinge_rescue_evidence"
        elif pair in selected_balanced:
            role = "balanced_domain_hinge_core_candidate"
            evidence = "domain_hinge_balanced_core_evidence"
        elif pair in selected_local:
            role = "local_support"
            evidence = "local_support"
        else:
            role = "runtime_role_aware_selected_pair"
            evidence = "role_aware_runtime_evidence"
        dynamic_pair_roles[pair] = {
            "domain_relation": "domain_hinge_context",
            "role_decision": role,
            "evidence_class": evidence,
            "selected": True,
            "separation_filter_applied": False,
            "fixed_residue_cutoff_used": False,
        }

    return {
        "protein": "4AKE",
        "artifact_status": "present_machine_readable_4ake_role_artifact" if positive else "present_machine_readable_4ake_artifact_no_selected_pairs",
        "artifact_path": str(path),
        "target_role": "domain_hinge_closure_object",
        "grammar_policy": "domain_hinge_dynamic_grammar_bridge_from_machine_readable_artifact",
        "topology_policy": payload.get("topology_policy") or payload.get("topology_mode") or "domain_hinge_closure",
        "chemical_policy": payload.get("chemical_policy") or "role_aware_chemical_policy_unknown_or_legacy",
        "separation_policy": "dynamic_contextual_role_assignment_no_fixed_residue_cutoff",
        "fixed_residue_cutoff_used": False,
        "selected_pairs": selected_pairs,
        "selected_strict_scaffold": selected_strict,
        "selected_balanced_core": selected_balanced,
        "selected_hinge_or_interdomain": selected_hinge,
        "selected_border_rescue": selected_rescue,
        "selected_local_support": selected_local,
        "dynamic_pair_roles": dynamic_pair_roles,
        "replica_support": payload.get("support_by_selected_pair") or payload.get("lane_vote_support") or {},
        "noise_added": payload.get("noise_added"),
        "long_range_evidence_polluted": payload.get("long_range_evidence_polluted"),
        "classification_coverage_ratio": payload.get("classification_coverage_ratio"),
        "claim_lock_status": "claim_locked_4ake_bridge_claim_disabled",
        "claim_lock_failed_checks": [],
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "positive_evidence_found": positive,
        "final_status": (
            "domain_hinge_object_machine_readable_role_evidence_bridged;claim_allowed=false"
            if positive else
            "domain_hinge_object_machine_readable_artifact_present_but_no_selected_role_evidence;claim_allowed=false"
        ),
    }


def build_bridge(
    *,
    role_cert_paths: list[Path],
    legacy_visual_dir: Path,
    legacy_cert_paths: list[Path],
    input_material_paths: list[Path],
) -> dict[str, Any]:
    for path in role_cert_paths:
        payload = _read_json(path)
        if payload is not None:
            row = _normalize_role_cert(path, payload)
            break
    else:
        visual_files = sorted([p for p in legacy_visual_dir.glob("*") if p.is_file()]) if legacy_visual_dir.exists() else []
        row = {
            "protein": "4AKE",
            "artifact_status": "present_legacy_visual_and_input_material_bridge_pending_machine_readable_role_artifact" if visual_files or _existing(input_material_paths) else "missing_machine_readable_grammar_artifact",
            "target_role": "domain_hinge_closure_object",
            "grammar_policy": "legacy_visual_input_bridge_pending_role_aware_machine_readable_artifact",
            "topology_policy": "domain_hinge_closure_pending_machine_readable_role_readout",
            "chemical_policy": None,
            "separation_policy": "dynamic_contextual_role_assignment_no_fixed_residue_cutoff",
            "fixed_residue_cutoff_used": False,
            "selected_pairs": [],
            "selected_strict_scaffold": [],
            "selected_balanced_core": [],
            "selected_hinge_or_interdomain": [],
            "selected_local_support": [],
            "dynamic_pair_roles": {},
            "replica_support": {},
            "noise_added": None,
            "long_range_evidence_polluted": None,
            "classification_coverage_ratio": None,
            "legacy_visual_files": [str(p) for p in visual_files],
            "legacy_global_certificates_present": _existing(legacy_cert_paths),
            "input_material_present": _existing(input_material_paths),
            "claim_lock_status": "bridge_pending_machine_readable_role_artifact",
            "claim_lock_failed_checks": ["machine_readable_4AKE_role_artifact_present"],
            "claim_allowed": False,
            "biological_transfer_claim_allowed": False,
            "positive_evidence_found": False,
            "final_status": "4AKE_legacy_visual_input_material_present_but_machine_readable_role_artifact_missing_claim_disabled",
        }

    return {
        "kind": "V15_4AKE_DYNAMIC_GRAMMAR_BRIDGE_v0",
        "run_mode": "postprocess_only_no_new_simulation",
        "claim_allowed": False,
        "biological_transfer_claim_allowed": False,
        "bridge_policy": "do_not_synthesize_4ake_evidence_from_pdb_or_visual_only_artifacts",
        "role_cert_candidates": [str(p) for p in role_cert_paths],
        "legacy_visual_dir": str(legacy_visual_dir),
        "legacy_cert_candidates": [str(p) for p in legacy_cert_paths],
        "input_material_candidates": [str(p) for p in input_material_paths],
        "protein_row": row,
    }


def _write_report(path: Path, bridge: dict[str, Any]) -> None:
    row = bridge["protein_row"]
    lines = [
        "# V15 4AKE Dynamic Grammar Bridge",
        "",
        "This bridge is postprocess-only and claim-disabled.",
        "",
        "It does not convert PDB-only or visual-only material into biological evidence. If a machine-readable role-aware 4AKE artifact is absent, it records bridge-pending status.",
        "",
        f"Artifact status: `{row.get('artifact_status')}`",
        f"Final status: `{row.get('final_status')}`",
        f"Claim allowed: `{row.get('claim_allowed')}`",
        "",
        "## Existing material",
    ]
    for key in ["legacy_visual_files", "legacy_global_certificates_present", "input_material_present"]:
        values = row.get(key) or []
        lines.append(f"### {key}")
        if values:
            for value in values:
                lines.append(f"- `{value}`")
        else:
            lines.append("- none")
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build V15 4AKE dynamic grammar bridge artifact.")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--legacy-visual-dir", default=str(DEFAULT_LEGACY_VISUAL_DIR))
    parser.add_argument("--role-cert", action="append", default=[], help="Optional machine-readable 4AKE role artifact path. Can be repeated.")
    args = parser.parse_args()

    role_paths = [Path(p) for p in args.role_cert] + DEFAULT_CANDIDATE_ROLE_CERTS
    bridge = build_bridge(
        role_cert_paths=role_paths,
        legacy_visual_dir=Path(args.legacy_visual_dir),
        legacy_cert_paths=DEFAULT_LEGACY_CERTS,
        input_material_paths=DEFAULT_INPUT_MATERIAL,
    )

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cert_path = out_dir / "v15_4ake_dynamic_grammar_bridge_certificate.json"
    report_path = out_dir / "V15_4AKE_DYNAMIC_GRAMMAR_BRIDGE_REPORT.md"
    cert_path.write_text(json.dumps(bridge, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(report_path, bridge)

    row = bridge["protein_row"]
    print(json.dumps({
        "kind": bridge["kind"],
        "certificate": str(cert_path),
        "report": str(report_path),
        "artifact_status": row.get("artifact_status"),
        "positive_evidence_found": row.get("positive_evidence_found"),
        "claim_allowed": False,
        "final_status": row.get("final_status"),
    }, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
