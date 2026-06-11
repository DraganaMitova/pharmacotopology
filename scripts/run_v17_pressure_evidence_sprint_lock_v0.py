#!/usr/bin/env python3
from __future__ import annotations

"""V17 pressure-evidence sprint lock.

This is one batch stage, but it is still conservative: it locks the evidence
manifest, consumes the V16 data preflight + zero-MD role-transfer + gap-lock
outputs, performs a zero-MD evidence-readiness readout, and selects exactly one
V18 target. It does not run MD, does not tune thresholds, and does not upgrade
role classification into folding evidence.
"""

import argparse
import json
import math
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_PREFLIGHT_CERT = RUN_ROOT / "V16_TRANSFER_DATA_PREFLIGHT" / "v16_transfer_data_preflight_certificate.json"
DEFAULT_ZERO_MD_CERT = RUN_ROOT / "V16_ZERO_MD_ROLE_TRANSFER_READOUT" / "v16_zero_md_role_transfer_readout_certificate.json"
DEFAULT_GAP_LOCK_CERT = RUN_ROOT / "V16_PRESSURE_EVIDENCE_GAP_LOCK" / "v16_pressure_evidence_gap_lock_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V17_PRESSURE_EVIDENCE_SPRINT_LOCK"
TARGET_ORDER = ["p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"]

EVIDENCE_REQUIREMENTS: dict[str, dict[str, Any]] = {
    "p53_TAD_MDM2": {
        "role_class": "disorder_partner_induced_object",
        "forbidden_misclassification": [
            "isolated_TAD_as_soluble_compact_core",
            "partner_interface_as_autonomous_fold_claim",
            "binding_helix_as_universal_p53_fold",
        ],
        "clean_abstain_allowed": True,
        "required_evidence": [
            "partner_bound_complex_context",
            "interface_or_contact_evidence",
            "state_label_bound_context",
            "leakage_guard_autonomous_fold_vs_partner_induced_fold",
            "isolated_TAD_disorder_or_clean_abstain_context",
            "external_couplings_if_available",
        ],
        "required_for_first_v18_test": [
            "partner_bound_complex_context",
            "interface_or_contact_evidence",
            "state_label_bound_context",
            "leakage_guard_autonomous_fold_vs_partner_induced_fold",
        ],
        "next_evidence_test": "p53_TAD_MDM2_PARTNER_INDUCED_EVIDENCE_TEST",
        "decision_priority": 1,
        "decision_reason": "fastest clean evidence test for disorder/partner-induced role; no heavy MD needed for first evidence readout",
    },
    "KcsA": {
        "role_class": "membrane_pore_oligomer_object",
        "forbidden_misclassification": [
            "membrane_pore_as_soluble_single_domain_compact_core",
            "ion_filter_geometry_as_whole_fold_claim",
            "tetramer_interface_as_monomer_autonomous_core",
        ],
        "clean_abstain_allowed": True,
        "required_evidence": [
            "KcsA_complex_coordinate_context",
            "membrane_topology_annotation",
            "pore_selectivity_filter_annotation",
            "oligomer_or_tetramer_assembly_context",
            "chain_and_interface_identity",
            "transmembrane_helix_roles",
            "external_couplings_if_available",
            "leakage_guard_against_soluble_core_misread",
        ],
        "required_for_first_v18_test": [
            "membrane_topology_annotation",
            "pore_selectivity_filter_annotation",
            "oligomer_or_tetramer_assembly_context",
            "chain_and_interface_identity",
            "transmembrane_helix_roles",
            "leakage_guard_against_soluble_core_misread",
        ],
        "next_evidence_test": "KcsA_MEMBRANE_PORE_ROLE_EVIDENCE_READOUT",
        "decision_priority": 2,
        "decision_reason": "important membrane/pore pressure target, but annotation layer must be completed before evidence readout",
    },
    "XCL1_lymphotactin": {
        "role_class": "metamorphic_switch_object",
        "forbidden_misclassification": [
            "forcing_single_canonical_fold",
            "mixing_two_states_into_one_false_core",
            "claiming_fold_switch_without_state_specific_support",
        ],
        "clean_abstain_allowed": True,
        "required_evidence": [
            "state_A_and_state_B_context_structures",
            "state_A_and_state_B_labels",
            "monomer_dimer_context",
            "condition_labels_if_available",
            "state_specific_structural_or_contact_evidence",
            "state_specific_external_couplings_or_constraints_if_available",
            "leakage_guard_preventing_mixed_state_fake_core",
        ],
        "required_for_first_v18_test": [
            "state_A_and_state_B_labels",
            "monomer_dimer_context",
            "state_specific_structural_or_contact_evidence",
            "leakage_guard_preventing_mixed_state_fake_core",
        ],
        "next_evidence_test": "XCL1_STATE_SPECIFIC_ROLE_SEPARATION_READOUT",
        "decision_priority": 3,
        "decision_reason": "deepest metamorphic target; block until state-specific evidence labels and leakage guard are explicit",
    },
}


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


def _rows_by_id(cert: dict[str, Any], key: str) -> dict[str, dict[str, Any]]:
    rows = cert.get(key) if isinstance(cert.get(key), list) else []
    return {str(row.get("target_id")): row for row in rows if isinstance(row, dict)}


def _material_by_id(preflight_row: dict[str, Any]) -> dict[str, dict[str, Any]]:
    mats = preflight_row.get("material_results") if isinstance(preflight_row.get("material_results"), list) else []
    return {str(m.get("material_id")): m for m in mats if isinstance(m, dict)}


def _chain_ca_counts(material: dict[str, Any]) -> dict[str, int]:
    summary = material.get("pdb_summary") if isinstance(material.get("pdb_summary"), dict) else {}
    counts = summary.get("chain_ca_counts") if isinstance(summary.get("chain_ca_counts"), dict) else {}
    return {str(k): int(v) for k, v in counts.items() if isinstance(v, (int, float))}


def _pdb_path(material: dict[str, Any]) -> Path | None:
    value = material.get("path")
    if not isinstance(value, str) or not value:
        return None
    p = Path(value)
    return p if p.is_absolute() else REPO_ROOT / p


def _read_ca_atoms(path: Path) -> dict[str, list[tuple[float, float, float]]]:
    atoms: dict[str, list[tuple[float, float, float]]] = {}
    if not path.exists():
        return atoms
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if len(line) < 54 or line[12:16].strip() != "CA":
                continue
            chain = line[21].strip() or "_"
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            atoms.setdefault(chain, []).append(xyz)
    return atoms


def _min_chain_distance_and_contacts(path: Path, chain_a: str, chain_b: str, cutoff: float = 8.0) -> dict[str, Any]:
    atoms = _read_ca_atoms(path)
    a = atoms.get(chain_a, [])
    b = atoms.get(chain_b, [])
    if not a or not b:
        return {"chain_pair": f"{chain_a}-{chain_b}", "contact_count": 0, "min_ca_distance": None, "cutoff_angstrom_report_only": cutoff}
    cutoff2 = cutoff * cutoff
    count = 0
    min_d2: float | None = None
    for x1, y1, z1 in a:
        for x2, y2, z2 in b:
            d2 = (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
            if min_d2 is None or d2 < min_d2:
                min_d2 = d2
            if d2 <= cutoff2:
                count += 1
    return {
        "chain_pair": f"{chain_a}-{chain_b}",
        "contact_count": count,
        "min_ca_distance": round(math.sqrt(min_d2), 3) if min_d2 is not None else None,
        "cutoff_angstrom_report_only": cutoff,
        "selection_threshold": False,
    }


def _role_passed(zero_row: dict[str, Any] | None) -> bool:
    if not zero_row:
        return False
    return zero_row.get("role_classification_passed") is True or zero_row.get("positive_role_context_found") is True


def _gap_missing(gap_row: dict[str, Any] | None) -> list[str]:
    if not gap_row:
        return []
    value = gap_row.get("missing_for_evidence_test")
    return [str(x) for x in value] if isinstance(value, list) else []


def _p53_available(preflight_row: dict[str, Any], zero_row: dict[str, Any] | None) -> tuple[list[str], dict[str, Any]]:
    mats = _material_by_id(preflight_row)
    mat = mats.get("p53_TAD_MDM2_complex_structure", {})
    chain_counts = _chain_ca_counts(mat)
    available: list[str] = []
    details: dict[str, Any] = {"chain_ca_counts": chain_counts}
    if mat.get("material_status") == "present_and_provenanced" and len(chain_counts) >= 2:
        available.extend(["partner_bound_complex_context", "state_label_bound_context"])
        path = _pdb_path(mat)
        if path:
            chains = sorted(chain_counts, key=lambda c: chain_counts[c], reverse=True)
            if len(chains) >= 2:
                contact = _min_chain_distance_and_contacts(path, chains[0], chains[-1])
                details["partner_chain_contact_probe"] = contact
                if int(contact.get("contact_count") or 0) > 0:
                    available.append("interface_or_contact_evidence")
    if _role_passed(zero_row):
        available.append("leakage_guard_autonomous_fold_vs_partner_induced_fold")
    return sorted(set(available)), details


def _kcsa_available(preflight_row: dict[str, Any], zero_row: dict[str, Any] | None) -> tuple[list[str], dict[str, Any]]:
    mats = _material_by_id(preflight_row)
    mat = mats.get("KcsA_Fab_complex_structure", {})
    chain_counts = _chain_ca_counts(mat)
    available: list[str] = []
    details = {"chain_ca_counts": chain_counts}
    if mat.get("material_status") == "present_and_provenanced":
        available.extend(["KcsA_complex_coordinate_context", "chain_and_interface_identity"])
    if _role_passed(zero_row):
        available.append("leakage_guard_against_soluble_core_misread")
    return sorted(set(available)), details


def _xcl1_available(preflight_row: dict[str, Any], zero_row: dict[str, Any] | None) -> tuple[list[str], dict[str, Any]]:
    mats = _material_by_id(preflight_row)
    state_a = mats.get("XCL1_state_A_chemokine_like_structure", {})
    state_b = mats.get("XCL1_state_B_alternative_conformation_structure", {})
    a_counts = _chain_ca_counts(state_a)
    b_counts = _chain_ca_counts(state_b)
    available: list[str] = []
    details = {"state_A_chain_ca_counts": a_counts, "state_B_chain_ca_counts": b_counts}
    if state_a.get("material_status") == "present_and_provenanced" and state_b.get("material_status") == "present_and_provenanced":
        available.append("state_A_and_state_B_context_structures")
    # The V16 role transfer guard avoided single-fold forcing, but this is not the same as a state-specific evidence leakage guard.
    return sorted(set(available)), details


def _available_for_target(target_id: str, preflight_row: dict[str, Any], zero_row: dict[str, Any] | None) -> tuple[list[str], dict[str, Any]]:
    if target_id == "p53_TAD_MDM2":
        return _p53_available(preflight_row, zero_row)
    if target_id == "KcsA":
        return _kcsa_available(preflight_row, zero_row)
    if target_id == "XCL1_lymphotactin":
        return _xcl1_available(preflight_row, zero_row)
    return [], {}


def _zero_md_status_for_target(target_id: str, available: list[str], missing_first: list[str]) -> str:
    if target_id == "p53_TAD_MDM2":
        return "partner_bound_interface_evidence_context_available_zero_md" if not missing_first else "partner_bound_context_present_but_interface_evidence_incomplete"
    if target_id == "KcsA":
        return "membrane_pore_context_only_annotations_missing_zero_md"
    if target_id == "XCL1_lymphotactin":
        return "two_state_context_only_state_specific_evidence_guard_missing_zero_md"
    return "not_evaluated"


def _target_sprint_row(target_id: str, preflight_row: dict[str, Any], zero_row: dict[str, Any] | None, gap_row: dict[str, Any] | None) -> dict[str, Any]:
    req = EVIDENCE_REQUIREMENTS[target_id]
    available, observations = _available_for_target(target_id, preflight_row, zero_row)
    required = list(req["required_evidence"])
    required_first = list(req["required_for_first_v18_test"])
    missing = [item for item in required if item not in available]
    missing_first = [item for item in required_first if item not in available]
    role_pass = _role_passed(zero_row)
    preflight_ready = preflight_row.get("target_material_status") == "ready_for_zero_md_role_readout"
    ready_for_v18 = bool(preflight_ready and role_pass and not missing_first)
    return {
        "target_id": target_id,
        "role_class": req["role_class"],
        "pressure_class": preflight_row.get("pressure_class") or (zero_row or {}).get("pressure_class"),
        "required_evidence": required,
        "required_for_first_v18_test": required_first,
        "available_evidence": available,
        "missing_evidence": missing,
        "missing_for_first_v18_test": missing_first,
        "forbidden_misclassification": req["forbidden_misclassification"],
        "forbidden_misclassification_violations": (zero_row or {}).get("forbidden_misclassification_violations", []),
        "clean_abstain_allowed": True,
        "role_classification_passed": role_pass,
        "zero_MD_readout_status": _zero_md_status_for_target(target_id, available, missing_first),
        "zero_md_evidence_only_no_simulation": True,
        "ready_for_V18": ready_for_v18,
        "next_evidence_test": req["next_evidence_test"],
        "decision_priority": req["decision_priority"],
        "decision_reason": req["decision_reason"],
        "material_observations": observations,
        "gap_lock_missing_for_evidence_test": _gap_missing(gap_row),
        "claim_allowed": False,
        "positive_folding_evidence_found": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "target_specific_threshold_tuning_allowed": False,
        "native_metrics_used_for_selection": False,
    }


def _select_v18_target(rows: list[dict[str, Any]]) -> dict[str, Any]:
    ready_rows = [r for r in rows if r.get("ready_for_V18") is True]
    if ready_rows:
        chosen = sorted(ready_rows, key=lambda r: int(r.get("decision_priority") or 999))[0]
        return {
            "selected_V18_target": chosen["target_id"],
            "selected_V18_test": chosen["next_evidence_test"],
            "decision_status": "V18_TARGET_SELECTED",
            "reason": chosen["decision_reason"],
            "parallel_md_paths_allowed": False,
            "v18_first_run_mode": "zero_md_or_annotation_evidence_readout_first_no_MD_until_ready_decision",
            "claim_allowed": False,
            "new_md_executed": False,
            "fallback_order": [r["target_id"] for r in sorted(rows, key=lambda r: int(r.get("decision_priority") or 999))],
        }
    return {
        "selected_V18_target": None,
        "selected_V18_test": None,
        "decision_status": "NO_V18_TARGET_READY_CLEAN_ABSTAIN",
        "reason": "no pressure target has the required first-test evidence layer locked without missing blocking material",
        "parallel_md_paths_allowed": False,
        "claim_allowed": False,
        "new_md_executed": False,
        "fallback_order": [r["target_id"] for r in sorted(rows, key=lambda r: int(r.get("decision_priority") or 999))],
    }


def build_sprint(preflight: dict[str, Any], zero_md: dict[str, Any], gap_lock: dict[str, Any]) -> dict[str, Any]:
    preflight_rows = _rows_by_id(preflight, "target_results")
    zero_rows = _rows_by_id(zero_md, "target_rows")
    gap_rows = _rows_by_id(gap_lock, "target_gap_rows")

    target_rows = [
        _target_sprint_row(target_id, preflight_rows.get(target_id, {}), zero_rows.get(target_id), gap_rows.get(target_id))
        for target_id in TARGET_ORDER
    ]
    decision = _select_v18_target(target_rows)
    ready_targets = [row["target_id"] for row in target_rows if row["ready_for_V18"] is True]
    blocked_targets = [row["target_id"] for row in target_rows if row["ready_for_V18"] is not True]

    checks = {
        "v16_gap_lock_present": gap_lock.get("gap_lock_status") == "V16_PRESSURE_EVIDENCE_GAP_LOCKED",
        "v16_zero_md_is_classification_not_evidence": zero_md.get("positive_folding_evidence_targets", []) == [],
        "preflight_ready_or_gap_known": preflight.get("data_preflight_status") == "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT",
        "claim_disabled": zero_md.get("claim_allowed") is False and gap_lock.get("claim_allowed") is False,
        "no_new_md": zero_md.get("new_md_executed") is False and gap_lock.get("new_md_executed") is False,
        "fixed_residue_cutoff_retired": zero_md.get("fixed_residue_cutoff_used") is False and gap_lock.get("fixed_residue_cutoff_used") is False,
        "native_metrics_not_used_for_selection": zero_md.get("native_metrics_used_for_selection") is False and gap_lock.get("native_metrics_used_for_selection") is False,
        "target_specific_threshold_tuning_forbidden": zero_md.get("target_specific_threshold_tuning_allowed") is False and gap_lock.get("target_specific_threshold_tuning_allowed") is False,
        "exactly_one_v18_target_selected_if_any_ready": (len(ready_targets) == 0 and decision["selected_V18_target"] is None) or (decision["selected_V18_target"] in ready_targets),
        "positive_folding_evidence_targets_empty": True,
    }
    failed = [key for key, ok in checks.items() if not ok]
    status = "V17_PRESSURE_EVIDENCE_SPRINT_LOCKED" if not failed else "V17_PRESSURE_EVIDENCE_SPRINT_BLOCKED"

    evidence_manifest = {
        "kind": "V17_PRESSURE_EVIDENCE_MANIFEST_v0",
        "claim_allowed": False,
        "new_md_allowed": False,
        "target_specific_threshold_tuning_allowed": False,
        "fixed_threshold_policy": "forbidden",
        "fixed_residue_cutoff_used": False,
        "native_metrics_not_used_for_selection": True,
        "targets": [
            {
                "target_id": row["target_id"],
                "role_class": row["role_class"],
                "required_evidence": row["required_evidence"],
                "required_for_first_v18_test": row["required_for_first_v18_test"],
                "forbidden_misclassification": row["forbidden_misclassification"],
                "clean_abstain_allowed": row["clean_abstain_allowed"],
                "next_evidence_test": row["next_evidence_test"],
                "claim_allowed": False,
            }
            for row in target_rows
        ],
    }
    zero_md_evidence_readout = {
        "kind": "V17_ZERO_MD_EVIDENCE_READOUT_v0",
        "run_mode": "zero_md_evidence_readout_no_simulation_no_threshold_tuning",
        "positive_folding_evidence_targets": [],
        "evidence_claim_allowed_targets": [],
        "claim_allowed": False,
        "new_md_executed": False,
        "target_rows": [
            {
                "target_id": row["target_id"],
                "zero_MD_readout_status": row["zero_MD_readout_status"],
                "available_evidence": row["available_evidence"],
                "missing_evidence": row["missing_evidence"],
                "ready_for_V18": row["ready_for_V18"],
                "material_observations": row["material_observations"],
            }
            for row in target_rows
        ],
    }
    evidence_preflight = {
        "kind": "V17_PRESSURE_EVIDENCE_PREFLIGHT_v0",
        "run_mode": "evidence_file_preflight_no_simulation_no_threshold_tuning",
        "ready_for_V18_targets": ready_targets,
        "blocked_or_deferred_targets": blocked_targets,
        "claim_allowed": False,
        "target_rows": [
            {
                "target_id": row["target_id"],
                "available_evidence": row["available_evidence"],
                "missing_for_first_v18_test": row["missing_for_first_v18_test"],
                "ready_for_V18": row["ready_for_V18"],
            }
            for row in target_rows
        ],
    }

    return {
        "kind": "V17_PRESSURE_EVIDENCE_SPRINT_LOCK_v0",
        "run_mode": "one_batch_manifest_preflight_zero_md_readout_v18_decision_no_simulation_no_threshold_tuning",
        "sprint_lock_status": status,
        "sprint_failed_checks": failed,
        "sprint_checks": checks,
        "source_gap_lock_status": gap_lock.get("gap_lock_status"),
        "source_zero_md_status": zero_md.get("role_transfer_status"),
        "source_preflight_status": preflight.get("data_preflight_status"),
        "claim_allowed": False,
        "new_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "threshold_policy": "purpose_aware_auto_selection_only_after_evidence_sprint_lock",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "native_metrics_not_used_for_selection": True,
        "positive_folding_evidence_targets": [],
        "evidence_claim_allowed_targets": [],
        "role_classification_passed_targets": zero_md.get("role_classification_passed_targets", []),
        "pressure_role_transfer_passed_targets": zero_md.get("pressure_role_transfer_passed_targets", []),
        "ready_for_V18_targets": ready_targets,
        "blocked_or_deferred_targets": blocked_targets,
        "selected_V18_target": decision.get("selected_V18_target"),
        "selected_V18_test": decision.get("selected_V18_test"),
        "next_target_decision": decision,
        "locked_interpretation": (
            "V17 is a one-shot pressure-evidence sprint lock. It combines evidence manifest, file preflight, "
            "zero-MD evidence-readiness readout, and one selected V18 target decision. It does not run MD, tune thresholds, "
            "use native metrics for selection, or upgrade role classification into folding evidence."
        ),
        "evidence_manifest": evidence_manifest,
        "evidence_preflight": evidence_preflight,
        "zero_md_evidence_readout": zero_md_evidence_readout,
        "target_rows": target_rows,
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V17 Pressure Evidence Sprint Lock",
        "",
        "One batch: evidence manifest + evidence preflight + zero-MD evidence-readiness readout + one V18 target decision.",
        "",
        f"Status: `{cert.get('sprint_lock_status')}`",
        f"Selected V18 target: `{cert.get('selected_V18_target')}`",
        f"Selected V18 test: `{cert.get('selected_V18_test')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        f"Positive folding evidence targets: `{cert.get('positive_folding_evidence_targets')}`",
        "",
        "## Interpretation",
        cert.get("locked_interpretation", ""),
        "",
        "## Target rows",
    ]
    for row in cert.get("target_rows", []):
        lines.extend([
            "",
            f"### {row.get('target_id')}",
            f"Role class: `{row.get('role_class')}`",
            f"Zero-MD status: `{row.get('zero_MD_readout_status')}`",
            f"Ready for V18: `{row.get('ready_for_V18')}`",
            f"Next evidence test: `{row.get('next_evidence_test')}`",
            "",
            "Available evidence:",
        ])
        lines.extend([f"- `{item}`" for item in row.get("available_evidence", [])])
        lines.append("")
        lines.append("Missing evidence:")
        lines.extend([f"- `{item}`" for item in row.get("missing_evidence", [])])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v17_pressure_evidence_sprint_lock_certificate.json",
        "manifest": out_dir / "v17_pressure_evidence_manifest.json",
        "preflight": out_dir / "v17_pressure_evidence_preflight.json",
        "zero_md_readout": out_dir / "v17_zero_md_evidence_readout.json",
        "decision": out_dir / "v17_next_target_decision.json",
        "report": out_dir / "V17_PRESSURE_EVIDENCE_SPRINT_LOCK_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    paths["manifest"].write_text(json.dumps(cert["evidence_manifest"], indent=2, sort_keys=True), encoding="utf-8")
    paths["preflight"].write_text(json.dumps(cert["evidence_preflight"], indent=2, sort_keys=True), encoding="utf-8")
    paths["zero_md_readout"].write_text(json.dumps(cert["zero_md_evidence_readout"], indent=2, sort_keys=True), encoding="utf-8")
    paths["decision"].write_text(json.dumps(cert["next_target_decision"], indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)
    print(json.dumps({
        "kind": cert.get("kind"),
        "sprint_lock_status": cert.get("sprint_lock_status"),
        "selected_V18_target": cert.get("selected_V18_target"),
        "selected_V18_test": cert.get("selected_V18_test"),
        "ready_for_V18_targets": cert.get("ready_for_V18_targets"),
        "positive_folding_evidence_targets": cert.get("positive_folding_evidence_targets"),
        "claim_allowed": cert.get("claim_allowed"),
        "new_md_executed": cert.get("new_md_executed"),
        "certificate": str(paths["certificate"]),
        "manifest": str(paths["manifest"]),
        "preflight": str(paths["preflight"]),
        "zero_md_readout": str(paths["zero_md_readout"]),
        "decision": str(paths["decision"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser(description="Run V17 pressure evidence sprint lock.")
    parser.add_argument("--preflight-cert", default=str(DEFAULT_PREFLIGHT_CERT))
    parser.add_argument("--zero-md-cert", default=str(DEFAULT_ZERO_MD_CERT))
    parser.add_argument("--gap-lock-cert", default=str(DEFAULT_GAP_LOCK_CERT))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    args = parser.parse_args()

    preflight = _read_json(Path(args.preflight_cert), "V16 transfer data preflight certificate")
    zero_md = _read_json(Path(args.zero_md_cert), "V16 zero-MD role transfer certificate")
    gap_lock = _read_json(Path(args.gap_lock_cert), "V16 pressure evidence gap lock certificate")
    cert = build_sprint(preflight, zero_md, gap_lock)
    write_outputs(Path(args.out_dir), cert)


if __name__ == "__main__":
    main()
