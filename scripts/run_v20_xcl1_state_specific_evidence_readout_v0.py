#!/usr/bin/env python3
from __future__ import annotations

"""V20 XCL1 state-specific metamorphic evidence readout.

Zero-MD, no-threshold-tuning readout for XCL1/lymphotactin. The target is a
metamorphic switch object, so the pass condition is not "one fold wins". The
readout keeps state A and state B in separate evidence buckets, forbids mixed
state contact pooling/fake core selection, and permits clean abstain if
state-specific evidence is insufficient.
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
DEFAULT_STATE_A_PDB = REPO_ROOT / "data" / "v16_pressure_targets" / "XCL1_lymphotactin" / "2HDM.pdb"
DEFAULT_STATE_B_PDB = REPO_ROOT / "data" / "v16_pressure_targets" / "XCL1_lymphotactin" / "2JP1.pdb"
DEFAULT_OUT_DIR = RUN_ROOT / "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT"
PROBE_RADII_ANGSTROM = [6.0, 8.0, 10.0]


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _target_row_by_id(cert: dict[str, Any], target_id: str) -> dict[str, Any]:
    for key in ["target_results", "target_rows", "targets"]:
        rows = cert.get(key) if isinstance(cert.get(key), list) else []
        for row in rows:
            if isinstance(row, dict) and row.get("target_id") == target_id:
                return row
    return {}


def _material_by_id(row: dict[str, Any]) -> dict[str, dict[str, Any]]:
    mats = row.get("material_results") if isinstance(row.get("material_results"), list) else []
    return {str(m.get("material_id")): m for m in mats if isinstance(m, dict)}


def _resolve_material_path(row: dict[str, Any], material_id: str, fallback: Path) -> Path:
    mat = _material_by_id(row).get(material_id, {})
    value = mat.get("path")
    if isinstance(value, str) and value:
        p = Path(value)
        return p if p.is_absolute() else REPO_ROOT / p
    return fallback


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _parse_pdb(path: Path) -> dict[str, Any]:
    atoms: dict[str, list[dict[str, Any]]] = {}
    helix_records: list[str] = []
    sheet_records: list[str] = []
    annotation_lines: list[str] = []
    model_count = 0
    current_model = 0
    if not path.exists():
        return {
            "exists": False,
            "model_count": 0,
            "atoms_by_chain": {},
            "helix_record_count": 0,
            "sheet_record_count": 0,
            "annotation_text": "",
        }
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            rec = line[:6].strip()
            if rec == "MODEL":
                model_count += 1
                try:
                    current_model = int(line[10:14].strip())
                except ValueError:
                    current_model = model_count
            elif rec in {"HEADER", "TITLE", "COMPND", "SOURCE", "KEYWDS", "REMARK", "EXPDTA"}:
                annotation_lines.append(line.strip())
            elif rec == "HELIX":
                helix_records.append(line.rstrip("\n"))
            elif rec == "SHEET":
                sheet_records.append(line.rstrip("\n"))
            elif rec == "ATOM" and len(line) >= 54 and line[12:16].strip() == "CA":
                chain = line[21].strip() or "_"
                try:
                    xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    resseq = int(line[22:26])
                except ValueError:
                    continue
                atoms.setdefault(chain, []).append({
                    "model": current_model or 1,
                    "residue_number": resseq,
                    "residue_name": line[17:20].strip().upper(),
                    "xyz": xyz,
                })
    if model_count == 0 and atoms:
        model_count = 1
    return {
        "exists": True,
        "model_count": model_count,
        "atoms_by_chain": atoms,
        "chain_ca_counts": {chain: len(items) for chain, items in atoms.items()},
        "chain_count_with_ca": len(atoms),
        "helix_record_count": len(helix_records),
        "sheet_record_count": len(sheet_records),
        "annotation_text": "\n".join(annotation_lines).upper(),
    }


def _chain_interface_probe(parsed: dict[str, Any]) -> dict[str, Any]:
    atoms = parsed.get("atoms_by_chain") if isinstance(parsed.get("atoms_by_chain"), dict) else {}
    chains = sorted(atoms)
    pair_counts: dict[str, dict[str, int]] = {}
    pair_min: dict[str, float | None] = {}
    for i, ca in enumerate(chains):
        for cb in chains[i + 1:]:
            key = f"{ca}-{cb}"
            pair_counts[key] = {str(r): 0 for r in PROBE_RADII_ANGSTROM}
            min_d: float | None = None
            for a in atoms.get(ca, []):
                for b in atoms.get(cb, []):
                    # Avoid reporting artificial contacts between different NMR models.
                    if a.get("model") != b.get("model"):
                        continue
                    d = _distance(a["xyz"], b["xyz"])
                    if min_d is None or d < min_d:
                        min_d = d
                    for r in PROBE_RADII_ANGSTROM:
                        if d <= r:
                            pair_counts[key][str(r)] += 1
            pair_min[key] = round(min_d, 3) if min_d is not None else None
    present = [pair for pair, counts in pair_counts.items() if sum(1 for v in counts.values() if v > 0) >= 2]
    return {
        "chain_ca_counts": parsed.get("chain_ca_counts", {}),
        "chain_count_with_ca": parsed.get("chain_count_with_ca", 0),
        "interface_pairs_present": present,
        "multi_radius_contact_counts_by_chain_pair": pair_counts,
        "min_ca_distance_by_chain_pair": pair_min,
        "contact_probe_policy": "state_specific_multi_radius_report_only_no_single_fixed_selection_threshold",
        "selection_threshold_used": False,
    }


def _state_readout(label: str, pdb_id: str, path: Path, role_context: str) -> dict[str, Any]:
    parsed = _parse_pdb(path)
    interface = _chain_interface_probe(parsed)
    text = str(parsed.get("annotation_text") or "")
    chain_count = int(parsed.get("chain_count_with_ca") or 0)
    ca_count = sum(int(v) for v in (parsed.get("chain_ca_counts") or {}).values())
    if label == "state_A":
        role_evidence = bool(parsed.get("exists") and ca_count > 0 and chain_count >= 1)
        role = "state_A_chemokine_monomer_support_context"
        expected_context = "chemokine_like_or_monomer_context"
    else:
        role_evidence = bool(parsed.get("exists") and ca_count > 0 and chain_count >= 2 and interface.get("interface_pairs_present"))
        role = "state_B_beta_sandwich_dimer_or_alternative_state_support_context"
        expected_context = "alternative_state_or_dimer_context"
    annotation_keywords = {
        "lymphotactin_or_xcl1_keyword": any(k in text for k in ["LYMPHOTACTIN", "XCL1", "CHEMOKINE"]),
        "dimer_keyword": "DIMER" in text,
        "monomer_keyword": "MONOMER" in text,
        "metamorphic_or_switch_keyword": any(k in text for k in ["METAMORPHIC", "SWITCH", "ALTERNATIVE"]),
    }
    return {
        "state_label": label,
        "pdb_id": pdb_id,
        "path": str(path),
        "role_context": role_context,
        "expected_context": expected_context,
        "state_role_bucket": role,
        "state_role_evidence_found": role_evidence,
        "state_material_present": bool(parsed.get("exists")),
        "model_count": parsed.get("model_count"),
        "chain_ca_counts": parsed.get("chain_ca_counts", {}),
        "chain_count_with_ca": chain_count,
        "helix_record_count": parsed.get("helix_record_count"),
        "sheet_record_count": parsed.get("sheet_record_count"),
        "annotation_keywords": annotation_keywords,
        "chain_interface_readout": interface,
        "state_policy": "state_specific_structure_context_only_no_cross_state_pooling_no_single_native_assumption",
        "selection_threshold_used": False,
    }


def build_evidence(preflight: dict[str, Any], zero_md: dict[str, Any], gap_lock: dict[str, Any], state_a_pdb: Path, state_b_pdb: Path) -> dict[str, Any]:
    target_id = "XCL1_lymphotactin"
    preflight_row = _target_row_by_id(preflight, target_id)
    zero_row = _target_row_by_id(zero_md, target_id)
    gap_row = _target_row_by_id(gap_lock, target_id)
    state_a_path = _resolve_material_path(preflight_row, "XCL1_state_A_chemokine_like_structure", state_a_pdb)
    state_b_path = _resolve_material_path(preflight_row, "XCL1_state_B_alternative_conformation_structure", state_b_pdb)

    state_a = _state_readout("state_A", "2HDM", state_a_path, "chemokine_like_monomer_context")
    state_b = _state_readout("state_B", "2JP1", state_b_path, "alternative_beta_or_dimer_context")

    preflight_ready = preflight.get("data_preflight_status") == "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT"
    role_classification_passed = bool(
        zero_row.get("role_classification_passed") is True
        or zero_row.get("positive_role_context_found") is True
        or target_id in (zero_md.get("role_classification_passed_targets") or [])
        or target_id in (zero_md.get("pressure_role_transfer_passed_targets") or [])
    )
    forbidden_from_zero = zero_row.get("forbidden_misclassification_violations") if isinstance(zero_row.get("forbidden_misclassification_violations"), list) else []
    state_labels_present = bool(state_a["state_material_present"] and state_b["state_material_present"])
    monomer_dimer_context = bool(state_a["chain_count_with_ca"] >= 1 and state_b["chain_count_with_ca"] >= 2)
    state_specific_evidence = bool(state_a["state_role_evidence_found"] and state_b["state_role_evidence_found"])
    state_specific_guard = bool(role_classification_passed and not forbidden_from_zero)

    mixed_state_contact_pooling_used = False
    mixed_state_fake_core_selected = False
    single_fold_claim_made = False
    fold_switch_claim_made = False
    forbidden_violations: list[str] = []
    if single_fold_claim_made:
        forbidden_violations.append("forcing_single_canonical_fold")
    if mixed_state_fake_core_selected or mixed_state_contact_pooling_used:
        forbidden_violations.append("mixing_two_states_into_one_false_core")
    if fold_switch_claim_made:
        forbidden_violations.append("claiming_fold_switch_without_state_specific_support")
    forbidden_violations.extend(str(x) for x in forbidden_from_zero)

    positive_pressure = bool(
        preflight_ready
        and role_classification_passed
        and state_labels_present
        and monomer_dimer_context
        and state_specific_evidence
        and state_specific_guard
        and not forbidden_violations
    )
    test_status = (
        "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED"
        if positive_pressure
        else "V20_XCL1_CLEAN_ABSTAIN_STATE_SPECIFIC_EVIDENCE_INSUFFICIENT_CLAIM_DISABLED"
    )
    missing: list[str] = []
    if not state_labels_present:
        missing.append("state_A_and_state_B_labels")
    if not monomer_dimer_context:
        missing.append("monomer_dimer_context")
    if not state_specific_evidence:
        missing.append("state_specific_structural_or_contact_evidence")
    if not state_specific_guard:
        missing.append("leakage_guard_preventing_mixed_state_fake_core")
    missing.extend(["condition_labels_if_available", "state_specific_external_couplings_or_constraints_if_available"])

    available = [
        item for item, present in {
            "state_A_and_state_B_context_structures": state_labels_present,
            "state_A_and_state_B_labels": state_labels_present,
            "monomer_dimer_context": monomer_dimer_context,
            "state_specific_structural_or_contact_evidence": state_specific_evidence,
            "leakage_guard_preventing_mixed_state_fake_core": state_specific_guard,
        }.items() if present
    ]
    role_buckets = []
    if state_a["state_role_evidence_found"]:
        role_buckets.append("state_A_chemokine_monomer_support")
    if state_b["state_role_evidence_found"]:
        role_buckets.append("state_B_beta_sandwich_dimer_support")
    if state_labels_present:
        role_buckets.append("state_specific_context_separation")
    if state_specific_guard:
        role_buckets.append("mixed_state_fake_core_guard")

    return {
        "kind": "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_v0",
        "run_mode": "zero_md_state_specific_evidence_readout_no_simulation_no_threshold_tuning",
        "test_status": test_status,
        "target_id": target_id,
        "role_class": "metamorphic_switch_object",
        "pressure_class": "metamorphic_fold_switching",
        "source_preflight_status": preflight.get("data_preflight_status"),
        "source_zero_md_status": zero_md.get("role_transfer_status"),
        "source_gap_lock_status": gap_lock.get("gap_lock_status"),
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "positive_folding_evidence_found": False,
        "positive_pressure_evidence_found": positive_pressure,
        "state_specific_role_evidence_found": positive_pressure,
        "state_A_role_evidence_found": bool(state_a["state_role_evidence_found"]),
        "state_B_role_evidence_found": bool(state_b["state_role_evidence_found"]),
        "clean_abstain_valid": not positive_pressure,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "selection_policy": "state_specific_role_readout_only_no_native_metric_selection_no_cross_state_pooling",
        "target_role": "metamorphic_switch_object",
        "state_A": "chemokine_like_monomer_context",
        "state_B": "beta_or_alternative_dimer_context",
        "state_labels_present": state_labels_present,
        "monomer_dimer_context_present": monomer_dimer_context,
        "mixed_state_fake_core_forbidden": True,
        "single_fold_forcing_forbidden": True,
        "mixed_state_contact_pooling_used": mixed_state_contact_pooling_used,
        "mixed_state_pollution": mixed_state_fake_core_selected or mixed_state_contact_pooling_used,
        "mixed_state_fake_core_selected": mixed_state_fake_core_selected,
        "single_fold_claim_made": single_fold_claim_made,
        "fold_switch_claim_made": fold_switch_claim_made,
        "forbidden_misclassification": [
            "forcing_single_canonical_fold",
            "mixing_two_states_into_one_false_core",
            "claiming_fold_switch_without_state_specific_support",
        ],
        "forbidden_misclassification_violations": forbidden_violations,
        "role_buckets_assigned": sorted(set(role_buckets)),
        "selected_core_or_clean_abstain": (
            "state_specific_role_separation_evidence_found_no_mixed_state_fake_core"
            if positive_pressure
            else "clean_abstain_state_specific_evidence_insufficient_no_single_fold_claim"
        ),
        "available_evidence": sorted(set(available)),
        "missing_evidence": sorted(set(missing)),
        "missing_for_folding_or_switch_claim": [
            "state_specific_external_couplings_or_constraints_if_available",
            "condition_labels_if_available",
            "independent_state_specific_validation_before_fold_switch_claim",
        ],
        "state_specific_readouts": {
            "state_A": state_a,
            "state_B": state_b,
        },
        "state_separation_guard": {
            "state_A_state_B_contact_pooling_used": mixed_state_contact_pooling_used,
            "mixed_state_fake_core_selected": mixed_state_fake_core_selected,
            "single_native_state_forced": single_fold_claim_made,
            "guard_policy": "state_specific_buckets_no_cross_state_pooling_no_single_native_assumption",
            "selection_threshold_used": False,
        },
        "limitations": [
            "zero_MD_only_no_fold_switch_dynamics_claim",
            "state_specific_external_couplings_later",
            "condition_or_oligomerization_context_later",
            "no_single_canonical_fold_claim",
        ],
        "next_recommended_panel": "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY",
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V20 XCL1 State-Specific Evidence Readout",
        "",
        "Zero-MD, no threshold tuning, no mixed-state contact pooling. This readout treats XCL1/lymphotactin as a metamorphic switch object and forbids forcing one canonical fold.",
        "",
        f"Status: `{cert.get('test_status')}`",
        f"Positive pressure evidence found: `{cert.get('positive_pressure_evidence_found')}`",
        f"Positive folding evidence found: `{cert.get('positive_folding_evidence_found')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        "",
        "## State evidence",
        f"State A evidence found: `{cert.get('state_A_role_evidence_found')}`",
        f"State B evidence found: `{cert.get('state_B_role_evidence_found')}`",
        f"State labels present: `{cert.get('state_labels_present')}`",
        f"Monomer/dimer context present: `{cert.get('monomer_dimer_context_present')}`",
        "",
        "## Guards",
        f"Mixed-state pollution: `{cert.get('mixed_state_pollution')}`",
        f"Single fold claim made: `{cert.get('single_fold_claim_made')}`",
        f"Forbidden violations: `{cert.get('forbidden_misclassification_violations')}`",
        "",
        "## Missing evidence",
        f"`{cert.get('missing_evidence')}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v20_xcl1_state_specific_evidence_readout_certificate.json",
        "state_a_readout": out_dir / "v20_xcl1_state_A_readout.json",
        "state_b_readout": out_dir / "v20_xcl1_state_B_readout.json",
        "state_separation_guard": out_dir / "v20_xcl1_state_separation_guard.json",
        "report": out_dir / "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    paths["state_a_readout"].write_text(json.dumps(cert.get("state_specific_readouts", {}).get("state_A", {}), indent=2, sort_keys=True), encoding="utf-8")
    paths["state_b_readout"].write_text(json.dumps(cert.get("state_specific_readouts", {}).get("state_B", {}), indent=2, sort_keys=True), encoding="utf-8")
    paths["state_separation_guard"].write_text(json.dumps(cert.get("state_separation_guard", {}), indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight-cert", type=Path, default=DEFAULT_PREFLIGHT_CERT)
    parser.add_argument("--zero-md-cert", type=Path, default=DEFAULT_ZERO_MD_CERT)
    parser.add_argument("--gap-lock-cert", type=Path, default=DEFAULT_GAP_LOCK_CERT)
    parser.add_argument("--state-a-pdb", type=Path, default=DEFAULT_STATE_A_PDB)
    parser.add_argument("--state-b-pdb", type=Path, default=DEFAULT_STATE_B_PDB)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    preflight = _read_json(args.preflight_cert, "V16 transfer preflight certificate")
    zero_md = _read_json(args.zero_md_cert, "V16 zero-MD role transfer certificate")
    gap_lock = _read_json(args.gap_lock_cert, "V16 pressure evidence gap lock certificate")
    cert = build_evidence(preflight, zero_md, gap_lock, args.state_a_pdb, args.state_b_pdb)
    write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "test_status": cert["test_status"],
        "positive_pressure_evidence_found": cert["positive_pressure_evidence_found"],
        "state_specific_role_evidence_found": cert["state_specific_role_evidence_found"],
        "state_A_role_evidence_found": cert["state_A_role_evidence_found"],
        "state_B_role_evidence_found": cert["state_B_role_evidence_found"],
        "positive_folding_evidence_found": cert["positive_folding_evidence_found"],
        "mixed_state_pollution": cert["mixed_state_pollution"],
        "single_fold_claim_made": cert["single_fold_claim_made"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "certificate": str(args.out_dir / "v20_xcl1_state_specific_evidence_readout_certificate.json"),
        "report": str(args.out_dir / "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_REPORT.md"),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
