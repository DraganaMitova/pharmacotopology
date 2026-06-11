#!/usr/bin/env python3
from __future__ import annotations

"""V22 external evidence and annotation acquisition panel.

This is a practical acquisition/readiness sprint, not a folding claim and not an MD
run. It consumes the locked V21 pressure-evidence panel and records which external
annotation/coupling layers are already available, which are still missing, and
which target should be selected for the first V23 real external-evidence contrast
test. It intentionally separates source-reference locks from raw local annotation
files and from evolutionary coupling/MSA availability.
"""

import argparse
import glob
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V21_CERT = RUN_ROOT / "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY" / "v21_pressure_evidence_panel_summary_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL"

SOURCE_REFS: dict[str, list[dict[str, str]]] = {
    "p53_TAD_MDM2": [
        {
            "source_id": "DisProt_DP00086",
            "kind": "disorder_annotation_reference",
            "url": "https://disprot.org/DP00086",
            "use_boundary": "reference_lock_for_isolated_TAD_disorder_context_not_raw_API_claim",
        },
        {
            "source_id": "RCSB_1YCR",
            "kind": "partner_bound_interface_structure_reference",
            "url": "https://www.rcsb.org/structure/1YCR",
            "use_boundary": "bound_p53_TAD_MDM2_complex_context_only_no_autonomous_fold_claim",
        },
    ],
    "KcsA": [
        {
            "source_id": "RCSB_1K4C",
            "kind": "KcsA_Fab_ion_filter_structure_reference",
            "url": "https://www.rcsb.org/structure/1K4C",
            "use_boundary": "pore_filter_ion_context_only_no_whole_channel_fold_claim",
        },
        {
            "source_id": "RCSB_1BL8",
            "kind": "KcsA_integral_membrane_potassium_channel_reference",
            "url": "https://www.rcsb.org/structure/1BL8",
            "use_boundary": "membrane_pore_context_only_no_tetramer_or_MD_claim",
        },
        {
            "source_id": "OPM_KcsA",
            "kind": "transmembrane_topology_reference",
            "url": "https://opm.phar.umich.edu/proteins/59",
            "use_boundary": "topology_annotation_reference_only_no_membrane_MD_claim",
        },
    ],
    "XCL1_lymphotactin": [
        {
            "source_id": "RCSB_2HDM",
            "kind": "state_A_context_structure_reference",
            "url": "https://www.rcsb.org/structure/2HDM",
            "use_boundary": "state_A_context_only_no_single_fold_forcing",
        },
        {
            "source_id": "RCSB_2JP1",
            "kind": "state_B_context_structure_reference",
            "url": "https://www.rcsb.org/structure/2JP1",
            "use_boundary": "state_B_context_only_no_mixed_state_pooling",
        },
        {
            "source_id": "RCSB_2N54",
            "kind": "XCL1_metamorphic_interconversion_reference",
            "url": "https://www.rcsb.org/structure/2N54",
            "use_boundary": "metamorphic_two_state_reference_only_no_fold_switch_claim",
        },
        {
            "source_id": "PMC_8017559",
            "kind": "metamorphic_fold_switching_literature_reference",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8017559/",
            "use_boundary": "literature_context_only_state_specific_guard_required",
        },
    ],
}

LOCAL_STRUCTURE_FILES: dict[str, list[tuple[str, Path]]] = {
    "p53_TAD_MDM2": [("p53_TAD_MDM2_complex_structure", REPO_ROOT / "data" / "v16_pressure_targets" / "p53_TAD_MDM2" / "1YCR.pdb")],
    "KcsA": [("KcsA_Fab_complex_structure", REPO_ROOT / "data" / "v16_pressure_targets" / "KcsA" / "1K4C.pdb")],
    "XCL1_lymphotactin": [
        ("XCL1_state_A_chemokine_like_structure", REPO_ROOT / "data" / "v16_pressure_targets" / "XCL1_lymphotactin" / "2HDM.pdb"),
        ("XCL1_state_B_alternative_conformation_structure", REPO_ROOT / "data" / "v16_pressure_targets" / "XCL1_lymphotactin" / "2JP1.pdb"),
    ],
}

RAW_ANNOTATION_PATTERNS: dict[str, list[str]] = {
    "p53_TAD_MDM2": [
        "data/v22_external_annotations/p53_TAD_MDM2/*disprot*.json",
        "data/v22_external_annotations/p53_TAD_MDM2/*p53*TAD*.json",
        "data/v22_external_annotations/p53_TAD_MDM2/*state*label*.json",
    ],
    "KcsA": [
        "data/v22_external_annotations/KcsA/*opm*.json",
        "data/v22_external_annotations/KcsA/*assembly*.json",
        "data/v22_external_annotations/KcsA/*tetramer*.json",
        "data/v22_external_annotations/KcsA/*filter*.json",
    ],
    "XCL1_lymphotactin": [
        "data/v22_external_annotations/XCL1_lymphotactin/*state*.json",
        "data/v22_external_annotations/XCL1_lymphotactin/*condition*.json",
        "data/v22_external_annotations/XCL1_lymphotactin/*monomer*dimer*.json",
    ],
}

COUPLING_PATTERNS: dict[str, list[str]] = {
    "p53_TAD_MDM2": [
        "external_msa/**/*p53*",
        "external_msa/**/*TP53*",
        "data/**/*p53*coupl*",
        "data/**/*TP53*coupl*",
    ],
    "KcsA": [
        "external_msa/**/*kcsa*",
        "external_msa/**/*KcsA*",
        "external_msa/**/*potassium*channel*",
        "data/**/*kcsa*coupl*",
        "data/**/*KcsA*coupl*",
    ],
    "XCL1_lymphotactin": [
        "external_msa/**/*xcl1*",
        "external_msa/**/*XCL1*",
        "external_msa/**/*lymphotactin*",
        "data/**/*xcl1*coupl*",
        "data/**/*XCL1*coupl*",
    ],
}

TARGET_REQUIREMENTS: dict[str, dict[str, Any]] = {
    "p53_TAD_MDM2": {
        "role_class": "disorder_partner_induced_object",
        "required_evidence": [
            "partner_bound_complex_context",
            "bound_interface_or_contact_evidence",
            "isolated_TAD_disorder_annotation_reference",
            "state_labels_isolated_vs_bound",
            "leakage_guard_autonomous_fold_vs_partner_induced_fold",
        ],
        "first_v23_nonblocking_optional": ["external_couplings_if_available", "raw_local_disorder_annotation_file"],
        "forbidden_misclassification": [
            "isolated_TAD_as_soluble_compact_core",
            "partner_interface_as_autonomous_fold_claim",
            "binding_helix_as_universal_p53_fold",
        ],
        "next_evidence_test": "p53_TAD_MDM2_EXTERNAL_EVIDENCE_CONTRAST_TEST",
    },
    "KcsA": {
        "role_class": "membrane_pore_oligomer_object",
        "required_evidence": [
            "KcsA_complex_coordinate_context",
            "pore_selectivity_filter_annotation_or_diagnostic",
            "transmembrane_topology_annotation_reference",
            "biological_tetramer_assembly_annotation_for_tetramer_claim",
            "leakage_guard_against_soluble_core_misread",
        ],
        "first_v23_nonblocking_optional": ["external_couplings_if_available"],
        "forbidden_misclassification": [
            "membrane_pore_as_soluble_single_domain_compact_core",
            "ion_filter_geometry_as_whole_fold_claim",
            "tetramer_interface_as_monomer_autonomous_core",
        ],
        "next_evidence_test": "KcsA_EXTERNAL_ANNOTATION_AND_COUPLING_READOUT",
    },
    "XCL1_lymphotactin": {
        "role_class": "metamorphic_switch_object",
        "required_evidence": [
            "state_A_and_state_B_context_structures",
            "state_A_and_state_B_labels",
            "monomer_dimer_context",
            "leakage_guard_preventing_mixed_state_fake_core",
            "state_specific_structural_or_contact_evidence",
        ],
        "first_v23_nonblocking_optional": ["condition_labels_if_available", "state_specific_external_couplings_or_constraints_if_available"],
        "forbidden_misclassification": [
            "forcing_single_canonical_fold",
            "mixing_two_states_into_one_false_core",
            "claiming_fold_switch_without_state_specific_support",
        ],
        "next_evidence_test": "XCL1_STATE_SPECIFIC_EXTERNAL_EVIDENCE_READOUT",
    },
}


def _read_json(path: Path, label: str) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing {label}: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be a JSON object: {path}")
    return data


def _glob_repo(root: Path, patterns: list[str], limit: int = 30) -> list[str]:
    hits: list[str] = []
    for pattern in patterns:
        for hit in glob.glob(str(root / pattern), recursive=True):
            p = Path(hit)
            if p.is_file():
                try:
                    hits.append(str(p.relative_to(root)))
                except ValueError:
                    hits.append(str(p))
    return sorted(set(hits))[:limit]


def _structure_status(target_id: str) -> dict[str, Any]:
    rows = []
    for material_id, path in LOCAL_STRUCTURE_FILES[target_id]:
        rows.append({
            "material_id": material_id,
            "path": str(path),
            "present": path.exists() and path.is_file() and path.stat().st_size > 0,
            "size_bytes": path.stat().st_size if path.exists() and path.is_file() else 0,
        })
    return {
        "all_required_structures_present": all(row["present"] for row in rows),
        "structures": rows,
    }


def _row_by_target(v21: dict[str, Any]) -> dict[str, dict[str, Any]]:
    rows = {}
    for row in v21.get("target_rows", []):
        if isinstance(row, dict) and row.get("target_id"):
            rows[row["target_id"]] = row
    return rows


def build_v22(v21: dict[str, Any], repo_root: Path = REPO_ROOT) -> dict[str, Any]:
    v21_rows = _row_by_target(v21)
    manifest_targets = []
    preflight_targets = []
    coupling_targets = []
    decision_rows = []

    failed_checks: list[str] = []
    if v21.get("summary_status") != "V21_PRESSURE_EVIDENCE_PANEL_SUMMARY_LOCKED":
        failed_checks.append("source_v21_summary_locked")
    if v21.get("positive_folding_evidence_targets") not in ([], None):
        failed_checks.append("source_v21_no_positive_folding_evidence")
    if v21.get("claim_allowed") is not False:
        failed_checks.append("source_v21_claim_disabled")
    if v21.get("new_md_executed") is not False:
        failed_checks.append("source_v21_no_new_md")

    for target_id in ["p53_TAD_MDM2", "KcsA", "XCL1_lymphotactin"]:
        req = TARGET_REQUIREMENTS[target_id]
        v21_row = v21_rows.get(target_id, {})
        structure = _structure_status(target_id)
        raw_annotation_hits = _glob_repo(repo_root, RAW_ANNOTATION_PATTERNS[target_id])
        coupling_hits = _glob_repo(repo_root, COUPLING_PATTERNS[target_id])
        refs = SOURCE_REFS[target_id]

        source_refs_locked = bool(refs)
        prior_pressure_evidence = v21_row.get("positive_pressure_evidence_found") is True
        no_forbidden_prior = not v21_row.get("forbidden_misclassification_violations")
        available_from_v21 = set(v21_row.get("available_evidence") or [])

        acquired_evidence: list[str] = []
        if structure["all_required_structures_present"]:
            acquired_evidence.append("required_context_structures_present")
        if source_refs_locked:
            acquired_evidence.append("external_annotation_source_references_locked")
        if prior_pressure_evidence:
            acquired_evidence.append("prior_zero_md_pressure_evidence_locked")
        if no_forbidden_prior:
            acquired_evidence.append("prior_forbidden_misclassification_guard_clean")
        if raw_annotation_hits:
            acquired_evidence.append("raw_local_external_annotation_files_present")
        if coupling_hits:
            acquired_evidence.append("external_couplings_or_msa_files_present")

        # Target-specific first V23 readiness. Couplings are intentionally not a blocker for
        # the first p53 external-evidence contrast test. They remain missing/optional.
        ready_for_v23 = False
        if target_id == "p53_TAD_MDM2":
            ready_for_v23 = (
                prior_pressure_evidence
                and no_forbidden_prior
                and structure["all_required_structures_present"]
                and source_refs_locked
                and "leakage_guard_autonomous_fold_vs_partner_induced_fold" in available_from_v21
            )
        elif target_id == "KcsA":
            ready_for_v23 = (
                prior_pressure_evidence
                and no_forbidden_prior
                and structure["all_required_structures_present"]
                and source_refs_locked
                and bool(raw_annotation_hits)
                and bool(coupling_hits)
            )
        elif target_id == "XCL1_lymphotactin":
            ready_for_v23 = (
                prior_pressure_evidence
                and no_forbidden_prior
                and structure["all_required_structures_present"]
                and source_refs_locked
                and bool(raw_annotation_hits)
                and bool(coupling_hits)
            )

        missing_evidence = []
        if not structure["all_required_structures_present"]:
            missing_evidence.append("required_context_structures_present")
        if not source_refs_locked:
            missing_evidence.append("external_annotation_source_references_locked")
        if not raw_annotation_hits:
            missing_evidence.append("raw_local_external_annotation_files_if_needed")
        if not coupling_hits:
            missing_evidence.append("external_couplings_or_msa_files_if_available")
        if target_id == "KcsA" and not raw_annotation_hits:
            missing_evidence.append("local_biological_assembly_or_OPM_membrane_annotation_file")
        if target_id == "XCL1_lymphotactin" and not raw_annotation_hits:
            missing_evidence.append("local_state_condition_or_monomer_dimer_annotation_file")

        manifest_targets.append({
            "target_id": target_id,
            "role_class": req["role_class"],
            "required_evidence": req["required_evidence"],
            "first_v23_nonblocking_optional": req["first_v23_nonblocking_optional"],
            "external_annotation_source_references": refs,
            "forbidden_misclassification": req["forbidden_misclassification"],
            "clean_abstain_allowed": True,
            "claim_allowed": False,
            "next_evidence_test": req["next_evidence_test"],
        })
        preflight_targets.append({
            "target_id": target_id,
            "role_class": req["role_class"],
            "structure_preflight": structure,
            "external_annotation_source_references_locked": source_refs_locked,
            "raw_local_external_annotation_files": raw_annotation_hits,
            "prior_pressure_evidence_locked": prior_pressure_evidence,
            "prior_forbidden_misclassification_guard_clean": no_forbidden_prior,
            "available_evidence": sorted(set((v21_row.get("available_evidence") or []) + acquired_evidence)),
            "missing_evidence": sorted(set(missing_evidence + (v21_row.get("missing_evidence") or []))),
            "ready_for_V23": ready_for_v23,
            "claim_allowed": False,
        })
        coupling_targets.append({
            "target_id": target_id,
            "external_couplings_or_msa_files_present": bool(coupling_hits),
            "matching_files": coupling_hits,
            "coupling_status": "present" if coupling_hits else "missing_or_not_yet_acquired",
            "coupling_required_for_selected_V23": False if target_id == "p53_TAD_MDM2" else True,
            "claim_allowed": False,
        })
        decision_rows.append({
            "target_id": target_id,
            "ready_for_V23": ready_for_v23,
            "candidate_next_test": req["next_evidence_test"],
            "decision_reason": (
                "strongest ready external evidence contrast: isolated disorder source-reference lock plus partner-bound interface evidence"
                if target_id == "p53_TAD_MDM2" and ready_for_v23 else
                "deferred until local annotation/coupling acquisition closes remaining gaps"
            ),
        })

    ready_targets = [r["target_id"] for r in decision_rows if r["ready_for_V23"]]
    selected = "p53_TAD_MDM2" if "p53_TAD_MDM2" in ready_targets else (ready_targets[0] if ready_targets else None)
    selected_test = TARGET_REQUIREMENTS[selected]["next_evidence_test"] if selected else None
    if selected is None:
        failed_checks.append("one_target_selected_for_v23")

    decision = {
        "decision_status": "V23_TARGET_SELECTED" if selected else "V23_TARGET_SELECTION_BLOCKED",
        "selected_V23_target": selected,
        "selected_V23_test": selected_test,
        "reason": (
            "strongest ready external evidence contrast: isolated disorder vs partner-induced bound interface; no MD required for first V23 test"
            if selected == "p53_TAD_MDM2" else "no target met V23 readiness"
        ),
        "ready_for_V23_targets": ready_targets,
        "blocked_or_deferred_targets": [r["target_id"] for r in decision_rows if not r["ready_for_V23"]],
        "parallel_md_paths_allowed": False,
        "v23_first_run_mode": "zero_md_external_annotation_evidence_contrast_no_MD_no_threshold_tuning",
        "decision_rows": decision_rows,
        "claim_allowed": False,
    }

    manifest = {
        "kind": "V22_EXTERNAL_EVIDENCE_MANIFEST_v0",
        "claim_allowed": False,
        "new_md_executed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "targets": manifest_targets,
    }
    annotation_preflight = {
        "kind": "V22_ANNOTATION_PREFLIGHT_v0",
        "claim_allowed": False,
        "new_md_executed": False,
        "targets": preflight_targets,
        "ready_for_V23_targets": ready_targets,
    }
    coupling_scan = {
        "kind": "V22_COUPLING_AVAILABILITY_SCAN_v0",
        "claim_allowed": False,
        "new_md_executed": False,
        "targets": coupling_targets,
        "targets_with_couplings_or_msa": [row["target_id"] for row in coupling_targets if row["external_couplings_or_msa_files_present"]],
    }

    status = "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL_LOCKED" if not failed_checks else "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL_BLOCKED"
    cert = {
        "kind": "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL_v0",
        "run_mode": "acquisition_manifest_preflight_coupling_scan_v23_decision_no_simulation_no_threshold_tuning",
        "panel_status": status,
        "source_v21_summary_status": v21.get("summary_status"),
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "positive_pressure_evidence_targets_from_v21": v21.get("positive_pressure_evidence_targets", []),
        "positive_folding_evidence_targets": [],
        "ready_for_V23_targets": ready_targets,
        "selected_V23_target": selected,
        "selected_V23_test": selected_test,
        "blocked_or_deferred_targets": decision["blocked_or_deferred_targets"],
        "annotation_source_refs_locked_targets": [r["target_id"] for r in preflight_targets if r["external_annotation_source_references_locked"]],
        "targets_with_local_raw_annotation_files": [r["target_id"] for r in preflight_targets if r["raw_local_external_annotation_files"]],
        "targets_with_couplings_or_msa": coupling_scan["targets_with_couplings_or_msa"],
        "panel_failed_checks": failed_checks,
        "manifest": manifest,
        "annotation_preflight": annotation_preflight,
        "coupling_availability_scan": coupling_scan,
        "next_target_decision": decision,
        "locked_interpretation": (
            "V22 acquires and locks external annotation references, checks local annotation/coupling availability, "
            "and selects one V23 target. p53/MDM2 is selected because it has the strongest ready external-evidence contrast: "
            "isolated disorder reference context plus partner-bound interface evidence. No MD, threshold tuning, native-metric selection, "
            "or folding claim is allowed."
        ),
    }
    return cert


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V22 External Evidence and Annotation Acquisition Panel",
        "",
        f"Panel status: `{cert.get('panel_status')}`",
        f"Selected V23 target: `{cert.get('selected_V23_target')}`",
        f"Selected V23 test: `{cert.get('selected_V23_test')}`",
        f"Ready for V23: `{cert.get('ready_for_V23_targets')}`",
        f"Positive folding evidence targets: `{cert.get('positive_folding_evidence_targets')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        "",
        "## Locked interpretation",
        cert.get("locked_interpretation", ""),
        "",
        "## Target preflight rows",
    ]
    for row in cert.get("annotation_preflight", {}).get("targets", []):
        lines += [
            f"### {row.get('target_id')}",
            f"- Ready for V23: `{row.get('ready_for_V23')}`",
            f"- Annotation source refs locked: `{row.get('external_annotation_source_references_locked')}`",
            f"- Raw local annotation files: `{row.get('raw_local_external_annotation_files')}`",
            f"- Missing evidence: `{row.get('missing_evidence')}`",
            "",
        ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v22_external_evidence_and_annotation_acquisition_panel_certificate.json",
        "manifest": out_dir / "v22_external_evidence_manifest.json",
        "annotation_preflight": out_dir / "v22_annotation_preflight.json",
        "coupling_availability_scan": out_dir / "v22_coupling_availability_scan.json",
        "next_target_decision": out_dir / "v22_next_target_decision.json",
        "report": out_dir / "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["manifest"].write_text(json.dumps(cert["manifest"], indent=2, sort_keys=True), encoding="utf-8")
    paths["annotation_preflight"].write_text(json.dumps(cert["annotation_preflight"], indent=2, sort_keys=True), encoding="utf-8")
    paths["coupling_availability_scan"].write_text(json.dumps(cert["coupling_availability_scan"], indent=2, sort_keys=True), encoding="utf-8")
    paths["next_target_decision"].write_text(json.dumps(cert["next_target_decision"], indent=2, sort_keys=True), encoding="utf-8")
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description="Run V22 external evidence and annotation acquisition panel")
    parser.add_argument("--v21-cert", type=Path, default=DEFAULT_V21_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    v21 = _read_json(args.v21_cert, "V21 pressure evidence panel summary certificate")
    cert = build_v22(v21, REPO_ROOT)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "panel_status": cert["panel_status"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "selected_V23_target": cert["selected_V23_target"],
        "selected_V23_test": cert["selected_V23_test"],
        "ready_for_V23_targets": cert["ready_for_V23_targets"],
        "positive_folding_evidence_targets": cert["positive_folding_evidence_targets"],
        "certificate": str(paths["certificate"]),
        "manifest": str(paths["manifest"]),
        "preflight": str(paths["annotation_preflight"]),
        "coupling_scan": str(paths["coupling_availability_scan"]),
        "decision": str(paths["next_target_decision"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
