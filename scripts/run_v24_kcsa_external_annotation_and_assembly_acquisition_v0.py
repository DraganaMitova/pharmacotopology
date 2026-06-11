#!/usr/bin/env python3
from __future__ import annotations

"""V24 KcsA external annotation and assembly acquisition sprint.

Practical zero-MD acquisition/readout layer for KcsA after the V19 membrane/pore
pressure-evidence win. It preserves the pore/filter, ion-shell, transmembrane
helix, and chain/interface evidence; locks external membrane/assembly annotation
references; scans coupling/MSA availability; and keeps tetramer/whole-channel
fold claims forbidden unless stronger biological assembly evidence is explicitly
available later.
"""

import argparse
import glob
import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run"
DEFAULT_V19_READOUT_CERT = RUN_ROOT / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT" / "v19_kcsa_membrane_pore_evidence_readout_certificate.json"
DEFAULT_V19_LOCK_CERT = RUN_ROOT / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCK" / "v19_kcsa_membrane_pore_evidence_lock_certificate.json"
DEFAULT_V22_CERT = RUN_ROOT / "V22_EXTERNAL_EVIDENCE_AND_ANNOTATION_ACQUISITION_PANEL" / "v22_external_evidence_and_annotation_acquisition_panel_certificate.json"
DEFAULT_OUT_DIR = RUN_ROOT / "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION"

KCSA_SOURCE_REFS = [
    {
        "source_id": "RCSB_1K4C",
        "kind": "KcsA_Fab_ion_filter_structure_reference",
        "url": "https://www.rcsb.org/structure/1K4C",
        "use_boundary": "pore_filter_ion_coordination_context_only_no_whole_channel_fold_claim",
    },
    {
        "source_id": "RCSB_1BL8",
        "kind": "KcsA_integral_membrane_potassium_channel_reference",
        "url": "https://www.rcsb.org/structure/1BL8",
        "use_boundary": "membrane_pore_architecture_reference_only_no_tetramer_solved_claim",
    },
    {
        "source_id": "OPM_KcsA",
        "kind": "transmembrane_topology_reference",
        "url": "https://opm.phar.umich.edu/proteins/59",
        "use_boundary": "membrane_topology_annotation_reference_only_no_membrane_MD_claim",
    },
]

RAW_ANNOTATION_PATTERNS = [
    "data/v22_external_annotations/KcsA/*opm*.json",
    "data/v22_external_annotations/KcsA/*assembly*.json",
    "data/v22_external_annotations/KcsA/*tetramer*.json",
    "data/v22_external_annotations/KcsA/*filter*.json",
    "data/v24_external_annotations/KcsA/*opm*.json",
    "data/v24_external_annotations/KcsA/*assembly*.json",
    "data/v24_external_annotations/KcsA/*tetramer*.json",
    "data/v24_external_annotations/KcsA/*filter*.json",
]

COUPLING_PATTERNS = [
    "external_msa/**/*kcsa*",
    "external_msa/**/*KcsA*",
    "external_msa/**/*potassium*channel*",
    "external_msa/**/*KcsA*",
    "data/**/*kcsa*coupl*",
    "data/**/*KcsA*coupl*",
    "data/**/*potassium*channel*coupl*",
]


def _read_json(path: Path, label: str, required: bool = True) -> dict[str, Any]:
    if not path.exists():
        if required:
            raise SystemExit(f"missing {label}: {path}")
        return {}
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SystemExit(f"{label} must be JSON object: {path}")
    return data


def _glob_repo(patterns: list[str], limit: int = 50) -> list[str]:
    hits: list[str] = []
    for pattern in patterns:
        for hit in glob.glob(str(REPO_ROOT / pattern), recursive=True):
            p = Path(hit)
            if p.is_file():
                try:
                    hits.append(str(p.relative_to(REPO_ROOT)))
                except ValueError:
                    hits.append(str(p))
    return sorted(set(hits))[:limit]


def _v22_kcsa_refs(v22: dict[str, Any]) -> list[dict[str, Any]]:
    manifest = v22.get("manifest") if isinstance(v22.get("manifest"), dict) else {}
    for row in manifest.get("targets", []):
        if isinstance(row, dict) and row.get("target_id") == "KcsA":
            refs = row.get("external_annotation_source_references")
            if isinstance(refs, list) and refs:
                return [r for r in refs if isinstance(r, dict)]
    return KCSA_SOURCE_REFS


def _has_ref(refs: list[dict[str, Any]], needle: str) -> bool:
    needle = needle.lower()
    for ref in refs:
        text = " ".join(str(ref.get(k, "")) for k in ("source_id", "kind", "url", "use_boundary")).lower()
        if needle in text:
            return True
    return False


def _v22_kcsa_couplings(v22: dict[str, Any]) -> list[str]:
    scan = v22.get("coupling_availability_scan") if isinstance(v22.get("coupling_availability_scan"), dict) else {}
    for row in scan.get("targets", []):
        if isinstance(row, dict) and row.get("target_id") == "KcsA":
            files = row.get("matching_files")
            return [str(f) for f in files] if isinstance(files, list) else []
    return []


def _raw_annotation_summary(files: list[str]) -> dict[str, Any]:
    return {
        "raw_local_external_annotation_files_present": bool(files),
        "raw_local_external_annotation_files": files,
        "raw_local_annotation_required_for_V24": False,
        "raw_local_annotation_note": (
            "source references and V19 structural diagnostics are enough for this V24 acquisition layer; "
            "raw local annotation files can strengthen future assembly/contact tests"
        ),
    }


def build_v24(v19_readout: dict[str, Any], v19_lock: dict[str, Any], v22: dict[str, Any] | None = None) -> dict[str, Any]:
    v22 = v22 or {}
    failed: list[str] = []

    source_refs = _v22_kcsa_refs(v22)
    raw_annotation_files = _glob_repo(RAW_ANNOTATION_PATTERNS)
    coupling_files = sorted(set(_v22_kcsa_couplings(v22) + _glob_repo(COUPLING_PATTERNS)))

    v19_passed = v19_readout.get("test_status") == "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED"
    v19_locked = v19_lock.get("lock_status") == "V19_KcsA_MEMBRANE_PORE_EVIDENCE_LOCKED"
    if not v19_passed:
        failed.append("source_v19_membrane_pore_readout_passed")
    if not v19_locked:
        failed.append("source_v19_membrane_pore_lock_present")

    pore_filter_preserved = bool(
        v19_readout.get("pore_filter_annotation_present") is True
        and v19_readout.get("membrane_pore_role_evidence_found") is True
        and v19_readout.get("positive_pressure_evidence_found") is True
    )
    if not pore_filter_preserved:
        failed.append("pore_filter_role_preserved")

    membrane_topology_ref_locked = _has_ref(source_refs, "opm") or _has_ref(source_refs, "transmembrane")
    membrane_topology_preserved = bool(v19_readout.get("transmembrane_role_present") is True and membrane_topology_ref_locked)
    if not membrane_topology_preserved:
        failed.append("membrane_topology_role_preserved")

    assembly_ref_locked = _has_ref(source_refs, "1bl8") or _has_ref(source_refs, "assembly") or _has_ref(source_refs, "tetramer")
    chain_interface_context = bool(v19_readout.get("oligomer_or_chain_interface_context_present") is True)
    biological_assembly_context_locked = bool(assembly_ref_locked and chain_interface_context)
    if not biological_assembly_context_locked:
        failed.append("biological_assembly_context_locked_for_context_only")

    soluble_guard = bool(
        v19_readout.get("soluble_core_misclassification_avoided") is True
        and v19_lock.get("soluble_core_misclassification_avoided") is True
        and v19_readout.get("forbidden_misclassification_violations") in ([], None)
    )
    if not soluble_guard:
        failed.append("soluble_core_misclassification_avoided")

    if v19_readout.get("whole_channel_fold_claim_made") is not False or v19_lock.get("whole_channel_fold_claim_made") is not False:
        failed.append("whole_channel_fold_claim_forbidden")
    if v19_readout.get("tetramer_claim_made") is not False or v19_lock.get("tetramer_claim_made") is not False:
        failed.append("tetramer_claim_forbidden")
    if v19_readout.get("new_md_executed") is not False or v19_lock.get("new_md_executed") is not False:
        failed.append("no_new_md")
    if v19_readout.get("membrane_md_executed") is not False or v19_lock.get("membrane_md_executed") is not False:
        failed.append("no_membrane_md")
    if v19_readout.get("fixed_residue_cutoff_used") is not False or v19_lock.get("fixed_residue_cutoff_used") is not False:
        failed.append("fixed_residue_cutoff_retired")
    if v19_readout.get("native_metrics_used_for_selection") is not False or v19_lock.get("native_metrics_used_for_selection") is not False:
        failed.append("native_metric_selection_forbidden")

    pore_filter_annotation = {
        "kind": "V24_KcsA_PORE_FILTER_ANNOTATION_v0",
        "target_id": "KcsA",
        "pore_filter_annotation_present": pore_filter_preserved,
        "selectivity_filter_role": "pore_filter_core" if pore_filter_preserved else "missing_or_not_preserved",
        "TVGYG_motif_detected": bool(
            ((v19_readout.get("pore_filter_readout") or {}).get("sequence_motif_probe") or {}).get("strong_TVGYG_motif_detected")
        ),
        "K_plus_ion_diagnostic_shell_present": bool(
            ((v19_readout.get("pore_filter_readout") or {}).get("potassium_ion_probe") or {}).get("ion_filter_diagnostic_shell_present")
        ),
        "potassium_ion_count": ((v19_readout.get("pore_filter_readout") or {}).get("potassium_ion_probe") or {}).get("potassium_ion_count"),
        "pore_filter_readout": v19_readout.get("pore_filter_readout", {}),
        "multi_radius_report_only": True,
        "selection_threshold_used": False,
        "claim_allowed": False,
    }

    membrane_topology_annotation = {
        "kind": "V24_KcsA_MEMBRANE_TOPOLOGY_ANNOTATION_v0",
        "target_id": "KcsA",
        "target_role": "membrane_pore_oligomer_object",
        "external_membrane_topology_reference_locked": membrane_topology_ref_locked,
        "source_references": [ref for ref in source_refs if _has_ref([ref], "opm") or _has_ref([ref], "transmembrane") or _has_ref([ref], "1bl8")],
        "transmembrane_helix_scaffold": "present" if v19_readout.get("transmembrane_role_present") is True else "missing_or_not_detected",
        "transmembrane_helix_readout": v19_readout.get("transmembrane_helix_readout", {}),
        "membrane_environment_context_present": v19_readout.get("membrane_environment_context_present") is True,
        "soluble_compact_core_misclassification": False,
        "soluble_core_misclassification_avoided": soluble_guard,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "claim_allowed": False,
    }

    assembly_preflight = {
        "kind": "V24_KcsA_ASSEMBLY_PREFLIGHT_v0",
        "target_id": "KcsA",
        "tetramer_context_annotation_present": biological_assembly_context_locked,
        "biological_assembly_context_locked": biological_assembly_context_locked,
        "biological_assembly_source_reference_locked": assembly_ref_locked,
        "oligomer_or_chain_interface_context_present": chain_interface_context,
        "chain_interface_readout": v19_readout.get("chain_interface_readout", {}),
        "local_biological_assembly_coordinate_or_annotation_file_present": any("assembly" in f.lower() or "tetramer" in f.lower() for f in raw_annotation_files),
        "tetramer_claim_allowed": False,
        "whole_channel_fold_claim_allowed": False,
        "tetramer_claim_made": False,
        "whole_channel_fold_claim_made": False,
        "missing_for_tetramer_claim": [
            item for item, missing in {
                "external_couplings_if_available": not bool(coupling_files),
                "expanded_biological_assembly_coordinate_or_author_defined_assembly_model_for_tetramer_claim": True,
            }.items() if missing
        ],
        "claim_allowed": False,
    }

    external_coupling_availability = {
        "kind": "V24_KcsA_EXTERNAL_COUPLING_AVAILABILITY_v0",
        "target_id": "KcsA",
        "external_couplings_present": bool(coupling_files),
        "external_coupling_or_msa_files": coupling_files,
        "external_couplings_required_for_V24": False,
        "external_couplings_required_for_future_MD_or_contact_test": True,
        "coupling_policy": "do_not_synthesize_couplings_missing_is_not_v24_failure",
        "claim_allowed": False,
    }

    raw_annotation = _raw_annotation_summary(raw_annotation_files)
    manifest = {
        "kind": "V24_KcsA_EXTERNAL_ANNOTATION_MANIFEST_v0",
        "target_id": "KcsA",
        "target_role": "membrane_pore_oligomer_object",
        "source_v19_lock_status": v19_lock.get("lock_status"),
        "source_v19_readout_status": v19_readout.get("test_status"),
        "external_annotation_source_references": source_refs,
        "external_annotation_source_references_locked": bool(source_refs),
        "required_roles": [
            "transmembrane_helix_scaffold",
            "pore_selectivity_filter_core",
            "tetramer_or_chain_interface_context_support",
            "membrane_environment_context",
            "ion_filter_diagnostic_shell",
        ],
        "forbidden_misclassification": [
            "soluble_single_domain_compact_core",
            "autonomous_soluble_fold_claim",
            "whole_channel_fold_solved_claim",
            "tetramer_solved_claim_without_locked_assembly_evidence",
        ],
        "raw_annotation_summary": raw_annotation,
        "claim_allowed": False,
        "new_md_executed": False,
    }

    positive_pressure = not failed
    status = (
        "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION_PASSED_CLAIM_DISABLED"
        if positive_pressure
        else "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION_BLOCKED_CLAIM_DISABLED"
    )

    next_decision = {
        "kind": "V24_KcsA_NEXT_DECISION_v0",
        "decision_status": "V25_READY_FOR_KcsA_COUPLING_AND_INTERFACE_EVIDENCE_READOUT" if positive_pressure else "V25_BLOCKED_PENDING_KcsA_ANNOTATION_GAPS",
        "selected_next_panel": "V25_KcsA_COUPLING_AND_INTERFACE_EVIDENCE_READOUT" if positive_pressure else None,
        "reason": (
            "KcsA external membrane/pore annotation layer is locked; pore/filter, TM scaffold, ion shell, and chain/interface context are preserved; no tetramer/folding claim allowed"
            if positive_pressure else "KcsA external annotation/assembly layer incomplete"
        ),
        "new_md_recommended": False,
        "membrane_md_recommended": False,
        "parallel_md_paths_allowed": False,
        "claim_allowed": False,
    }

    cert = {
        "kind": "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION_v0",
        "run_mode": "zero_md_external_annotation_assembly_acquisition_no_simulation_no_threshold_tuning",
        "test_status": status,
        "target_id": "KcsA",
        "target_role": "membrane_pore_oligomer_object",
        "source_v19_lock_status": v19_lock.get("lock_status"),
        "source_v19_readout_status": v19_readout.get("test_status"),
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "positive_pressure_evidence_found": positive_pressure,
        "positive_folding_evidence_found": False,
        "membrane_pore_role_evidence_found": positive_pressure,
        "pore_filter_annotation_present": pore_filter_preserved,
        "TVGYG_motif_detected": pore_filter_annotation["TVGYG_motif_detected"],
        "K_plus_ion_diagnostic_shell_present": pore_filter_annotation["K_plus_ion_diagnostic_shell_present"],
        "transmembrane_role_present": membrane_topology_preserved,
        "membrane_environment_context_present": v19_readout.get("membrane_environment_context_present") is True,
        "oligomer_or_chain_interface_context_present": chain_interface_context,
        "tetramer_context_annotation_present": biological_assembly_context_locked,
        "biological_assembly_context_locked": biological_assembly_context_locked,
        "tetramer_claim_allowed": False,
        "whole_channel_fold_claim_allowed": False,
        "tetramer_claim_made": False,
        "whole_channel_fold_claim_made": False,
        "soluble_core_misclassification_avoided": soluble_guard,
        "forbidden_misclassification_violations": [],
        "external_couplings_present": bool(coupling_files),
        "external_couplings_required_for_V24": False,
        "external_couplings_required_for_future_MD_or_contact_test": True,
        "available_evidence": [
            item for item, present in {
                "KcsA_external_annotation_source_references_locked": bool(source_refs),
                "pore_selectivity_filter_annotation_or_diagnostic": pore_filter_preserved,
                "TVGYG_selectivity_filter_motif": pore_filter_annotation["TVGYG_motif_detected"],
                "K_plus_ion_diagnostic_shell": pore_filter_annotation["K_plus_ion_diagnostic_shell_present"],
                "transmembrane_helix_scaffold_context": membrane_topology_preserved,
                "membrane_environment_context": v19_readout.get("membrane_environment_context_present") is True,
                "chain_or_interface_context": chain_interface_context,
                "biological_assembly_context_reference_locked_for_context_only": biological_assembly_context_locked,
                "leakage_guard_against_soluble_core_misread": soluble_guard,
            }.items() if present
        ],
        "missing_evidence": [
            item for item, present in {
                "external_couplings_if_available": bool(coupling_files),
                "expanded_biological_assembly_coordinate_or_author_defined_assembly_model_for_tetramer_claim": False,
            }.items() if not present
        ],
        "v24_failed_checks": failed,
        "external_annotation_manifest": manifest,
        "assembly_preflight": assembly_preflight,
        "pore_filter_annotation": pore_filter_annotation,
        "membrane_topology_annotation": membrane_topology_annotation,
        "external_coupling_availability": external_coupling_availability,
        "next_decision": next_decision,
        "locked_interpretation": (
            "V24 locks KcsA external membrane/pore annotation and assembly-context evidence: pore/filter diagnostic, TVGYG motif, potassium-ion shell, "
            "transmembrane helix scaffold, membrane context, and chain/interface context are preserved while soluble-core, tetramer-solved, "
            "whole-channel fold, membrane-MD, fixed-threshold, and native-metric claims remain forbidden. External couplings are not required for V24 "
            "but are marked as required/conditional for future V25 coupling/contact evidence tests."
        ),
    }
    return cert


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    lines = [
        "# V24 KcsA External Annotation and Assembly Acquisition",
        "",
        f"Status: `{cert.get('test_status')}`",
        f"Target role: `{cert.get('target_role')}`",
        f"Positive pressure evidence found: `{cert.get('positive_pressure_evidence_found')}`",
        f"Positive folding evidence found: `{cert.get('positive_folding_evidence_found')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        f"Membrane MD executed: `{cert.get('membrane_md_executed')}`",
        "",
        "## Locked interpretation",
        cert.get("locked_interpretation", ""),
        "",
        "## Evidence",
        f"Available: `{cert.get('available_evidence')}`",
        f"Missing: `{cert.get('missing_evidence')}`",
        f"Failed checks: `{cert.get('v24_failed_checks')}`",
        "",
        "## Next decision",
        f"`{cert.get('next_decision')}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v24_kcsa_external_annotation_and_assembly_acquisition_certificate.json",
        "manifest": out_dir / "v24_kcsa_external_annotation_manifest.json",
        "assembly_preflight": out_dir / "v24_kcsa_assembly_preflight.json",
        "pore_filter_annotation": out_dir / "v24_kcsa_pore_filter_annotation.json",
        "membrane_topology_annotation": out_dir / "v24_kcsa_membrane_topology_annotation.json",
        "external_coupling_availability": out_dir / "v24_kcsa_external_coupling_availability.json",
        "next_decision": out_dir / "v24_kcsa_next_decision.json",
        "report": out_dir / "V24_KcsA_EXTERNAL_ANNOTATION_AND_ASSEMBLY_ACQUISITION_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["manifest"].write_text(json.dumps(cert["external_annotation_manifest"], indent=2, sort_keys=True), encoding="utf-8")
    paths["assembly_preflight"].write_text(json.dumps(cert["assembly_preflight"], indent=2, sort_keys=True), encoding="utf-8")
    paths["pore_filter_annotation"].write_text(json.dumps(cert["pore_filter_annotation"], indent=2, sort_keys=True), encoding="utf-8")
    paths["membrane_topology_annotation"].write_text(json.dumps(cert["membrane_topology_annotation"], indent=2, sort_keys=True), encoding="utf-8")
    paths["external_coupling_availability"].write_text(json.dumps(cert["external_coupling_availability"], indent=2, sort_keys=True), encoding="utf-8")
    paths["next_decision"].write_text(json.dumps(cert["next_decision"], indent=2, sort_keys=True), encoding="utf-8")
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v19-readout-cert", type=Path, default=DEFAULT_V19_READOUT_CERT)
    parser.add_argument("--v19-lock-cert", type=Path, default=DEFAULT_V19_LOCK_CERT)
    parser.add_argument("--v22-cert", type=Path, default=DEFAULT_V22_CERT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    v19_readout = _read_json(args.v19_readout_cert, "V19 KcsA readout certificate")
    v19_lock = _read_json(args.v19_lock_cert, "V19 KcsA lock certificate")
    v22 = _read_json(args.v22_cert, "V22 external evidence acquisition certificate", required=False)
    cert = build_v24(v19_readout, v19_lock, v22)
    paths = write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "test_status": cert["test_status"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "positive_pressure_evidence_found": cert["positive_pressure_evidence_found"],
        "positive_folding_evidence_found": cert["positive_folding_evidence_found"],
        "tetramer_context_annotation_present": cert["tetramer_context_annotation_present"],
        "tetramer_claim_allowed": cert["tetramer_claim_allowed"],
        "external_couplings_present": cert["external_couplings_present"],
        "next_decision": str(paths["next_decision"]),
        "certificate": str(paths["certificate"]),
        "manifest": str(paths["manifest"]),
        "assembly_preflight": str(paths["assembly_preflight"]),
        "pore_filter_annotation": str(paths["pore_filter_annotation"]),
        "membrane_topology_annotation": str(paths["membrane_topology_annotation"]),
        "external_coupling_availability": str(paths["external_coupling_availability"]),
        "report": str(paths["report"]),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
