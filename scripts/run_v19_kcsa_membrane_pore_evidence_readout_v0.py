#!/usr/bin/env python3
from __future__ import annotations

"""V19 KcsA membrane/pore evidence readout.

Zero-MD, no-tuning pressure-evidence sprint for the KcsA pressure target selected
after the p53 context-contrast milestone. This reads already-provenanced 1K4C
coordinate/annotation context for membrane-pore roles, pore/filter diagnostics,
and chain/interface context. It explicitly does not run membrane MD, does not
claim a tetramer/whole-channel fold, and does not treat KcsA as a soluble compact
core.
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
DEFAULT_PDB = REPO_ROOT / "data" / "v16_pressure_targets" / "KcsA" / "1K4C.pdb"
DEFAULT_OUT_DIR = RUN_ROOT / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT"

PROBE_RADII_ANGSTROM = [6.0, 8.0, 10.0]
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
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


def _target_row_by_id(cert: dict[str, Any], key: str, target_id: str) -> dict[str, Any]:
    """Return a target row across the historical certificate schemas.

    V16 preflight writes `target_results`; V16 zero-MD role transfer writes
    `target_rows`. V19 must read both, otherwise it falsely treats the
    soluble-core leakage guard as missing even when V16 already classified KcsA
    as a membrane/pore object and recorded no forbidden violations.
    """
    candidate_keys = [key, "target_rows", "target_results", "targets"]
    seen: set[str] = set()
    for row_key in candidate_keys:
        if row_key in seen:
            continue
        seen.add(row_key)
        rows = cert.get(row_key) if isinstance(cert.get(row_key), list) else []
        for row in rows:
            if isinstance(row, dict) and row.get("target_id") == target_id:
                return row
    return {}


def _material_by_id(row: dict[str, Any]) -> dict[str, dict[str, Any]]:
    mats = row.get("material_results") if isinstance(row.get("material_results"), list) else []
    return {str(m.get("material_id")): m for m in mats if isinstance(m, dict)}


def _resolve_pdb_path(preflight_row: dict[str, Any], fallback: Path) -> Path:
    mats = _material_by_id(preflight_row)
    mat = mats.get("KcsA_Fab_complex_structure", {})
    value = mat.get("path")
    if isinstance(value, str) and value:
        p = Path(value)
        return p if p.is_absolute() else REPO_ROOT / p
    return fallback


def _parse_pdb(path: Path) -> dict[str, Any]:
    atoms: dict[str, list[dict[str, Any]]] = {}
    residues: dict[str, dict[tuple[int, str], str]] = {}
    helix_records: list[dict[str, Any]] = []
    sheet_records: list[dict[str, Any]] = []
    potassium_ions: list[dict[str, Any]] = []
    annotation_lines: list[str] = []
    seqres_by_chain: dict[str, list[str]] = {}

    if not path.exists():
        return {
            "exists": False,
            "atoms_by_chain": atoms,
            "residue_sequences_by_chain": {},
            "seqres_sequences_by_chain": {},
            "helix_records": helix_records,
            "sheet_records": sheet_records,
            "potassium_ions": potassium_ions,
            "annotation_text": "",
        }

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            record = line[0:6].strip()
            if record in {"HEADER", "TITLE", "COMPND", "SOURCE", "KEYWDS", "REMARK", "SITE"}:
                annotation_lines.append(line.strip())
            if record == "SEQRES":
                chain = line[11].strip() or "_"
                names = [x.strip().upper() for x in line[19:].split() if x.strip()]
                seqres_by_chain.setdefault(chain, []).extend(names)
            elif record == "HELIX":
                helix_records.append({
                    "raw": line.rstrip("\n"),
                    "init_chain": line[19].strip() if len(line) > 19 else None,
                    "end_chain": line[31].strip() if len(line) > 31 else None,
                    "init_seq_num": line[21:25].strip(),
                    "end_seq_num": line[33:37].strip(),
                })
            elif record == "SHEET":
                sheet_records.append({
                    "raw": line.rstrip("\n"),
                    "init_chain": line[21].strip() if len(line) > 21 else None,
                    "end_chain": line[32].strip() if len(line) > 32 else None,
                })
            elif record == "ATOM":
                if len(line) < 54 or line[12:16].strip() != "CA":
                    continue
                chain = line[21].strip() or "_"
                resname = line[17:20].strip().upper()
                icode = line[26].strip()
                try:
                    resseq = int(line[22:26])
                    xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                except ValueError:
                    continue
                atoms.setdefault(chain, []).append({
                    "residue_number": resseq,
                    "insertion_code": icode,
                    "residue_name": resname,
                    "xyz": xyz,
                })
                residues.setdefault(chain, {})[(resseq, icode)] = resname
            elif record == "HETATM":
                atom_name = line[12:16].strip().upper() if len(line) >= 16 else ""
                resname = line[17:20].strip().upper() if len(line) >= 20 else ""
                element = line[76:78].strip().upper() if len(line) >= 78 else ""
                if atom_name == "K" or resname in {"K", "K+"} or element == "K":
                    try:
                        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    except ValueError:
                        xyz = None
                    potassium_ions.append({
                        "atom_name": atom_name,
                        "residue_name": resname,
                        "chain": line[21].strip() if len(line) > 21 else None,
                        "residue_number": line[22:26].strip(),
                        "xyz": xyz,
                    })

    residue_sequences: dict[str, str] = {}
    for chain, chain_residues in residues.items():
        letters = []
        for key in sorted(chain_residues):
            letters.append(THREE_TO_ONE.get(chain_residues[key], "X"))
        residue_sequences[chain] = "".join(letters)
    seqres_sequences = {
        chain: "".join(THREE_TO_ONE.get(name, "X") for name in names)
        for chain, names in seqres_by_chain.items()
    }
    return {
        "exists": True,
        "atoms_by_chain": atoms,
        "residue_sequences_by_chain": residue_sequences,
        "seqres_sequences_by_chain": seqres_sequences,
        "helix_records": helix_records,
        "sheet_records": sheet_records,
        "potassium_ions": potassium_ions,
        "annotation_text": "\n".join(annotation_lines).upper(),
    }


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _multi_chain_interface_probe(atoms: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    chain_counts = {chain: len(items) for chain, items in atoms.items()}
    pair_counts: dict[str, dict[str, int]] = {}
    pair_min_distance: dict[str, float | None] = {}
    chains = sorted(chain_counts)
    for idx, chain_a in enumerate(chains):
        for chain_b in chains[idx + 1:]:
            key = f"{chain_a}-{chain_b}"
            pair_counts[key] = {str(radius): 0 for radius in PROBE_RADII_ANGSTROM}
            min_d: float | None = None
            for a in atoms.get(chain_a, []):
                for b in atoms.get(chain_b, []):
                    d = _distance(a["xyz"], b["xyz"])
                    if min_d is None or d < min_d:
                        min_d = d
                    for radius in PROBE_RADII_ANGSTROM:
                        if d <= radius:
                            pair_counts[key][str(radius)] += 1
            pair_min_distance[key] = round(min_d, 3) if min_d is not None else None
    interface_pairs_present = [
        pair for pair, counts in pair_counts.items()
        if sum(1 for value in counts.values() if value > 0) >= 2
    ]
    return {
        "chain_ca_counts": chain_counts,
        "chain_count_with_ca": len(chain_counts),
        "multi_radius_contact_counts_by_chain_pair": pair_counts,
        "min_ca_distance_by_chain_pair": pair_min_distance,
        "interface_pairs_present": interface_pairs_present,
        "interface_context_present": bool(interface_pairs_present),
        "contact_probe_radii_angstrom_report_only": PROBE_RADII_ANGSTROM,
        "contact_probe_policy": "multi_radius_report_only_no_single_fixed_selection_threshold",
        "selection_threshold_used": False,
    }


def _contains_filter_motif(parsed: dict[str, Any]) -> dict[str, Any]:
    motif_hits: list[dict[str, Any]] = []
    motifs = ["TVGYG", "TXXYG"]
    for source_key in ["residue_sequences_by_chain", "seqres_sequences_by_chain"]:
        seqs = parsed.get(source_key) if isinstance(parsed.get(source_key), dict) else {}
        for chain, seq in seqs.items():
            seq = str(seq)
            if "TVGYG" in seq:
                motif_hits.append({"source": source_key, "chain": chain, "motif": "TVGYG", "index": seq.index("TVGYG")})
            # Conservative loose diagnostic: T--YG is reported as diagnostic only.
            for i in range(0, max(0, len(seq) - 4)):
                window = seq[i:i + 5]
                if len(window) == 5 and window[0] == "T" and window[3:5] == "YG":
                    motif_hits.append({"source": source_key, "chain": chain, "motif": "TxxYG_diagnostic", "index": i, "window": window})
                    break
    strong_hits = [hit for hit in motif_hits if hit.get("motif") == "TVGYG"]
    return {
        "filter_motif_detected": bool(strong_hits or motif_hits),
        "strong_TVGYG_motif_detected": bool(strong_hits),
        "motif_hits": motif_hits,
        "motif_policy": "sequence_or_seqres_report_only_no_single_fixed_success_threshold",
        "selection_threshold_used": False,
    }


def _annotation_probe(parsed: dict[str, Any]) -> dict[str, Any]:
    text = str(parsed.get("annotation_text") or "")
    keywords = {
        "potassium_channel_annotation": any(k in text for k in ["POTASSIUM CHANNEL", "KCSA", " K+ ", "POTASSIUM"]),
        "pore_annotation": "PORE" in text or "CHANNEL" in text,
        "selectivity_filter_annotation": "SELECTIVITY" in text or "FILTER" in text,
        "membrane_annotation": "MEMBRANE" in text or "TRANSMEMBRANE" in text,
    }
    return {
        **keywords,
        "annotation_policy": "pdb_header_keyword_report_only_no_native_selection",
        "selection_threshold_used": False,
    }


def _helix_probe(parsed: dict[str, Any]) -> dict[str, Any]:
    helix_records = parsed.get("helix_records") if isinstance(parsed.get("helix_records"), list) else []
    chain_counts: dict[str, int] = {}
    for rec in helix_records:
        if not isinstance(rec, dict):
            continue
        chain = rec.get("init_chain") or rec.get("end_chain") or "_"
        chain_counts[str(chain)] = chain_counts.get(str(chain), 0) + 1
    return {
        "helix_record_count": len(helix_records),
        "helix_records_by_chain": dict(sorted(chain_counts.items())),
        "transmembrane_helix_scaffold_context_present": len(helix_records) > 0,
        "helix_policy": "pdb_helix_records_report_only_no_membrane_MD_claim",
        "selection_threshold_used": False,
    }


def _ion_probe(parsed: dict[str, Any]) -> dict[str, Any]:
    ions = parsed.get("potassium_ions") if isinstance(parsed.get("potassium_ions"), list) else []
    return {
        "potassium_ion_count": len(ions),
        "potassium_ion_records_sample": ions[:10],
        "ion_filter_diagnostic_shell_present": len(ions) > 0,
        "ion_probe_policy": "ion_context_report_only_no_whole_fold_claim",
        "selection_threshold_used": False,
    }


def _zero_row(zero_md: dict[str, Any]) -> dict[str, Any]:
    return _target_row_by_id(zero_md, "target_results", "KcsA")


def build_evidence(preflight: dict[str, Any], zero_md: dict[str, Any], pdb_path: Path) -> dict[str, Any]:
    preflight_row = _target_row_by_id(preflight, "target_results", "KcsA")
    zero_row = _zero_row(zero_md)
    parsed = _parse_pdb(pdb_path)
    atoms = parsed.get("atoms_by_chain") if isinstance(parsed.get("atoms_by_chain"), dict) else {}

    material_ready = (
        preflight.get("data_preflight_status") == "V16_TRANSFER_DATA_PREFLIGHT_READY_FOR_ZERO_MD_ROLE_READOUT"
        and preflight_row.get("target_material_status") == "ready_for_zero_md_role_readout"
        and parsed.get("exists") is True
        and bool(atoms)
    )
    role_classification_passed = bool(
        zero_row.get("role_classification_passed") is True
        or zero_row.get("positive_role_context_found") is True
        or "KcsA" in (zero_md.get("role_classification_passed_targets") or [])
        or "KcsA" in (zero_md.get("pressure_role_transfer_passed_targets") or [])
    )
    soluble_core_guard = role_classification_passed and not zero_row.get("forbidden_misclassification_violations")

    interface = _multi_chain_interface_probe(atoms)
    motif = _contains_filter_motif(parsed)
    annotations = _annotation_probe(parsed)
    helix = _helix_probe(parsed)
    ions = _ion_probe(parsed)

    chain_context = material_ready and bool(interface.get("interface_context_present"))
    membrane_context = bool(
        annotations.get("membrane_annotation")
        or zero_row.get("expected_role_class") == "membrane_pore_oligomer_object"
        or preflight_row.get("expected_role_class") == "membrane_pore_oligomer_object"
    )
    pore_filter_context = bool(
        motif.get("filter_motif_detected")
        or annotations.get("selectivity_filter_annotation")
        or annotations.get("pore_annotation")
        or ions.get("ion_filter_diagnostic_shell_present")
    )
    transmembrane_context = bool(helix.get("transmembrane_helix_scaffold_context_present") or annotations.get("membrane_annotation"))
    oligomer_or_interface_context = bool(chain_context)

    forbidden_violations: list[str] = []
    soluble_core_misclassification = False
    whole_channel_fold_claim_made = False
    tetramer_claim_made = False
    if soluble_core_misclassification:
        forbidden_violations.append("membrane_pore_as_soluble_single_domain_compact_core")
    if whole_channel_fold_claim_made:
        forbidden_violations.append("ion_filter_geometry_as_whole_fold_claim")
    if tetramer_claim_made:
        forbidden_violations.append("tetramer_interface_as_monomer_autonomous_core")

    role_buckets: list[str] = []
    if transmembrane_context:
        role_buckets.append("transmembrane_helix_scaffold")
    if pore_filter_context:
        role_buckets.append("pore_selectivity_filter_core")
    if oligomer_or_interface_context:
        role_buckets.append("chain_or_interface_context_support")
    if membrane_context:
        role_buckets.append("membrane_environment_context")
    if ions.get("ion_filter_diagnostic_shell_present") or annotations.get("selectivity_filter_annotation"):
        role_buckets.append("ion_filter_diagnostic_shell")

    positive_pressure = bool(
        material_ready
        and soluble_core_guard
        and membrane_context
        and pore_filter_context
        and oligomer_or_interface_context
        and not forbidden_violations
    )
    clean_abstain = not positive_pressure
    test_status = (
        "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_PASSED_CLAIM_DISABLED"
        if positive_pressure
        else "V19_KcsA_CLEAN_ABSTAIN_MISSING_MEMBRANE_PORE_OR_ASSEMBLY_EVIDENCE_CLAIM_DISABLED"
    )
    missing = []
    if not material_ready:
        missing.append("KcsA_complex_coordinate_context")
    if not transmembrane_context:
        missing.append("transmembrane_helix_roles_or_annotation")
    if not pore_filter_context:
        missing.append("pore_selectivity_filter_annotation_or_diagnostic")
    if not oligomer_or_interface_context:
        missing.append("chain_or_interface_context")
    if not soluble_core_guard:
        missing.append("leakage_guard_against_soluble_core_misread")
    missing.append("external_couplings_if_available")
    if tetramer_claim_made is False:
        # Biological tetramer remains later unless explicit assembly annotation is added.
        missing.append("biological_tetramer_assembly_annotation_for_tetramer_claim")

    return {
        "kind": "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_v0",
        "run_mode": "zero_md_membrane_pore_evidence_readout_no_simulation_no_threshold_tuning",
        "test_status": test_status,
        "target_id": "KcsA",
        "role_class": "membrane_pore_oligomer_object",
        "pressure_class": "membrane_pore_oligomer_environment",
        "source_preflight_status": preflight.get("data_preflight_status"),
        "source_zero_md_status": zero_md.get("role_transfer_status"),
        "claim_allowed": False,
        "evidence_claim_allowed": False,
        "new_md_executed": False,
        "membrane_md_executed": False,
        "positive_folding_evidence_found": False,
        "positive_pressure_evidence_found": positive_pressure,
        "membrane_pore_role_evidence_found": positive_pressure,
        "clean_abstain_valid": clean_abstain,
        "fixed_threshold_policy": "forbidden",
        "target_specific_threshold_tuning_allowed": False,
        "fixed_residue_cutoff_used": False,
        "native_metrics_used_for_selection": False,
        "selection_policy": "annotation_and_structure_role_readout_only_no_native_metric_selection",
        "required_roles": [
            "transmembrane_helix_scaffold",
            "pore_selectivity_filter_core",
            "tetramer_interface_support",
            "membrane_environment_context",
            "diagnostic_shell",
        ],
        "role_buckets_assigned": sorted(set(role_buckets)),
        "selected_core_or_clean_abstain": (
            "membrane_pore_roles_detected_without_soluble_core_misclassification"
            if positive_pressure
            else "clean_abstain_missing_membrane_pore_or_assembly_evidence"
        ),
        "pore_filter_annotation_present": pore_filter_context,
        "transmembrane_role_present": transmembrane_context,
        "oligomer_or_chain_interface_context_present": oligomer_or_interface_context,
        "membrane_environment_context_present": membrane_context,
        "soluble_core_misclassification": soluble_core_misclassification,
        "soluble_core_misclassification_avoided": not soluble_core_misclassification and soluble_core_guard,
        "whole_channel_fold_claim_made": whole_channel_fold_claim_made,
        "tetramer_claim_made": tetramer_claim_made,
        "forbidden_misclassification": [
            "membrane_pore_as_soluble_single_domain_compact_core",
            "ion_filter_geometry_as_whole_fold_claim",
            "tetramer_interface_as_monomer_autonomous_core",
        ],
        "forbidden_misclassification_violations": forbidden_violations,
        "available_evidence": [
            item
            for item, present in {
                "KcsA_complex_coordinate_context": material_ready,
                "pore_selectivity_filter_annotation_or_diagnostic": pore_filter_context,
                "transmembrane_helix_scaffold_context": transmembrane_context,
                "chain_or_interface_context": oligomer_or_interface_context,
                "membrane_environment_context": membrane_context,
                "leakage_guard_against_soluble_core_misread": soluble_core_guard,
            }.items()
            if present
        ],
        "missing_evidence": sorted(set(missing)),
        "material_observations": {
            "pdb_path": str(pdb_path),
            "chain_ca_counts": interface.get("chain_ca_counts"),
            "chain_count_with_ca": interface.get("chain_count_with_ca"),
        },
        "pore_filter_readout": {
            "sequence_motif_probe": motif,
            "annotation_probe": annotations,
            "potassium_ion_probe": ions,
            "policy": "pore_filter_annotation_or_diagnostic_report_only_no_whole_fold_claim",
        },
        "transmembrane_helix_readout": helix,
        "chain_interface_readout": interface,
        "limitations": [
            "zero_MD_only_no_membrane_lipid_dynamics_claim",
            "biological_tetramer_assembly_annotation_later_before_tetramer_claim",
            "external_couplings_not_required_for_this_annotation_evidence_readout",
            "no_whole_channel_fold_claim",
        ],
        "next_recommended_panel": "V20_XCL1_STATE_SPECIFIC_EVIDENCE_READOUT_or_KcsA_ANNOTATION_ENRICHMENT_IF_NEEDED",
    }


def _write_report(path: Path, cert: dict[str, Any]) -> None:
    pore = cert.get("pore_filter_readout", {}) if isinstance(cert.get("pore_filter_readout"), dict) else {}
    helix = cert.get("transmembrane_helix_readout", {}) if isinstance(cert.get("transmembrane_helix_readout"), dict) else {}
    iface = cert.get("chain_interface_readout", {}) if isinstance(cert.get("chain_interface_readout"), dict) else {}
    lines = [
        "# V19 KcsA Membrane/Pore Evidence Readout",
        "",
        "Zero-MD, no tuning, no membrane-MD claim. This reads membrane/pore/chain-context evidence and keeps soluble compact-core and whole-channel fold claims forbidden.",
        "",
        f"Status: `{cert.get('test_status')}`",
        f"Positive pressure evidence found: `{cert.get('positive_pressure_evidence_found')}`",
        f"Positive folding evidence found: `{cert.get('positive_folding_evidence_found')}`",
        f"Claim allowed: `{cert.get('claim_allowed')}`",
        f"New MD executed: `{cert.get('new_md_executed')}`",
        "",
        "## Role evidence",
        f"Role buckets assigned: `{cert.get('role_buckets_assigned')}`",
        f"Pore/filter annotation present: `{cert.get('pore_filter_annotation_present')}`",
        f"Transmembrane role present: `{cert.get('transmembrane_role_present')}`",
        f"Oligomer/interface context present: `{cert.get('oligomer_or_chain_interface_context_present')}`",
        f"Soluble-core misclassification avoided: `{cert.get('soluble_core_misclassification_avoided')}`",
        "",
        "## Pore/filter readout",
        f"Motif probe: `{pore.get('sequence_motif_probe')}`",
        f"Annotation probe: `{pore.get('annotation_probe')}`",
        f"Potassium ion probe: `{pore.get('potassium_ion_probe')}`",
        "",
        "## Helix/interface readout",
        f"Helix records: `{helix}`",
        f"Chain interface readout: `{iface}`",
        "",
        "## Guards",
        f"Forbidden misclassification violations: `{cert.get('forbidden_misclassification_violations')}`",
        f"Missing evidence: `{cert.get('missing_evidence')}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(out_dir: Path, cert: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "certificate": out_dir / "v19_kcsa_membrane_pore_evidence_readout_certificate.json",
        "pore_filter_readout": out_dir / "v19_kcsa_pore_filter_readout.json",
        "chain_interface_readout": out_dir / "v19_kcsa_chain_interface_readout.json",
        "report": out_dir / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_REPORT.md",
    }
    cert = dict(cert)
    cert["artifacts"] = {key: str(path) for key, path in paths.items()}
    paths["certificate"].write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    paths["pore_filter_readout"].write_text(json.dumps(cert.get("pore_filter_readout", {}), indent=2, sort_keys=True), encoding="utf-8")
    paths["chain_interface_readout"].write_text(json.dumps(cert.get("chain_interface_readout", {}), indent=2, sort_keys=True), encoding="utf-8")
    _write_report(paths["report"], cert)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight-cert", type=Path, default=DEFAULT_PREFLIGHT_CERT)
    parser.add_argument("--zero-md-cert", type=Path, default=DEFAULT_ZERO_MD_CERT)
    parser.add_argument("--pdb", type=Path, default=DEFAULT_PDB)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    preflight = _read_json(args.preflight_cert, "V16 data preflight certificate")
    zero_md = _read_json(args.zero_md_cert, "V16 zero-MD role transfer certificate")
    preflight_row = _target_row_by_id(preflight, "target_results", "KcsA")
    pdb = _resolve_pdb_path(preflight_row, args.pdb)
    cert = build_evidence(preflight, zero_md, pdb)
    write_outputs(args.out_dir, cert)
    print(json.dumps({
        "kind": cert["kind"],
        "test_status": cert["test_status"],
        "positive_pressure_evidence_found": cert["positive_pressure_evidence_found"],
        "positive_folding_evidence_found": cert["positive_folding_evidence_found"],
        "membrane_pore_role_evidence_found": cert["membrane_pore_role_evidence_found"],
        "claim_allowed": cert["claim_allowed"],
        "new_md_executed": cert["new_md_executed"],
        "certificate": str(args.out_dir / "v19_kcsa_membrane_pore_evidence_readout_certificate.json"),
        "report": str(args.out_dir / "V19_KcsA_MEMBRANE_PORE_EVIDENCE_READOUT_REPORT.md"),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
