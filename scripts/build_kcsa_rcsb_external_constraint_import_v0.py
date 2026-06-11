#!/usr/bin/env python3
from __future__ import annotations

"""Build real external KcsA source-import files from RCSB PDB 1BL8.

This is a V32 source-import/acquisition helper, not a folding solver.
It downloads the external RCSB PDB structure 1BL8, extracts two concrete
constraint files, and rewrites the V32 import manifest so the preflight can
select KcsA for the next claim-disabled readout.

Scientific boundary:
- Source is external coordinate-derived evidence from RCSB PDB 1BL8.
- claim_allowed stays false.
- new_MD_allowed stays false.
- Do not use this as a de novo native-contact benchmark for KcsA.
"""

import argparse
import csv
import json
import math
import os
import sys
import urllib.request
from datetime import date
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(os.environ.get("REPO_ROOT", Path(__file__).resolve().parents[1])).resolve()
PDB_ID = "1BL8"
RCSB_PAGE = "https://www.rcsb.org/structure/1BL8"
PDB_URL = "https://files.rcsb.org/download/1BL8.pdb"
SOURCE_NAME = "RCSB PDB 1BL8 / Doyle et al. Science 1998"
SOURCE_CITATION = "RCSB PDB 1BL8; Doyle DA et al. Science 280:69-77 (1998), DOI:10.1126/science.280.5360.69"
SOURCE_DATE = f"accessed_{date.today().isoformat()}"

RAW_DIR = REPO_ROOT / "data" / "external_constraints" / "KcsA" / "raw_sources"
PORE_DIR = REPO_ROOT / "data" / "external_constraints" / "KcsA" / "pore_filter"
IFACE_DIR = REPO_ROOT / "data" / "external_constraints" / "KcsA" / "assembly_interface"
OUT_DIR = REPO_ROOT / "first_contact_clean_pharmacotopology_layer_run" / "KCSA_RCSB_EXTERNAL_SOURCE_IMPORT_BUILD"
MANIFEST = REPO_ROOT / "data" / "external_constraints" / "v32_external_constraint_source_import_manifest.json"
PDB_PATH = RAW_DIR / "1BL8.pdb"
PORE_CSV = PORE_DIR / "kcsa_1bl8_pore_filter_external_contacts.csv"
IFACE_CSV = IFACE_DIR / "kcsa_1bl8_assembly_interface_external_contacts.csv"
CERT_PATH = OUT_DIR / "kcsa_rcsb_external_source_import_build_certificate.json"

AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
WATERS = {"HOH", "WAT", "DOD"}
ION_NAMES = {"K", "K+", "POT"}
OXYGEN_PREFIXES = ("O",)


def rel(path: Path) -> str:
    return str(path.resolve().relative_to(REPO_ROOT)).replace("\\", "/")


def download_pdb(force: bool = False) -> None:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    if PDB_PATH.exists() and not force:
        return
    print(f"Downloading external RCSB source: {PDB_URL}")
    try:
        with urllib.request.urlopen(PDB_URL, timeout=45) as response:
            data = response.read()
    except Exception as exc:  # noqa: BLE001 - user-facing CLI helper
        raise SystemExit(
            "Could not download RCSB 1BL8. Check internet, then rerun. "
            f"URL: {PDB_URL}\nOriginal error: {exc}"
        ) from exc
    if not data.startswith((b"HEADER", b"TITLE", b"REMARK", b"ATOM", b"HETATM")):
        raise SystemExit("Downloaded file does not look like a PDB file; refusing to build source import.")
    PDB_PATH.write_bytes(data)


def parse_pdb(path: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        rec = line[0:6].strip()
        if rec not in {"ATOM", "HETATM"}:
            continue
        try:
            serial = int(line[6:11])
            atom = line[12:16].strip()
            altloc = line[16:17].strip()
            resname = line[17:20].strip()
            chain = line[21:22].strip() or "?"
            resseq = int(line[22:26])
            icode = line[26:27].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = (line[76:78].strip() if len(line) >= 78 else "") or atom[0]
        except Exception:
            continue
        if altloc not in {"", "A"}:
            continue
        atoms.append(
            {
                "record": rec,
                "serial": serial,
                "atom": atom,
                "resname": resname,
                "chain": chain,
                "resseq": resseq,
                "icode": icode,
                "x": x,
                "y": y,
                "z": z,
                "element": element.upper(),
            }
        )
    if not atoms:
        raise SystemExit(f"No atoms parsed from {path}")
    return atoms


def is_protein_atom(atom: dict[str, object]) -> bool:
    return atom["record"] == "ATOM" and str(atom["resname"]) not in WATERS and str(atom["element"]).upper() != "H"


def is_potassium(atom: dict[str, object]) -> bool:
    return atom["record"] == "HETATM" and (str(atom["resname"]).strip().upper() in ION_NAMES or str(atom["element"]).upper() == "K")


def dist(a: dict[str, object], b: dict[str, object]) -> float:
    return math.sqrt(
        (float(a["x"]) - float(b["x"])) ** 2
        + (float(a["y"]) - float(b["y"])) ** 2
        + (float(a["z"]) - float(b["z"])) ** 2
    )


def residue_key(atom: dict[str, object]) -> tuple[str, int, str, str]:
    return (str(atom["chain"]), int(atom["resseq"]), str(atom["icode"]), str(atom["resname"]))


def group_residues(atoms: Iterable[dict[str, object]]) -> dict[tuple[str, int, str, str], list[dict[str, object]]]:
    residues: dict[tuple[str, int, str, str], list[dict[str, object]]] = {}
    for atom in atoms:
        if not is_protein_atom(atom):
            continue
        residues.setdefault(residue_key(atom), []).append(atom)
    return residues


def detect_tvgyg_filter(residues: dict[tuple[str, int, str, str], list[dict[str, object]]]) -> list[tuple[str, int, str, str]]:
    by_chain: dict[str, list[tuple[str, int, str, str]]] = {}
    for key in residues:
        chain, _resseq, _icode, _resname = key
        by_chain.setdefault(chain, []).append(key)
    found: list[tuple[str, int, str, str]] = []
    for chain, keys in by_chain.items():
        ordered = sorted(keys, key=lambda k: (k[1], k[2]))
        seq = "".join(AA1.get(k[3], "X") for k in ordered)
        start = seq.find("TVGYG")
        if start >= 0:
            found.extend(ordered[start : start + 5])
    if found:
        return found

    # Conservative fallback for canonical KcsA numbering used in the literature.
    fallback = [key for key in residues if 75 <= key[1] <= 79]
    if fallback:
        return sorted(fallback, key=lambda k: (k[0], k[1], k[2]))
    raise SystemExit("Could not detect TVGYG selectivity filter or canonical residues 75-79 in 1BL8.")


def min_atom_distance(a_atoms: list[dict[str, object]], b_atoms: list[dict[str, object]]) -> tuple[float, str, str]:
    best = 999.0
    best_a = ""
    best_b = ""
    for a in a_atoms:
        if str(a["element"]).upper() == "H":
            continue
        for b in b_atoms:
            if str(b["element"]).upper() == "H":
                continue
            d = dist(a, b)
            if d < best:
                best = d
                best_a = str(a["atom"])
                best_b = str(b["atom"])
    return best, best_a, best_b


def write_pore_filter_csv(atoms: list[dict[str, object]], residues: dict[tuple[str, int, str, str], list[dict[str, object]]], filter_keys: list[tuple[str, int, str, str]]) -> list[dict[str, object]]:
    PORE_DIR.mkdir(parents=True, exist_ok=True)
    ions = [a for a in atoms if is_potassium(a)]
    rows: list[dict[str, object]] = []

    # Primary evidence rows: filter oxygen atoms coordinating external K+ ions.
    for key in filter_keys:
        chain, resseq, icode, resname = key
        oxy_atoms = [a for a in residues[key] if str(a["atom"]).upper().startswith(OXYGEN_PREFIXES)]
        if not oxy_atoms:
            oxy_atoms = residues[key]
        for ion in ions:
            best = min((dist(a, ion), a) for a in oxy_atoms)
            d, atom = best
            if d <= 3.5:
                rows.append(
                    {
                        "constraint_id": f"KCSA_1BL8_PORE_{len(rows)+1:04d}",
                        "chain": chain,
                        "residue_number": resseq,
                        "icode": icode,
                        "residue_name": resname,
                        "filter_motif": "TVGYG",
                        "atom_name": atom["atom"],
                        "ion_serial": ion["serial"],
                        "ion_name": ion["resname"],
                        "ion_chain": ion["chain"],
                        "ion_residue_number": ion["resseq"],
                        "min_distance_angstrom": f"{d:.3f}",
                        "constraint_class": "pore_filter_potassium_coordination_contact",
                        "source_id": PDB_ID,
                        "source_url": RCSB_PAGE,
                        "notes": "External RCSB coordinate-derived pore/filter contact; source-import only, not de-novo folding claim.",
                    }
                )

    # Fallback evidence rows: local TVGYG filter residue geometry if K+ contacts are not parsed.
    # This keeps the file source-derived and non-empty, while the certificate records row counts.
    if not rows:
        ordered = sorted(filter_keys, key=lambda k: (k[0], k[1], k[2]))
        for chain in sorted({k[0] for k in ordered}):
            chain_keys = [k for k in ordered if k[0] == chain]
            for a_key, b_key in zip(chain_keys, chain_keys[1:]):
                d, atom_a, atom_b = min_atom_distance(residues[a_key], residues[b_key])
                rows.append(
                    {
                        "constraint_id": f"KCSA_1BL8_PORE_{len(rows)+1:04d}",
                        "chain": chain,
                        "residue_i": a_key[1],
                        "residue_name_i": a_key[3],
                        "atom_i": atom_a,
                        "residue_j": b_key[1],
                        "residue_name_j": b_key[3],
                        "atom_j": atom_b,
                        "filter_motif": "TVGYG",
                        "min_distance_angstrom": f"{d:.3f}",
                        "constraint_class": "pore_filter_residue_geometry_contact",
                        "source_id": PDB_ID,
                        "source_url": RCSB_PAGE,
                        "notes": "External RCSB coordinate-derived TVGYG pore/filter geometry; source-import only, not de-novo folding claim.",
                    }
                )

    fieldnames = sorted({key for row in rows for key in row.keys()})
    with PORE_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_interface_csv(residues: dict[tuple[str, int, str, str], list[dict[str, object]]], cutoff: float, max_rows: int) -> list[dict[str, object]]:
    IFACE_DIR.mkdir(parents=True, exist_ok=True)
    keys = sorted(residues.keys(), key=lambda k: (k[0], k[1], k[2]))
    rows_raw: list[tuple[float, tuple[str, int, str, str], tuple[str, int, str, str], str, str]] = []
    for idx, a_key in enumerate(keys):
        for b_key in keys[idx + 1 :]:
            if a_key[0] == b_key[0]:
                continue
            d, atom_a, atom_b = min_atom_distance(residues[a_key], residues[b_key])
            if d <= cutoff:
                rows_raw.append((d, a_key, b_key, atom_a, atom_b))
    rows_raw.sort(key=lambda row: row[0])
    rows: list[dict[str, object]] = []
    for d, a_key, b_key, atom_a, atom_b in rows_raw[:max_rows]:
        rows.append(
            {
                "constraint_id": f"KCSA_1BL8_IFACE_{len(rows)+1:04d}",
                "chain_a": a_key[0],
                "residue_i": a_key[1],
                "icode_i": a_key[2],
                "residue_name_i": a_key[3],
                "atom_i": atom_a,
                "chain_b": b_key[0],
                "residue_j": b_key[1],
                "icode_j": b_key[2],
                "residue_name_j": b_key[3],
                "atom_j": atom_b,
                "min_distance_angstrom": f"{d:.3f}",
                "constraint_class": "assembly_interface_heavy_atom_contact",
                "source_id": PDB_ID,
                "source_url": RCSB_PAGE,
                "notes": "External RCSB coordinate-derived inter-chain/tetramer interface contact; source-import only, not de-novo folding claim.",
            }
        )
    if not rows:
        raise SystemExit("No inter-chain interface contacts found; refusing to write empty assembly/interface source file.")
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with IFACE_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def write_manifest() -> None:
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "kind": "V32_EXTERNAL_CONSTRAINT_SOURCE_IMPORT_MANIFEST_v0",
        "claim_allowed": False,
        "instructions": [
            "Generated by scripts/build_kcsa_rcsb_external_constraint_import_v0.py from external RCSB PDB 1BL8.",
            "This enables KcsA source-import/readout gating only; it is not a de novo folding benchmark.",
            "Do not treat coordinate-derived constraints as proof of universal folding prediction.",
        ],
        "rows": [
            {
                "target": "KcsA",
                "file_path": rel(PORE_CSV),
                "evidence_type": "pore_filter_coupling",
                "state_or_context": "TVGYG_selectivity_filter_or_pore_filter_external_coordinate_contact",
                "source_name": SOURCE_NAME,
                "source_url_or_citation": SOURCE_CITATION,
                "source_date_or_version": SOURCE_DATE,
                "provenance_notes": "Real external RCSB PDB coordinate-derived pore/filter contact file. Claim-disabled source import only; not internal runtime report.",
            },
            {
                "target": "KcsA",
                "file_path": rel(IFACE_CSV),
                "evidence_type": "assembly_interface_constraint",
                "state_or_context": "tetramer_chain_interface_context_external_coordinate_contact",
                "source_name": SOURCE_NAME,
                "source_url_or_citation": SOURCE_CITATION,
                "source_date_or_version": SOURCE_DATE,
                "provenance_notes": "Real external RCSB PDB coordinate-derived assembly/interface contact file. Claim-disabled source import only; not internal runtime report.",
            },
        ],
    }
    MANIFEST.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def write_certificate(pore_rows: list[dict[str, object]], iface_rows: list[dict[str, object]], filter_keys: list[tuple[str, int, str, str]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cert = {
        "kind": "KCSA_RCSB_EXTERNAL_SOURCE_IMPORT_BUILD_v0",
        "claim_allowed": False,
        "new_MD_allowed": False,
        "new_md_executed": False,
        "source_boundary": "external_RCSB_coordinate_derived_constraints_for_source_import_only_not_de_novo_folding_claim",
        "source_name": SOURCE_NAME,
        "source_url": RCSB_PAGE,
        "source_citation": SOURCE_CITATION,
        "source_date_or_version": SOURCE_DATE,
        "artifacts": {
            "raw_pdb": rel(PDB_PATH),
            "pore_filter_csv": rel(PORE_CSV),
            "assembly_interface_csv": rel(IFACE_CSV),
            "manifest": rel(MANIFEST),
        },
        "counts": {
            "detected_filter_residue_count": len(filter_keys),
            "pore_filter_contact_rows": len(pore_rows),
            "assembly_interface_contact_rows": len(iface_rows),
        },
        "detected_filter_residues": [
            {"chain": k[0], "residue_number": k[1], "icode": k[2], "residue_name": k[3]} for k in filter_keys
        ],
        "next_commands": [
            "python3 scripts/run_v32_external_constraint_source_import_preflight_v0.py",
            "python3 scripts/print_v32_external_constraint_source_import_preflight.py",
        ],
    }
    CERT_PATH.write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build KcsA RCSB 1BL8 external source-import files for V32.")
    parser.add_argument("--force-download", action="store_true", help="redownload 1BL8.pdb even if a local copy exists")
    parser.add_argument("--interface-cutoff", type=float, default=4.5, help="heavy-atom inter-chain contact cutoff in Å")
    parser.add_argument("--max-interface-rows", type=int, default=250, help="maximum assembly/interface contact rows to write")
    args = parser.parse_args()

    download_pdb(force=args.force_download)
    atoms = parse_pdb(PDB_PATH)
    residues = group_residues(atoms)
    filter_keys = detect_tvgyg_filter(residues)
    pore_rows = write_pore_filter_csv(atoms, residues, filter_keys)
    iface_rows = write_interface_csv(residues, cutoff=args.interface_cutoff, max_rows=args.max_interface_rows)
    write_manifest()
    write_certificate(pore_rows, iface_rows, filter_keys)

    print(json.dumps({
        "ok": True,
        "kind": "KCSA_RCSB_EXTERNAL_SOURCE_IMPORT_BUILD_v0",
        "claim_allowed": False,
        "new_MD_allowed": False,
        "raw_pdb": rel(PDB_PATH),
        "pore_filter_csv": rel(PORE_CSV),
        "pore_filter_rows": len(pore_rows),
        "assembly_interface_csv": rel(IFACE_CSV),
        "assembly_interface_rows": len(iface_rows),
        "manifest": rel(MANIFEST),
        "certificate": rel(CERT_PATH),
        "next": [
            "python3 scripts/run_v32_external_constraint_source_import_preflight_v0.py",
            "python3 scripts/print_v32_external_constraint_source_import_preflight.py",
        ],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
