from __future__ import annotations

"""Chai-1/Boltz-style conformational ensemble contact extraction for 4AKE.

This module does not run a neural model by itself.  It consumes a directory of
PDB/mmCIF structures produced by a single-sequence conformational sampler, or by
a bounded caller script, and converts the ensemble into auditable contact maps.

The native 4AKE structure is never read here.  Callers may attach native metrics
after the selected contact map is frozen.
"""

from dataclasses import asdict, dataclass
from math import sqrt
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from pharmacotopology.folding_native_contact_eval import ContactPair, contact_map_hash, normalized_contact_pairs
from pharmacotopology.folding_real_coordinate_visual_benchmark import CONTACT_CUTOFF_ANGSTROM, MIN_SEQUENCE_SEPARATION


CONFORMATIONAL_ENSEMBLE_KIND = "msa_free_conformational_ensemble_contact_source_v0"
OPEN_ADK_PROXY_RULE = (
    "select_conformer_with_shortest_D118_K136_CA_distance_when_residues_exist;"
    "fallback_to_largest_contact_count;native_truth_after_audit_only"
)


@dataclass(frozen=True)
class ParsedConformer:
    path: str
    format: str
    chain_id: str
    ca_count: int
    contact_count: int
    long_range_contact_count: int
    d118_k136_ca_distance_angstrom: float | None
    contact_map_hash: str
    selected_by_open_state_proxy: bool
    contacts: tuple[ContactPair, ...]

    def safe_dict(self) -> dict[str, object]:
        data = asdict(self)
        data.pop("contacts", None)
        return data


@dataclass(frozen=True)
class ConformationalEnsemblePacket:
    kind: str
    source_id: str
    structure_file_count: int
    parsed_conformer_count: int
    selected_conformer_path: str | None
    selected_contact_count: int
    selected_long_range_contact_count: int
    selected_contact_map_hash: str
    selection_rule: str
    contacts: tuple[ContactPair, ...]
    conformers: tuple[ParsedConformer, ...]
    alphafold_like_input_rejected: bool = False
    msa_used_by_this_module: bool = False
    native_truth_used_before_selection: bool = False
    coordinate_truth_used_before_selection: bool = False
    raw_sequence_exposed: bool = False

    def to_report_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "source_id": self.source_id,
            "structure_file_count": self.structure_file_count,
            "parsed_conformer_count": self.parsed_conformer_count,
            "selected_conformer_path": self.selected_conformer_path,
            "selected_contact_count": self.selected_contact_count,
            "selected_long_range_contact_count": self.selected_long_range_contact_count,
            "selected_contact_map_hash": self.selected_contact_map_hash,
            "selection_rule": self.selection_rule,
            "alphafold_like_input_rejected": self.alphafold_like_input_rejected,
            "msa_used_by_this_module": self.msa_used_by_this_module,
            "native_truth_used_before_selection": self.native_truth_used_before_selection,
            "coordinate_truth_used_before_selection": self.coordinate_truth_used_before_selection,
            "raw_sequence_exposed": self.raw_sequence_exposed,
            "conformers": [conformer.safe_dict() for conformer in self.conformers],
        }


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def _distance(left: Sequence[float], right: Sequence[float]) -> float:
    return sqrt(sum((float(a) - float(b)) ** 2 for a, b in zip(left, right)))


def _looks_alphafold_like(path: Path, source_id: str) -> bool:
    probe = f"{path.name} {source_id}".lower()
    return any(token in probe for token in ("alphafold", "af-", "af_", "colabfold"))


def _parse_pdb_ca(path: Path, *, chain_id: str | None) -> dict[int, tuple[float, float, float]]:
    coords: dict[int, tuple[float, float, float]] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom = line[12:16].strip()
            if atom != "CA":
                continue
            chain = line[21].strip() or "A"
            if chain_id and chain != chain_id:
                continue
            try:
                resseq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            coords.setdefault(resseq, (x, y, z))
    return coords


def _parse_cif_ca(path: Path, *, chain_id: str | None) -> dict[int, tuple[float, float, float]]:
    try:
        from Bio.PDB.MMCIFParser import MMCIFParser
    except Exception:
        return {}
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(path.stem, str(path))
    except Exception:
        return {}
    coords: dict[int, tuple[float, float, float]] = {}
    for model in structure:
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue
            for residue in chain:
                if "CA" not in residue:
                    continue
                hetflag, resseq, icode = residue.id
                if hetflag.strip():
                    continue
                atom = residue["CA"]
                xyz = atom.coord
                coords.setdefault(int(resseq), (float(xyz[0]), float(xyz[1]), float(xyz[2])))
        break
    return coords


def parse_ca_coordinates(path: Path, *, chain_id: str | None = "A") -> dict[int, tuple[float, float, float]]:
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return _parse_cif_ca(path, chain_id=chain_id)
    return _parse_pdb_ca(path, chain_id=chain_id)


def contacts_from_ca_coordinates(
    coords: Mapping[int, Sequence[float]],
    *,
    contact_cutoff_angstrom: float = CONTACT_CUTOFF_ANGSTROM,
    minimum_sequence_separation: int = MIN_SEQUENCE_SEPARATION,
) -> tuple[ContactPair, ...]:
    pairs: list[ContactPair] = []
    keys = sorted(int(key) for key in coords)
    for left_index, i in enumerate(keys):
        left = coords[i]
        for j in keys[left_index + 1 :]:
            if j - i < minimum_sequence_separation:
                continue
            if _distance(left, coords[j]) <= contact_cutoff_angstrom:
                pairs.append((i, j))
    return normalized_contact_pairs(pairs)


def _d118_k136(coords: Mapping[int, Sequence[float]]) -> float | None:
    if 118 not in coords or 136 not in coords:
        return None
    return round(_distance(coords[118], coords[136]), 6)


def _structure_paths(directory: Path) -> tuple[Path, ...]:
    allowed = {".pdb", ".ent", ".cif", ".mmcif"}
    return tuple(sorted(path for path in directory.rglob("*") if path.is_file() and path.suffix.lower() in allowed))


def build_conformational_ensemble_packet(
    *,
    structure_dir: Path,
    source_id: str,
    chain_id: str | None = "A",
    reject_alphafold_like_inputs: bool = True,
    minimum_ca_count: int = 150,
) -> ConformationalEnsemblePacket:
    paths = _structure_paths(structure_dir)
    if reject_alphafold_like_inputs and any(_looks_alphafold_like(path, source_id) for path in paths):
        return ConformationalEnsemblePacket(
            kind=CONFORMATIONAL_ENSEMBLE_KIND,
            source_id=source_id,
            structure_file_count=len(paths),
            parsed_conformer_count=0,
            selected_conformer_path=None,
            selected_contact_count=0,
            selected_long_range_contact_count=0,
            selected_contact_map_hash=contact_map_hash(()),
            selection_rule=OPEN_ADK_PROXY_RULE,
            contacts=(),
            conformers=(),
            alphafold_like_input_rejected=True,
        )

    conformers: list[ParsedConformer] = []
    for path in paths:
        coords = parse_ca_coordinates(path, chain_id=chain_id)
        if len(coords) < minimum_ca_count:
            continue
        contacts = contacts_from_ca_coordinates(coords)
        conformers.append(
            ParsedConformer(
                path=str(path),
                format="mmcif" if path.suffix.lower() in {".cif", ".mmcif"} else "pdb",
                chain_id=chain_id or "*",
                ca_count=len(coords),
                contact_count=len(contacts),
                long_range_contact_count=sum(1 for i, j in contacts if j - i >= 24),
                d118_k136_ca_distance_angstrom=_d118_k136(coords),
                contact_map_hash=contact_map_hash(contacts),
                selected_by_open_state_proxy=False,
                contacts=contacts,
            )
        )
    if not conformers:
        return ConformationalEnsemblePacket(
            kind=CONFORMATIONAL_ENSEMBLE_KIND,
            source_id=source_id,
            structure_file_count=len(paths),
            parsed_conformer_count=0,
            selected_conformer_path=None,
            selected_contact_count=0,
            selected_long_range_contact_count=0,
            selected_contact_map_hash=contact_map_hash(()),
            selection_rule=OPEN_ADK_PROXY_RULE,
            contacts=(),
            conformers=(),
        )

    with_d118 = [conformer for conformer in conformers if conformer.d118_k136_ca_distance_angstrom is not None]
    if with_d118:
        selected = min(with_d118, key=lambda conformer: (float(conformer.d118_k136_ca_distance_angstrom or 999.0), -conformer.contact_count, conformer.path))
    else:
        selected = max(conformers, key=lambda conformer: (conformer.contact_count, conformer.long_range_contact_count, conformer.path))
    tagged = tuple(
        ParsedConformer(
            path=conformer.path,
            format=conformer.format,
            chain_id=conformer.chain_id,
            ca_count=conformer.ca_count,
            contact_count=conformer.contact_count,
            long_range_contact_count=conformer.long_range_contact_count,
            d118_k136_ca_distance_angstrom=conformer.d118_k136_ca_distance_angstrom,
            contact_map_hash=conformer.contact_map_hash,
            selected_by_open_state_proxy=(conformer.path == selected.path),
            contacts=conformer.contacts,
        )
        for conformer in conformers
    )
    return ConformationalEnsemblePacket(
        kind=CONFORMATIONAL_ENSEMBLE_KIND,
        source_id=source_id,
        structure_file_count=len(paths),
        parsed_conformer_count=len(tagged),
        selected_conformer_path=selected.path,
        selected_contact_count=len(selected.contacts),
        selected_long_range_contact_count=sum(1 for i, j in selected.contacts if j - i >= 24),
        selected_contact_map_hash=contact_map_hash(selected.contacts),
        selection_rule=OPEN_ADK_PROXY_RULE,
        contacts=tuple(selected.contacts),
        conformers=tagged,
    )
