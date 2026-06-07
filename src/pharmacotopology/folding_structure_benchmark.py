from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, dataclass
from math import sqrt
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from pharmacotopology.folding_reference_loader import load_folding_reference_dataset
from pharmacotopology.folding_topology import (
    BROAD_FOLD_CLASSES,
    CONTACT_PROXY_DIMENSIONS,
    DISORDER_PROMOTING_AMINO_ACIDS,
    FOLDING_TOPOLOGY_DIMENSIONS,
    HYDROPHOBIC_AMINO_ACIDS,
    FoldingReferenceExample,
    FoldingTopologySignature,
    classify_broad_fold_topology,
    contact_map_proxy_similarity,
    normalize_fold_class,
    normalize_sequence,
    predict_topology_signature,
    sequence_features,
    signature_to_dict,
)


STRUCTURE_BENCHMARK_KIND = "structure_derived_topology_benchmark"
LABEL_BENCHMARK_KIND_V0 = "real_external_label_benchmark_v0"
STRUCTURE_SIGNATURE_KIND = "coordinate_contact_graph_extracted"
DISORDER_SIGNATURE_KIND = "disorder_reference_extracted_no_coordinates"

ORDER_CONTROL_KINDS: tuple[str, ...] = (
    "composition_shuffle",
    "reversed_sequence",
    "local_window_shuffle",
    "hydrophobic_cluster_destroyed",
    "charge_pattern_scrambled",
)

AA3_TO_1: dict[str, str] = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


@dataclass(frozen=True)
class StructureResidue:
    index: int
    residue_number: int
    residue_name: str
    amino_acid: str
    x: float
    y: float
    z: float
    b_factor: float
    helix: bool = False
    sheet: bool = False


@dataclass(frozen=True)
class ContactGraphEvidence:
    residues: tuple[StructureResidue, ...]
    contact_pairs: tuple[tuple[int, int], ...]
    helix_fraction: float
    sheet_fraction: float
    bfactor_mean: float


@dataclass(frozen=True)
class StructureEvidenceRow:
    protein_id: str
    source_accession: str
    reference_structure_source: str
    label_fold_class: str
    evidence_kind: str
    structure_topology_signature_kind: str
    structure_fold_class: str
    structure_topology_signature: FoldingTopologySignature
    structure_features: dict[str, float]
    structure_notes: tuple[str, ...]
    clinical_use_allowed: bool = False
    drug_design_created: bool = False
    molecule_generated: bool = False
    protein_sequence_design_created: bool = False
    folding_solution_claim_created: bool = False
    folding_problem_solved: bool = False

    def to_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["structure_topology_signature"] = signature_to_dict(
            self.structure_topology_signature
        )
        return data


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return max(lower, min(upper, value))


def _rounded(value: float) -> float:
    return round(_clamp(value), 6)


def _mean(values: Iterable[float]) -> float:
    values = tuple(values)
    if not values:
        return 0.0
    return sum(values) / len(values)


def _distance(left: StructureResidue, right: StructureResidue) -> float:
    return sqrt(
        (left.x - right.x) ** 2
        + (left.y - right.y) ** 2
        + (left.z - right.z) ** 2
    )


def _parse_secondary_ranges(
    pdb_text: str,
    *,
    chain_id: str,
) -> tuple[set[int], set[int]]:
    helix: set[int] = set()
    sheet: set[int] = set()
    for line in pdb_text.splitlines():
        fields = line.split()
        try:
            if line.startswith("HELIX") and len(fields) >= 9:
                start_chain = fields[4]
                start = int(fields[5])
                end_chain = fields[7]
                end = int(fields[8])
                if start_chain == chain_id and end_chain == chain_id:
                    helix.update(range(min(start, end), max(start, end) + 1))
            elif line.startswith("SHEET") and len(fields) >= 10:
                start_chain = fields[5]
                start = int(fields[6])
                end_chain = fields[8]
                end = int(fields[9])
                if start_chain == chain_id and end_chain == chain_id:
                    sheet.update(range(min(start, end), max(start, end) + 1))
        except (IndexError, ValueError):
            continue
    return helix, sheet


def parse_pdb_contact_graph(
    pdb_text: str,
    *,
    chain_id: str,
    contact_cutoff_angstrom: float = 8.0,
    minimum_sequence_separation: int = 3,
) -> ContactGraphEvidence:
    helix_numbers, sheet_numbers = _parse_secondary_ranges(
        pdb_text,
        chain_id=chain_id,
    )
    residues: list[StructureResidue] = []
    seen: set[tuple[str, int, str]] = set()
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL") and residues:
            break
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        alternate = line[16:17].strip()
        chain = line[21:22].strip()
        if atom_name != "CA" or chain != chain_id or alternate not in ("", "A"):
            continue
        try:
            residue_number = int(line[22:26])
            insertion_code = line[26:27].strip()
            residue_name = line[17:20].strip()
            amino_acid = AA3_TO_1.get(residue_name, "X")
            key = (chain, residue_number, insertion_code)
            if key in seen or amino_acid == "X":
                continue
            seen.add(key)
            residues.append(
                StructureResidue(
                    index=len(residues),
                    residue_number=residue_number,
                    residue_name=residue_name,
                    amino_acid=amino_acid,
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    b_factor=float((line[60:66] or "0").strip() or 0.0),
                    helix=residue_number in helix_numbers,
                    sheet=residue_number in sheet_numbers,
                )
            )
        except ValueError:
            continue

    contact_pairs: list[tuple[int, int]] = []
    for left_index, left in enumerate(residues):
        for right_index in range(left_index + minimum_sequence_separation, len(residues)):
            right = residues[right_index]
            if _distance(left, right) <= contact_cutoff_angstrom:
                contact_pairs.append((left_index, right_index))

    residue_count = max(len(residues), 1)
    return ContactGraphEvidence(
        residues=tuple(residues),
        contact_pairs=tuple(contact_pairs),
        helix_fraction=_rounded(
            sum(1 for residue in residues if residue.helix) / residue_count
        ),
        sheet_fraction=_rounded(
            sum(1 for residue in residues if residue.sheet) / residue_count
        ),
        bfactor_mean=round(_mean(residue.b_factor for residue in residues), 6),
    )


def _radius_of_gyration(residues: Sequence[StructureResidue]) -> float:
    if not residues:
        return 0.0
    center_x = _mean(residue.x for residue in residues)
    center_y = _mean(residue.y for residue in residues)
    center_z = _mean(residue.z for residue in residues)
    return sqrt(
        _mean(
            (residue.x - center_x) ** 2
            + (residue.y - center_y) ** 2
            + (residue.z - center_z) ** 2
            for residue in residues
        )
    )


def _contact_crossing_fraction(contact_pairs: Sequence[tuple[int, int]]) -> float:
    if len(contact_pairs) < 2:
        return 0.0
    crossings = 0
    reviewed = 0
    for pair_index, (left_a, left_b) in enumerate(contact_pairs):
        for right_a, right_b in contact_pairs[pair_index + 1 :]:
            reviewed += 1
            if left_a < right_a < left_b < right_b:
                crossings += 1
    return crossings / reviewed if reviewed else 0.0


def structure_signature_from_contact_graph(
    graph: ContactGraphEvidence,
    *,
    expected_sequence_length: int,
    reference_sequence: str,
) -> tuple[FoldingTopologySignature, dict[str, float]]:
    sequence = normalize_sequence(reference_sequence)
    features = sequence_features(sequence)
    residues = graph.residues
    residue_count = len(residues)
    possible_contacts = max(
        (residue_count - 3) * (residue_count - 2) / 2,
        1.0,
    )
    contact_count = len(graph.contact_pairs)
    contact_density = contact_count / possible_contacts
    long_contacts = [
        right - left
        for left, right in graph.contact_pairs
        if right - left >= 24
    ]
    average_contact_order = _mean(
        (right - left) / max(residue_count, 1)
        for left, right in graph.contact_pairs
    )
    hydrophobic_indices = {
        residue.index
        for residue in residues
        if residue.amino_acid in HYDROPHOBIC_AMINO_ACIDS
    }
    hydrophobic_contact_count = sum(
        1
        for left, right in graph.contact_pairs
        if left in hydrophobic_indices and right in hydrophobic_indices
    )
    hydrophobic_possible = max(len(hydrophobic_indices) ** 2 / 2, 1.0)
    hydrophobic_contact_density = hydrophobic_contact_count / hydrophobic_possible
    contact_degrees = {index: 0 for index in range(residue_count)}
    for left, right in graph.contact_pairs:
        contact_degrees[left] += 1
        contact_degrees[right] += 1
    hydrophobic_core_fraction = (
        sum(
            1
            for index in hydrophobic_indices
            if contact_degrees.get(index, 0) >= 4
        )
        / max(len(hydrophobic_indices), 1)
    )
    midpoint = residue_count // 2
    left_internal = sum(
        1 for left, right in graph.contact_pairs if left < midpoint and right < midpoint
    )
    right_internal = sum(
        1 for left, right in graph.contact_pairs if left >= midpoint and right >= midpoint
    )
    bridge_contacts = sum(
        1 for left, right in graph.contact_pairs if left < midpoint <= right
    )
    internal_contacts = max(left_internal + right_internal, 1)
    bridge_ratio = bridge_contacts / internal_contacts
    missing_fraction = max(expected_sequence_length - residue_count, 0) / max(
        expected_sequence_length,
        1,
    )
    bfactor_pressure = _clamp((graph.bfactor_mean - 20.0) / 60.0)
    loop_fraction = _clamp(1.0 - graph.helix_fraction - graph.sheet_fraction)
    radius_of_gyration = _radius_of_gyration(residues)
    compactness = _clamp(1.0 - (radius_of_gyration / max(residue_count, 1) ** (1 / 3)) / 6.0)
    domain_boundary_signal = _clamp(1.0 - bridge_ratio * 1.8)
    crossing_fraction = _contact_crossing_fraction(graph.contact_pairs)

    structure_features = {
        "residue_count": float(residue_count),
        "coordinate_coverage": _rounded(residue_count / max(expected_sequence_length, 1)),
        "contact_count": float(contact_count),
        "contact_density": round(contact_density, 6),
        "long_range_contact_fraction": _rounded(
            len(long_contacts) / max(contact_count, 1)
        ),
        "average_contact_order": round(average_contact_order, 6),
        "hydrophobic_contact_density": round(hydrophobic_contact_density, 6),
        "hydrophobic_core_fraction": _rounded(hydrophobic_core_fraction),
        "domain_boundary_signal": _rounded(domain_boundary_signal),
        "helix_fraction": graph.helix_fraction,
        "sheet_fraction": graph.sheet_fraction,
        "loop_fraction": _rounded(loop_fraction),
        "bfactor_mean": graph.bfactor_mean,
        "bfactor_flexibility_pressure": _rounded(bfactor_pressure),
        "missing_residue_fraction": _rounded(missing_fraction),
        "compactness": _rounded(compactness),
        "radius_of_gyration": round(radius_of_gyration, 6),
        "contact_crossing_fraction": round(crossing_fraction, 6),
    }
    signature = FoldingTopologySignature(
        sequence_complexity=float(features["sequence_complexity"]),
        secondary_structure_balance=_rounded(
            1.0 - abs(graph.helix_fraction - graph.sheet_fraction)
        ),
        contact_map_closure=_rounded(contact_density / 0.12),
        hydrophobic_core_closure=_rounded(
            hydrophobic_contact_density / 0.12 * 0.45
            + hydrophobic_core_fraction * 0.55
        ),
        loop_disorder_pressure=_rounded(
            loop_fraction * 0.45 + bfactor_pressure * 0.25 + missing_fraction * 0.30
        ),
        domain_boundary_stability=_rounded(1.0 - domain_boundary_signal),
        long_range_contact_order=_rounded(
            average_contact_order * 1.35
            + len(long_contacts) / max(contact_count, 1) * 0.35
        ),
        conformational_flexibility=_rounded(
            loop_fraction * 0.35 + bfactor_pressure * 0.35 + missing_fraction * 0.30
        ),
        knot_or_entanglement_signature=_rounded(
            crossing_fraction * 3.0
            + structure_features["long_range_contact_fraction"] * 0.25
            + average_contact_order * 0.35
        ),
        uncertainty_radius=_rounded(
            0.18 + missing_fraction * 0.35 + bfactor_pressure * 0.20
        ),
    )
    return signature, structure_features


def structure_fold_class_from_evidence(
    signature: FoldingTopologySignature,
    features: Mapping[str, float],
    *,
    evidence_kind: str,
) -> str:
    if evidence_kind == DISORDER_SIGNATURE_KIND:
        return "disordered_flexible"
    helix_fraction = float(features.get("helix_fraction", 0.0))
    sheet_fraction = float(features.get("sheet_fraction", 0.0))
    boundary_signal = float(features.get("domain_boundary_signal", 0.0))
    compactness = float(features.get("compactness", 0.0))
    sequence_length = float(features.get("residue_count", 0.0))
    if helix_fraction >= 0.42 and sheet_fraction < 0.20 and compactness >= 0.48:
        return "alpha_rich"
    if sheet_fraction >= 0.30 and helix_fraction < 0.25:
        return "beta_rich"
    if helix_fraction >= 0.15 and sheet_fraction >= 0.15:
        return "alpha_beta_mixed"
    if (
        boundary_signal >= 0.70
        and sequence_length >= 120
        and compactness < 0.50
    ):
        return "multidomain_boundary"
    return classify_broad_fold_topology(signature)


def structure_evidence_from_pdb_text(
    reference: FoldingReferenceExample,
    pdb_text: str,
    *,
    chain_id: str,
    coordinate_source: str,
) -> StructureEvidenceRow:
    graph = parse_pdb_contact_graph(pdb_text, chain_id=chain_id)
    signature, features = structure_signature_from_contact_graph(
        graph,
        expected_sequence_length=reference.sequence_length
        or len(normalize_sequence(reference.sequence)),
        reference_sequence=reference.sequence,
    )
    structure_class = structure_fold_class_from_evidence(
        signature,
        features,
        evidence_kind=STRUCTURE_SIGNATURE_KIND,
    )
    return StructureEvidenceRow(
        protein_id=reference.protein_id,
        source_accession=reference.source_accession,
        reference_structure_source=reference.reference_structure_source,
        label_fold_class=normalize_fold_class(reference.reference_fold_class),
        evidence_kind="coordinate_contact_graph",
        structure_topology_signature_kind=STRUCTURE_SIGNATURE_KIND,
        structure_fold_class=structure_class,
        structure_topology_signature=signature,
        structure_features=features,
        structure_notes=(
            f"derived from C-alpha contacts in {coordinate_source}",
            "labels are not used to compute the coordinate-derived signature",
        ),
    )


def disorder_evidence_from_reference(
    reference: FoldingReferenceExample,
) -> StructureEvidenceRow:
    sequence = normalize_sequence(reference.sequence)
    features = sequence_features(sequence)
    disorder_fraction = sum(
        1 for residue in sequence if residue in DISORDER_PROMOTING_AMINO_ACIDS
    ) / len(sequence)
    signature = FoldingTopologySignature(
        sequence_complexity=float(features["sequence_complexity"]),
        secondary_structure_balance=0.28,
        contact_map_closure=0.06,
        hydrophobic_core_closure=0.08,
        loop_disorder_pressure=_rounded(0.72 + disorder_fraction * 0.28),
        domain_boundary_stability=0.20,
        long_range_contact_order=0.05,
        conformational_flexibility=0.90,
        knot_or_entanglement_signature=0.03,
        uncertainty_radius=0.65,
    )
    return StructureEvidenceRow(
        protein_id=reference.protein_id,
        source_accession=reference.source_accession,
        reference_structure_source=reference.reference_structure_source,
        label_fold_class=normalize_fold_class(reference.reference_fold_class),
        evidence_kind="disorder_reference",
        structure_topology_signature_kind=DISORDER_SIGNATURE_KIND,
        structure_fold_class="disordered_flexible",
        structure_topology_signature=signature,
        structure_features={
            "residue_count": float(len(sequence)),
            "coordinate_coverage": 0.0,
            "disprot_disorder_reference": 1.0,
            "disorder_promoting_fraction": _rounded(disorder_fraction),
            "contact_density": 0.0,
            "domain_boundary_signal": 0.0,
            "helix_fraction": 0.0,
            "sheet_fraction": 0.0,
            "loop_fraction": 1.0,
            "compactness": 0.0,
        },
        structure_notes=(
            "derived from DisProt disorder reference channel, not coordinates",
            "kept separate from folded PDB coordinate rows",
        ),
    )


def pdb_id_and_chain(reference: FoldingReferenceExample) -> tuple[str, str]:
    source = reference.reference_structure_source
    if source.lower().startswith("pdb:"):
        parts = source.split(":")
        if len(parts) >= 3:
            return parts[1].upper(), parts[2]
    accession_parts = reference.source_accession.split(":")
    if len(accession_parts) >= 2:
        return accession_parts[0].upper(), accession_parts[1]
    raise ValueError(f"{reference.protein_id} does not include a PDB chain id.")


def build_structure_evidence_rows(
    references: Sequence[FoldingReferenceExample],
    *,
    pdb_dir: Path,
) -> tuple[StructureEvidenceRow, ...]:
    rows: list[StructureEvidenceRow] = []
    for reference in references:
        if reference.reference_structure_source.lower().startswith("disprot:"):
            rows.append(disorder_evidence_from_reference(reference))
            continue
        pdb_id, chain_id = pdb_id_and_chain(reference)
        pdb_path = pdb_dir / f"{pdb_id}.pdb"
        if not pdb_path.exists():
            raise FileNotFoundError(f"Missing PDB coordinate file: {pdb_path}")
        rows.append(
            structure_evidence_from_pdb_text(
                reference,
                pdb_path.read_text(encoding="utf-8"),
                chain_id=chain_id,
                coordinate_source=str(pdb_path),
            )
        )
    return tuple(rows)


def evidence_rows_to_dicts(
    rows: Sequence[StructureEvidenceRow],
) -> list[dict[str, object]]:
    return [row.to_dict() for row in rows]


def write_structure_evidence_file(
    rows: Sequence[StructureEvidenceRow],
    output_path: Path,
    *,
    source_benchmark_file: Path,
) -> Path:
    payload = {
        "benchmark_kind": STRUCTURE_BENCHMARK_KIND,
        "source_benchmark_file": str(source_benchmark_file),
        "reference_topology_source": "coordinate_contact_graph_or_disorder_reference",
        "coordinate_rows": sum(
            1 for row in rows if row.evidence_kind == "coordinate_contact_graph"
        ),
        "disorder_reference_rows": sum(
            1 for row in rows if row.evidence_kind == "disorder_reference"
        ),
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "references": evidence_rows_to_dicts(rows),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output_path


def load_structure_evidence_rows(path: Path) -> tuple[StructureEvidenceRow, ...]:
    data = json.loads(path.read_text(encoding="utf-8"))
    raw_rows = data.get("references") if isinstance(data, Mapping) else None
    if not isinstance(raw_rows, list):
        raise ValueError("Structure evidence file must include a references list.")
    rows: list[StructureEvidenceRow] = []
    for row in raw_rows:
        if not isinstance(row, Mapping):
            raise ValueError("Each structure evidence row must be an object.")
        signature_data = row.get("structure_topology_signature")
        if not isinstance(signature_data, Mapping):
            raise ValueError("structure_topology_signature must be a mapping.")
        signature = FoldingTopologySignature(
            **{
                dimension: float(signature_data[dimension])
                for dimension in FOLDING_TOPOLOGY_DIMENSIONS
            }
        )
        features = row.get("structure_features", {})
        notes = row.get("structure_notes", ())
        rows.append(
            StructureEvidenceRow(
                protein_id=str(row["protein_id"]),
                source_accession=str(row.get("source_accession", "")),
                reference_structure_source=str(row["reference_structure_source"]),
                label_fold_class=str(row["label_fold_class"]),
                evidence_kind=str(row["evidence_kind"]),
                structure_topology_signature_kind=str(
                    row["structure_topology_signature_kind"]
                ),
                structure_fold_class=str(row["structure_fold_class"]),
                structure_topology_signature=signature,
                structure_features={
                    str(key): float(value)
                    for key, value in dict(features).items()
                    if isinstance(value, (int, float))
                },
                structure_notes=tuple(str(note) for note in notes),
            )
        )
    return tuple(rows)


def _signature_delta_mean(
    left: FoldingTopologySignature,
    right: FoldingTopologySignature,
) -> float:
    left_values = signature_to_dict(left)
    right_values = signature_to_dict(right)
    return round(
        sum(
            abs(left_values[dimension] - right_values[dimension])
            for dimension in FOLDING_TOPOLOGY_DIMENSIONS
        )
        / len(FOLDING_TOPOLOGY_DIMENSIONS),
        6,
    )


def _control_sequence(sequence: str, control_kind: str, *, protein_id: str) -> str:
    normalized = normalize_sequence(sequence)
    if control_kind == "reversed_sequence":
        return normalized[::-1]
    if control_kind == "composition_shuffle":
        return "".join(
            residue
            for _, residue in sorted(
                (
                    hashlib.sha256(
                        f"{protein_id}:composition:{index}:{residue}".encode("utf-8")
                    ).hexdigest(),
                    residue,
                )
                for index, residue in enumerate(normalized)
            )
        )
    if control_kind == "local_window_shuffle":
        windows = []
        for start in range(0, len(normalized), 9):
            window = normalized[start : start + 9]
            windows.append(
                "".join(
                    residue
                    for _, residue in sorted(
                        (
                            hashlib.sha256(
                                f"{protein_id}:window:{start}:{index}:{residue}".encode(
                                    "utf-8"
                                )
                            ).hexdigest(),
                            residue,
                        )
                        for index, residue in enumerate(window)
                    )
                )
            )
        return "".join(windows)
    if control_kind == "hydrophobic_cluster_destroyed":
        hydrophobic = [
            residue for residue in normalized if residue in HYDROPHOBIC_AMINO_ACIDS
        ]
        other = [
            residue for residue in normalized if residue not in HYDROPHOBIC_AMINO_ACIDS
        ]
        output: list[str] = []
        while hydrophobic or other:
            if other:
                output.append(other.pop(0))
            if hydrophobic:
                output.append(hydrophobic.pop(0))
        return "".join(output)
    if control_kind == "charge_pattern_scrambled":
        charged = [residue for residue in normalized if residue in "DEKRH"]
        other = [residue for residue in normalized if residue not in "DEKRH"]
        output = []
        for index in range(len(normalized)):
            if index % 3 == 0 and charged:
                output.append(charged.pop(0))
            elif other:
                output.append(other.pop(0))
            elif charged:
                output.append(charged.pop(0))
        return "".join(output)
    raise ValueError(f"Unknown order control kind: {control_kind}")


def sequence_order_control_rows(
    references: Sequence[FoldingReferenceExample],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for reference in references:
        real_sequence = normalize_sequence(reference.sequence)
        real_signature = predict_topology_signature(real_sequence)
        real_fold_class = classify_broad_fold_topology(real_signature)
        for control_kind in ORDER_CONTROL_KINDS:
            control_sequence = _control_sequence(
                real_sequence,
                control_kind,
                protein_id=reference.protein_id,
            )
            control_signature = predict_topology_signature(control_sequence)
            control_fold_class = classify_broad_fold_topology(control_signature)
            delta = _signature_delta_mean(real_signature, control_signature)
            rows.append(
                {
                    "protein_id": reference.protein_id,
                    "control_kind": control_kind,
                    "sequence_length": len(real_sequence),
                    "same_composition": sorted(real_sequence) == sorted(control_sequence),
                    "real_predicted_fold_class": real_fold_class,
                    "control_predicted_fold_class": control_fold_class,
                    "fold_class_changed": real_fold_class != control_fold_class,
                    "signature_delta_mean": delta,
                    "control_sequence_written": False,
                    "control_purpose": (
                        "internal falsification control; sequence order is perturbed "
                        "without emitting a designed protein sequence"
                    ),
                }
            )
    return rows


def structure_benchmark_rows(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
) -> list[dict[str, object]]:
    evidence_by_id = {row.protein_id: row for row in evidence_rows}
    rows: list[dict[str, object]] = []
    for reference in references:
        evidence = evidence_by_id[reference.protein_id]
        predicted = predict_topology_signature(reference.sequence)
        predicted_class = classify_broad_fold_topology(predicted)
        label_class = normalize_fold_class(reference.reference_fold_class)
        structure_similarity = contact_map_proxy_similarity(
            predicted,
            evidence.structure_topology_signature,
            CONTACT_PROXY_DIMENSIONS,
        )
        label_similarity = contact_map_proxy_similarity(
            predicted,
            reference.reference_topology_signature,
            CONTACT_PROXY_DIMENSIONS,
        )
        structure_label_similarity = contact_map_proxy_similarity(
            evidence.structure_topology_signature,
            reference.reference_topology_signature,
            CONTACT_PROXY_DIMENSIONS,
        )
        structure_class_match = predicted_class == evidence.structure_fold_class
        label_class_match = predicted_class == label_class
        structure_label_agreement = evidence.structure_fold_class == label_class
        failure_reason = ""
        if not structure_class_match:
            failure_reason = "prediction_structure_class_mismatch"
        elif structure_similarity < 0.70:
            failure_reason = "low_prediction_structure_similarity"
        rows.append(
            {
                "protein_id": reference.protein_id,
                "sequence_length": len(normalize_sequence(reference.sequence)),
                "source_accession": reference.source_accession,
                "reference_structure_source": reference.reference_structure_source,
                "evidence_kind": evidence.evidence_kind,
                "predicted_fold_class": predicted_class,
                "structure_fold_class": evidence.structure_fold_class,
                "label_fold_class": label_class,
                "prediction_vs_structure_score": structure_similarity,
                "prediction_vs_label_score": label_similarity,
                "structure_vs_label_score": structure_label_similarity,
                "prediction_structure_class_match": structure_class_match,
                "prediction_label_class_match": label_class_match,
                "structure_vs_label_agreement": structure_label_agreement,
                "predicted_topology_signature": signature_to_dict(predicted),
                "structure_topology_signature": signature_to_dict(
                    evidence.structure_topology_signature
                ),
                "label_topology_signature": signature_to_dict(
                    reference.reference_topology_signature
                ),
                "structure_features": evidence.structure_features,
                "structure_topology_signature_kind": (
                    evidence.structure_topology_signature_kind
                ),
                "failure_reason": failure_reason,
                "folding_problem_solved": False,
                "folding_solution_claim_created": False,
            }
        )
    return rows


def sequence_order_sensitivity_score(
    control_rows: Sequence[Mapping[str, object]],
) -> float:
    if not control_rows:
        return 0.0
    return round(
        sum(float(row.get("signature_delta_mean", 0.0)) for row in control_rows)
        / len(control_rows),
        6,
    )


def build_structure_benchmark_report(
    references: Sequence[FoldingReferenceExample],
    evidence_rows: Sequence[StructureEvidenceRow],
    *,
    source_benchmark_file: Path,
    structure_evidence_file: Path,
) -> dict[str, object]:
    rows = structure_benchmark_rows(references, evidence_rows)
    control_rows = sequence_order_control_rows(references)
    structure_matches = sum(
        1 for row in rows if bool(row["prediction_structure_class_match"])
    )
    label_matches = sum(1 for row in rows if bool(row["prediction_label_class_match"]))
    structure_label_agreements = sum(
        1 for row in rows if bool(row["structure_vs_label_agreement"])
    )
    sensitivity = sequence_order_sensitivity_score(control_rows)
    coordinate_rows = sum(1 for row in rows if row["evidence_kind"] == "coordinate_contact_graph")
    disorder_rows = sum(1 for row in rows if row["evidence_kind"] == "disorder_reference")
    revision_required = (
        structure_matches < len(rows)
        or sensitivity < 0.03
        or any(row["failure_reason"] for row in rows)
    )
    return {
        "benchmark_kind": STRUCTURE_BENCHMARK_KIND,
        "source_label_benchmark_kind": LABEL_BENCHMARK_KIND_V0,
        "source_benchmark_file": str(source_benchmark_file),
        "structure_evidence_file": str(structure_evidence_file),
        "reference_topology_source": "coordinate_contact_graph_or_disorder_reference",
        "benchmark_size": len(rows),
        "coordinate_rows": coordinate_rows,
        "disorder_reference_rows": disorder_rows,
        "prediction_vs_structure_match_count": structure_matches,
        "prediction_vs_structure_mismatch_count": len(rows) - structure_matches,
        "prediction_vs_structure_accuracy": round(structure_matches / max(len(rows), 1), 6),
        "prediction_vs_label_match_count": label_matches,
        "prediction_vs_label_accuracy": round(label_matches / max(len(rows), 1), 6),
        "structure_vs_label_agreement_count": structure_label_agreements,
        "structure_vs_label_agreement_rate": round(
            structure_label_agreements / max(len(rows), 1),
            6,
        ),
        "mean_prediction_vs_structure_score": round(
            _mean(float(row["prediction_vs_structure_score"]) for row in rows),
            6,
        ),
        "mean_prediction_vs_label_score": round(
            _mean(float(row["prediction_vs_label_score"]) for row in rows),
            6,
        ),
        "sequence_order_sensitivity_score": sensitivity,
        "composition_only_warning": sensitivity < 0.03,
        "contradiction_against_stability_seen": revision_required,
        "revision_candidate_created": revision_required,
        "revision_required": revision_required,
        "claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "boundary_statement": (
            "This report compares sequence-derived topology predictions with "
            "coordinate-derived structure evidence where PDB coordinates are "
            "available, and with disorder-reference evidence for DisProt rows. "
            "It is not a folding solution."
        ),
        "structure_rows": rows,
    }


def build_falsification_report(
    structure_report: Mapping[str, object],
    control_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    unchanged_controls = sum(
        1 for row in control_rows if float(row.get("signature_delta_mean", 0.0)) == 0.0
    )
    sensitivity = sequence_order_sensitivity_score(control_rows)
    revision_required = bool(structure_report.get("revision_required", True))
    return {
        "benchmark_kind": "knot_style_folding_falsification_summary",
        "source_benchmark_kind": STRUCTURE_BENCHMARK_KIND,
        "lambda_local_signal": True,
        "delta_comparison": True,
        "tau_transfer_falsification": sensitivity >= 0.03 and unchanged_controls == 0,
        "sigma_stability": False,
        "chi_contradiction_seen": revision_required,
        "rho_revision_needed": revision_required,
        "omega_terminal_claim_allowed": False,
        "sequence_order_sensitivity_score": sensitivity,
        "order_controls_reviewed": len(control_rows),
        "controls_with_no_signature_change": unchanged_controls,
        "composition_only_warning": sensitivity < 0.03,
        "contradiction_against_stability_seen": revision_required,
        "revision_candidate_created": revision_required,
        "claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "recommended_revision": (
            "Add sequence-order, local-window, and contact-aware terms before "
            "treating broad fold-class matches as stability evidence."
        ),
    }


def write_structure_outputs(
    *,
    report: Mapping[str, object],
    control_rows: Sequence[Mapping[str, object]],
    report_path: Path,
    rows_path: Path,
    dashboard_path: Path,
    order_controls_path: Path,
    falsification_report_path: Path,
) -> tuple[Path, Path, Path, Path, Path]:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    rows = list(report.get("structure_rows", []))
    _write_csv_rows(rows, rows_path, nested_json_fields={
        "predicted_topology_signature",
        "structure_topology_signature",
        "label_topology_signature",
        "structure_features",
    })
    _write_csv_rows(control_rows, order_controls_path, nested_json_fields=set())
    falsification = build_falsification_report(report, control_rows)
    falsification_report_path.write_text(
        json.dumps(falsification, ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    dashboard_path.write_text(
        render_structure_dashboard(report, falsification),
        encoding="utf-8",
    )
    return (
        report_path,
        rows_path,
        dashboard_path,
        order_controls_path,
        falsification_report_path,
    )


def _write_csv_rows(
    rows: Sequence[Mapping[str, object]],
    path: Path,
    *,
    nested_json_fields: set[str],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as file:
        if not rows:
            file.write("")
            return path
        fieldnames = list(rows[0])
        writer = csv.DictWriter(file, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            output = {}
            for key, value in row.items():
                if key in nested_json_fields:
                    output[key] = json.dumps(value, sort_keys=True, separators=(",", ":"))
                else:
                    output[key] = value
            writer.writerow(output)
    return path


def _escape(value: object) -> str:
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def render_structure_dashboard(
    report: Mapping[str, object],
    falsification_report: Mapping[str, object],
) -> str:
    rows = [
        row for row in report.get("structure_rows", []) if isinstance(row, Mapping)
    ]
    body_rows = []
    for row in rows:
        body_rows.append(
            "<tr>"
            f"<td>{_escape(row.get('protein_id', ''))}</td>"
            f"<td>{_escape(row.get('evidence_kind', ''))}</td>"
            f"<td>{_escape(row.get('predicted_fold_class', ''))}</td>"
            f"<td>{_escape(row.get('structure_fold_class', ''))}</td>"
            f"<td>{_escape(row.get('label_fold_class', ''))}</td>"
            f"<td>{_escape(row.get('prediction_vs_structure_score', ''))}</td>"
            f"<td>{_escape(row.get('failure_reason', ''))}</td>"
            "</tr>"
        )
    metric_cards = {
        "Rows": report.get("benchmark_size", 0),
        "Coordinate Rows": report.get("coordinate_rows", 0),
        "Pred vs Structure": report.get("prediction_vs_structure_accuracy", 0.0),
        "Pred vs Label": report.get("prediction_vs_label_accuracy", 0.0),
        "Order Sensitivity": report.get("sequence_order_sensitivity_score", 0.0),
        "Claim Allowed": report.get("claim_allowed", False),
    }
    metrics = "".join(
        "<div class=\"metric\">"
        f"<span>{_escape(label)}</span><strong>{_escape(value)}</strong>"
        "</div>"
        for label, value in metric_cards.items()
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Structure-Derived Folding Topology Benchmark</title>
<style>
body {{ margin: 0; font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #17211c; background: #fff; }}
header {{ padding: 28px clamp(18px, 4vw, 54px); background: #eef5f2; border-bottom: 1px solid #cad6d0; }}
main {{ padding: 24px clamp(18px, 4vw, 54px) 52px; }}
h1 {{ margin: 0 0 10px; font-size: clamp(28px, 4vw, 44px); letter-spacing: 0; }}
h2 {{ font-size: 22px; letter-spacing: 0; }}
p {{ max-width: 980px; line-height: 1.55; }}
.warning {{ display: inline-block; padding: 8px 10px; border: 1px solid #9b3d2e; color: #9b3d2e; background: #fff8f6; font-weight: 700; }}
.metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 10px; margin-top: 18px; }}
.metric {{ border: 1px solid #cad6d0; background: #f7faf8; padding: 12px; min-height: 82px; }}
.metric span {{ display: block; color: #5d6762; font-size: 13px; }}
.metric strong {{ display: block; margin-top: 8px; font-size: 22px; }}
table {{ width: 100%; border-collapse: collapse; font-size: 14px; }}
th, td {{ border: 1px solid #cad6d0; padding: 8px; text-align: left; vertical-align: top; }}
th {{ background: #eef5f2; }}
code {{ background: #f0f3f1; padding: 1px 4px; }}
</style>
</head>
<body>
<header>
  <p class="warning">NOT A FOLDING SOLUTION / INTERNAL FALSIFICATION REVIEW</p>
  <h1>Structure-Derived Folding Topology Benchmark</h1>
  <p>This dashboard separates sequence-derived prediction, coordinate-derived
  topology evidence, and external label agreement. It keeps the first 10-row
  result as a contradiction/revision signal rather than a terminal claim.</p>
  <div class="metrics">{metrics}</div>
</header>
<main>
  <section>
    <h2>Falsification Summary</h2>
    <p><code>composition_only_warning</code>: {_escape(falsification_report.get('composition_only_warning', False))}</p>
    <p><code>rho_revision_needed</code>: {_escape(falsification_report.get('rho_revision_needed', True))}</p>
    <p><code>omega_terminal_claim_allowed</code>: {_escape(falsification_report.get('omega_terminal_claim_allowed', False))}</p>
  </section>
  <section>
    <h2>Prediction vs Structure Rows</h2>
    <table>
      <thead><tr><th>protein_id</th><th>evidence</th><th>predicted</th><th>structure</th><th>label</th><th>score</th><th>failure</th></tr></thead>
      <tbody>{''.join(body_rows)}</tbody>
    </table>
  </section>
</main>
</body>
</html>
"""


def load_real10_with_structure_evidence(
    benchmark_file: Path,
    structure_evidence_file: Path,
) -> tuple[tuple[FoldingReferenceExample, ...], tuple[StructureEvidenceRow, ...]]:
    dataset = load_folding_reference_dataset(benchmark_file, require_external=True)
    evidence = load_structure_evidence_rows(structure_evidence_file)
    expected = {reference.protein_id for reference in dataset.references}
    observed = {row.protein_id for row in evidence}
    missing = sorted(expected.difference(observed))
    if missing:
        raise ValueError("Missing structure evidence rows: " + ", ".join(missing))
    unknown = sorted(observed.difference(expected))
    if unknown:
        raise ValueError("Unknown structure evidence rows: " + ", ".join(unknown))
    return dataset.references, evidence
