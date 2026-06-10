from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_contact_topology import sha256_sequence  # noqa: E402
from pharmacotopology.folding_native_contact_eval import contact_map_hash  # noqa: E402
from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    CONTACT_CUTOFF_ANGSTROM,
    MIN_SEQUENCE_SEPARATION,
    REAL_COORDINATE_NATIVE_KIND,
    REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
    coordinate_native_contact_pairs,
    coordinate_trace_hash,
    parse_pdb_ca_coordinate_points,
)
from pharmacotopology.folding_structure_benchmark import AA3_TO_1  # noqa: E402
from pharmacotopology.folding_topology import normalize_sequence  # noqa: E402


DEFAULT_PDB_FILE = Path("data/rcsb_pdb/1UBQ.pdb")
DEFAULT_OUTPUT = Path("data/folding_real_coordinate_holdout_1ubq.locked.json")


def _axis_for_reference_class(reference_fold_class: str) -> dict[str, str]:
    if reference_fold_class == "alpha_rich":
        secondary = "alpha_rich"
        architecture = "compact_single_domain"
    elif reference_fold_class == "beta_rich":
        secondary = "beta_rich"
        architecture = "compact_single_domain"
    elif reference_fold_class == "alpha_beta_mixed":
        secondary = "alpha_beta_mixed"
        architecture = "compact_single_domain"
    elif reference_fold_class == "multidomain_boundary":
        secondary = "alpha_beta_mixed"
        architecture = "multidomain_or_segmented"
    else:
        secondary = "weak_or_unknown"
        architecture = "unknown"
    return {
        "secondary_structure_axis": secondary,
        "architecture_axis": architecture,
        "order_axis": "ordered",
        "environment_axis": "soluble_like",
    }


def _sequence_from_pdb_ca(pdb_text: str, *, chain_id: str) -> str:
    residues: list[str] = []
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
        residue_name = line[17:20].strip()
        residue = AA3_TO_1.get(residue_name, "")
        if not residue:
            continue
        try:
            residue_number = int(line[22:26])
        except ValueError:
            continue
        insertion_code = line[26:27].strip()
        key = (chain, residue_number, insertion_code)
        if key in seen:
            continue
        seen.add(key)
        residues.append(residue)
    return normalize_sequence("".join(residues))


def build_rcsb_single_protein_coordinate_holdout(
    *,
    pdb_file: Path,
    output: Path,
    pdb_id: str,
    chain_id: str,
    source_id: str,
    reference_fold_class: str,
) -> Path:
    pdb_text = pdb_file.read_text(encoding="utf-8")
    sequence = _sequence_from_pdb_ca(pdb_text, chain_id=chain_id)
    points = parse_pdb_ca_coordinate_points(pdb_text, chain_id=chain_id)
    native_pairs = coordinate_native_contact_pairs(points)
    accession = f"{pdb_id.upper()}:{chain_id}"
    row_id = f"coord_holdout_{pdb_id.lower()}_{chain_id.lower()}"
    payload: dict[str, object] = {
        "benchmark_kind": REAL_COORDINATE_VISUAL_BENCHMARK_KIND,
        "holdout_split": "known_protein_single_holdout_v0",
        "source_coordinate_database": "RCSB_PDB",
        "coordinate_fixture_kind": "locked_c_alpha_coordinate_trace_v1",
        "native_contact_derivation_kind": REAL_COORDINATE_NATIVE_KIND,
        "native_contact_derivation_rule": (
            "C-alpha distance <= 8.0 angstrom with sequence separation >= 3"
        ),
        "contact_cutoff_angstrom": CONTACT_CUTOFF_ANGSTROM,
        "minimum_sequence_separation": MIN_SEQUENCE_SEPARATION,
        "benchmark_size": 1,
        "locked_after_generation": True,
        "blind_prediction_before_coordinate_scoring": True,
        "toy_locked_contact_targets_used": False,
        "coarse_ca_only": True,
        "full_atomic_folding_available": False,
        "mechanism_discovery_claim_allowed": False,
        "global_folding_claim_allowed": False,
        "folding_problem_solved": False,
        "folding_solution_claim_created": False,
        "clinical_use_allowed": False,
        "drug_design_created": False,
        "molecule_generated": False,
        "protein_sequence_design_created": False,
        "references": [
            {
                "row_id": row_id,
                "source_id": source_id,
                "source_accession": accession,
                "source_database": "RCSB_PDB",
                "reference_structure_source": f"pdb:{accession}",
                "reference_fold_class": reference_fold_class,
                "sequence": sequence,
                "sequence_sha256": sha256_sequence(sequence),
                "sequence_length": len(sequence),
                "coordinate_source_url": (
                    f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
                ),
                "coordinate_points": [point.to_dict() for point in points],
                "coordinate_trace_hash": coordinate_trace_hash(points),
                "coordinate_residue_count": len(points),
                "coordinate_coverage": round(len(points) / max(len(sequence), 1), 6),
                "native_contact_map_hash": contact_map_hash(native_pairs),
                "truth_axes": _axis_for_reference_class(reference_fold_class),
            }
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a one-protein RCSB coordinate holdout lock from a local PDB file."
        )
    )
    parser.add_argument("--pdb-file", default=str(DEFAULT_PDB_FILE))
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--pdb-id", default="1UBQ")
    parser.add_argument("--chain-id", default="A")
    parser.add_argument("--source-id", default="ubiquitin_1ubq_a")
    parser.add_argument("--reference-fold-class", default="alpha_beta_mixed")
    args = parser.parse_args()

    print(
        build_rcsb_single_protein_coordinate_holdout(
            pdb_file=Path(args.pdb_file),
            output=Path(args.output),
            pdb_id=args.pdb_id,
            chain_id=args.chain_id,
            source_id=args.source_id,
            reference_fold_class=args.reference_fold_class,
        )
    )


if __name__ == "__main__":
    main()
