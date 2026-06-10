#!/usr/bin/env python3
from __future__ import annotations

"""Prepare independent target PDB segments for V13 transfer runs.

This script trims downloaded AlphaFold model PDBs to the exact benchmark sequence
segment and renumbers the segment to 1..N. It is a provenance/format repair
step only: it does not read native coordinates and it does not tune selector
rules.
"""

import argparse
import hashlib
import json
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    RealCoordinateVisualRow,
    load_real_coordinate_visual_rows,
)

AA3_TO_1 = {
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
    "SEC": "U",
    "PYL": "O",
}


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _load_row(benchmark_file: Path, source_accession: str) -> RealCoordinateVisualRow:
    for row in load_real_coordinate_visual_rows(benchmark_file):
        if row.source_accession == source_accession:
            return row
    raise SystemExit(f"no row found for source_accession={source_accession!r} in {benchmark_file}")


def _iter_first_model_atom_lines(pdb_path: Path) -> Iterable[str]:
    """Yield ATOM/HETATM lines from the first model, or all lines if no MODEL tags.

    Some RCSB structures are NMR multi-model PDBs. The target-preparation step
    only needs one coordinate model; otherwise duplicate CA atoms make the
    prepared segment look corrupt.
    """
    in_first_model = False
    seen_model_tag = False
    first_model_done = False
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("MODEL"):
            seen_model_tag = True
            if first_model_done:
                break
            in_first_model = True
            continue
        if line.startswith("ENDMDL") and in_first_model:
            first_model_done = True
            in_first_model = False
            continue
        if seen_model_tag and not in_first_model:
            continue
        if line.startswith("ATOM") or line.startswith("HETATM"):
            yield line


def _parse_residue_order(pdb_path: Path, chain_id: str = "A") -> tuple[list[tuple[str, int, str, str]], str]:
    """Return ordered residues and one-letter sequence for a chain.

    Residue key is (chain, residue_number, insertion_code, residue_name).
    The sequence is based only on residue order, not on coordinates from the
    benchmark/native structure.
    """
    if not pdb_path.is_file():
        raise SystemExit(f"raw_target_pdb_missing: {pdb_path}")
    residues: "OrderedDict[tuple[str, int, str, str], str]" = OrderedDict()
    for line in _iter_first_model_atom_lines(pdb_path):
        if not line.startswith("ATOM"):
            continue
        chain = line[21].strip() or "A"
        if chain_id and chain != chain_id:
            continue
        resname = line[17:20].strip().upper()
        icode = line[26].strip()
        try:
            resseq = int(line[22:26].strip())
        except Exception:
            continue
        aa = AA3_TO_1.get(resname, "X")
        key = (chain, resseq, icode, resname)
        residues.setdefault(key, aa)
    if not residues:
        raise SystemExit(f"no_residues_found_in_raw_target_pdb_chain_{chain_id}: {pdb_path}")
    ordered_keys = list(residues.keys())
    sequence = "".join(residues[key] for key in ordered_keys)
    return ordered_keys, sequence


def _find_exact_segment(raw_sequence: str, target_sequence: str) -> int:
    start = raw_sequence.find(target_sequence)
    if start >= 0:
        return start
    # A common mature-protein case: the AF model has an initiator Met before the
    # benchmark segment. Keep this as a format repair, not a scoring rule.
    if raw_sequence.startswith("M" + target_sequence):
        return 1
    raise SystemExit(
        "target_sequence_not_found_in_raw_af_model: "
        f"target_len={len(target_sequence)} raw_len={len(raw_sequence)} "
        f"raw_prefix={raw_sequence[:30]!r} target_prefix={target_sequence[:30]!r}"
    )


def _renumber_atom_line(line: str, *, new_resseq: int, new_chain: str = "A") -> str:
    # PDB columns: chain id at 22 (0-index 21), resseq columns 23-26 (22:26).
    if len(line) < 54:
        line = line.rstrip("\n").ljust(54)
    return f"{line[:21]}{new_chain}{new_resseq:4d}{line[26:]}"


def _write_segment_pdb(
    *,
    raw_pdb: Path,
    out_pdb: Path,
    segment_keys: list[tuple[str, int, str, str]],
) -> int:
    key_to_new_index = {key: idx for idx, key in enumerate(segment_keys, start=1)}
    out_lines: list[str] = []
    ca_count = 0
    for line in _iter_first_model_atom_lines(raw_pdb):
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        chain = line[21].strip() or "A"
        resname = line[17:20].strip().upper()
        icode = line[26].strip()
        try:
            resseq = int(line[22:26].strip())
        except Exception:
            continue
        key = (chain, resseq, icode, resname)
        new_index = key_to_new_index.get(key)
        if new_index is None:
            continue
        out_lines.append(_renumber_atom_line(line, new_resseq=new_index, new_chain="A"))
        if line[12:16].strip() == "CA":
            ca_count += 1
    if ca_count != len(segment_keys):
        raise SystemExit(
            f"segment_ca_count_mismatch: out={out_pdb} ca_count={ca_count} expected={len(segment_keys)}"
        )
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdb.write_text("\n".join(out_lines + ["TER", "END", ""]), encoding="utf-8")
    return ca_count


def prepare_one(
    *,
    raw_pdb: Path,
    out_pdb: Path,
    benchmark_file: Path,
    source_accession: str,
    chain_id: str,
) -> dict[str, object]:
    row = _load_row(benchmark_file, source_accession)
    residue_keys, raw_sequence = _parse_residue_order(raw_pdb, chain_id=chain_id)
    start = _find_exact_segment(raw_sequence, row.sequence)
    end = start + len(row.sequence)
    segment_keys = residue_keys[start:end]
    if len(segment_keys) != row.sequence_length:
        raise SystemExit(
            f"segment_length_mismatch: got={len(segment_keys)} expected={row.sequence_length} source={source_accession}"
        )
    ca_count = _write_segment_pdb(raw_pdb=raw_pdb, out_pdb=out_pdb, segment_keys=segment_keys)
    provenance = {
        "kind": "v13_independent_target_segment_provenance_v0",
        "source_accession": source_accession,
        "row_id": row.row_id,
        "benchmark_file": str(benchmark_file),
        "raw_target_pdb": str(raw_pdb),
        "prepared_target_pdb": str(out_pdb),
        "raw_chain_id": chain_id,
        "raw_sequence_length": len(raw_sequence),
        "segment_start_residue_order_index_1based": start + 1,
        "segment_end_residue_order_index_1based": end,
        "prepared_ca_count": ca_count,
        "expected_sequence_length": row.sequence_length,
        "sequence_exact_match": True,
        "sequence_sha256": _sha256_text(row.sequence),
        "coordinate_truth_used_to_prepare_target": False,
        "native_truth_used_to_prepare_target": False,
        "selector_rules_modified": False,
        "claim_allowed": False,
    }
    out_pdb.with_suffix(out_pdb.suffix + ".provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare V13 independent target PDB segment files.")
    parser.add_argument("--repo-root", default=str(REPO_ROOT))
    parser.add_argument("--chain-id", default="A")
    parser.add_argument(
        "--ubq-raw-pdb",
        default="data/independent_contact_sources/raw/AF-P62988-F1-model_v4.pdb",
    )
    parser.add_argument(
        "--ubq-out-pdb",
        default="data/independent_contact_sources/AF-P62988-F1-model_v4_1UBQ_mature_segment.pdb",
    )
    parser.add_argument(
        "--cll-raw-pdb",
        default="data/independent_contact_sources/raw/AF-P0DP23-F1-model_v4.pdb",
    )
    parser.add_argument(
        "--cll-out-pdb",
        default="data/independent_contact_sources/AF-P0DP23-F1-model_v4_1CLL_calmodulin_segment.pdb",
    )
    parser.add_argument(
        "--ubq-benchmark-file",
        default="data/folding_real_coordinate_holdout_1ubq.locked.json",
    )
    parser.add_argument(
        "--cll-benchmark-file",
        default="data/folding_real_coordinate_visual_8.locked.json",
    )
    parser.add_argument("--only", choices=("all", "1UBQ", "1CLL"), default="all")
    args = parser.parse_args()

    repo_root = Path(args.repo_root)
    results: list[dict[str, object]] = []
    if args.only in {"all", "1UBQ"}:
        results.append(
            prepare_one(
                raw_pdb=repo_root / args.ubq_raw_pdb,
                out_pdb=repo_root / args.ubq_out_pdb,
                benchmark_file=repo_root / args.ubq_benchmark_file,
                source_accession="1UBQ:A",
                chain_id=args.chain_id,
            )
        )
    if args.only in {"all", "1CLL"}:
        results.append(
            prepare_one(
                raw_pdb=repo_root / args.cll_raw_pdb,
                out_pdb=repo_root / args.cll_out_pdb,
                benchmark_file=repo_root / args.cll_benchmark_file,
                source_accession="1CLL:A",
                chain_id=args.chain_id,
            )
        )
    print(json.dumps({"prepared": results}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
