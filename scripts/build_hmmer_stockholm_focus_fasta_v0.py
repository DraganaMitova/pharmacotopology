from __future__ import annotations

import argparse
import gzip
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    load_real_coordinate_visual_rows,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")


def _open_text(path: Path):
    if path.read_bytes()[:2] == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _read_stockholm(path: Path) -> tuple[OrderedDict[str, str], str]:
    records: OrderedDict[str, list[str]] = OrderedDict()
    reference_fragments: list[str] = []
    with _open_text(path) as file:
        for line in file:
            stripped = line.strip()
            if not stripped or stripped == "//":
                continue
            parts = stripped.split()
            if len(parts) < 2:
                continue
            if stripped.startswith("#=GC RF") and len(parts) >= 3:
                reference_fragments.append("".join(parts[2:]))
                continue
            if stripped.startswith("#"):
                continue
            name, fragment = parts[0], parts[1]
            records.setdefault(name, []).append(fragment)
    return OrderedDict((name, "".join(parts)) for name, parts in records.items()), "".join(
        reference_fragments
    )


def _focus_alignment_from_rf(
    *,
    rf: str,
    sequence: str,
    row_id: str,
    domain_start: int,
    domain_end: int,
) -> str:
    if domain_start < 1 or domain_end < domain_start or domain_end > len(sequence):
        raise ValueError(
            f"{row_id} invalid domain range {domain_start}-{domain_end} "
            f"for sequence length {len(sequence)}"
        )
    prefix = sequence[: domain_start - 1].lower()
    domain_sequence = sequence[domain_start - 1 : domain_end]
    suffix = sequence[domain_end:].lower()
    match_columns = [index for index, char in enumerate(rf) if char == "x"]
    if len(match_columns) != len(domain_sequence):
        raise ValueError(
            f"{row_id} has {len(match_columns)} RF match columns, "
            f"but domain length is {len(domain_sequence)}"
        )
    focus = ["."] * len(rf)
    for residue, index in zip(domain_sequence, match_columns):
        focus[index] = residue
    return prefix + "".join(focus) + suffix


def _focus_alignment_from_exact_record(
    *,
    records: OrderedDict[str, str],
    sequence: str,
    row_id: str,
    domain_start: int,
    domain_end: int,
    focus_record_name: Optional[str],
) -> str:
    domain_sequence = sequence[domain_start - 1 : domain_end]
    selected_name = ""
    selected_alignment = ""
    for name, alignment in records.items():
        if focus_record_name and name != focus_record_name:
            continue
        letters = "".join(char for char in alignment if char.isalpha()).upper()
        if letters == domain_sequence:
            selected_name = name
            selected_alignment = alignment
            break
    if not selected_alignment:
        label = focus_record_name or "any exact alignment record"
        raise ValueError(
            f"{row_id} could not find {label} matching domain "
            f"{domain_start}-{domain_end}"
        )
    prefix = sequence[: domain_start - 1].lower()
    suffix = sequence[domain_end:].lower()
    domain_focus = "".join(
        char.upper() if char.isalpha() else char for char in selected_alignment
    )
    return prefix + domain_focus + suffix


def build_hmmer_stockholm_focus_fasta_v0(
    *,
    benchmark_file: Path,
    row_id: str,
    stockholm_file: Path,
    output: Path,
    focus_id: str,
    max_records: int,
    domain_start: int = 1,
    domain_end: Optional[int] = None,
    focus_record_name: Optional[str] = None,
) -> Path:
    rows = tuple(load_real_coordinate_visual_rows(benchmark_file))
    row_by_id = {row.row_id: row for row in rows}
    row = row_by_id[row_id]
    domain_end = domain_end or row.sequence_length
    records, rf = _read_stockholm(stockholm_file)
    if rf:
        focus_alignment = _focus_alignment_from_rf(
            rf=rf,
            sequence=row.sequence,
            row_id=row_id,
            domain_start=domain_start,
            domain_end=domain_end,
        )
    else:
        focus_alignment = _focus_alignment_from_exact_record(
            records=records,
            sequence=row.sequence,
            row_id=row_id,
            domain_start=domain_start,
            domain_end=domain_end,
            focus_record_name=focus_record_name,
        )
    prefix_gap = "-" * (domain_start - 1)
    suffix_gap = "-" * (row.sequence_length - domain_end)
    output.parent.mkdir(parents=True, exist_ok=True)
    lines = [f">{focus_id}", focus_alignment]
    for index, (name, sequence) in enumerate(records.items()):
        if index >= max_records:
            break
        padded_sequence = prefix_gap + sequence + suffix_gap
        if len(padded_sequence) != len(focus_alignment):
            raise ValueError(
                f"{name} padded alignment width {len(padded_sequence)} does not match "
                f"focus width {len(focus_alignment)}"
            )
        lines.extend((f">{name}", padded_sequence))
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build a focus-first FASTA from an EBI HMMER Stockholm export by "
            "placing the frozen benchmark sequence onto #=GC RF match columns."
        )
    )
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--row-id", required=True)
    parser.add_argument("--stockholm-file", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--focus-id", required=True)
    parser.add_argument("--max-records", type=int, default=2000)
    parser.add_argument("--domain-start", type=int, default=1)
    parser.add_argument("--domain-end", type=int)
    parser.add_argument("--focus-record-name")
    args = parser.parse_args()
    output = build_hmmer_stockholm_focus_fasta_v0(
        benchmark_file=Path(args.benchmark_file),
        row_id=args.row_id,
        stockholm_file=Path(args.stockholm_file),
        output=Path(args.output),
        focus_id=args.focus_id,
        max_records=args.max_records,
        domain_start=args.domain_start,
        domain_end=args.domain_end,
        focus_record_name=args.focus_record_name,
    )
    print(output)


if __name__ == "__main__":
    main()
