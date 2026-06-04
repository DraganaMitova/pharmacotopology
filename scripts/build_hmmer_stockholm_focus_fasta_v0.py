from __future__ import annotations

import argparse
import gzip
import sys
from collections import OrderedDict
from pathlib import Path

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


def _focus_alignment_from_rf(*, rf: str, sequence: str, row_id: str) -> str:
    match_columns = [index for index, char in enumerate(rf) if char == "x"]
    if len(match_columns) != len(sequence):
        raise ValueError(
            f"{row_id} has {len(match_columns)} RF match columns, "
            f"but frozen sequence length is {len(sequence)}"
        )
    focus = ["."] * len(rf)
    for residue, index in zip(sequence, match_columns):
        focus[index] = residue
    return "".join(focus)


def build_hmmer_stockholm_focus_fasta_v0(
    *,
    benchmark_file: Path,
    row_id: str,
    stockholm_file: Path,
    output: Path,
    focus_id: str,
    max_records: int,
) -> Path:
    rows = tuple(load_real_coordinate_visual_rows(benchmark_file))
    row_by_id = {row.row_id: row for row in rows}
    row = row_by_id[row_id]
    records, rf = _read_stockholm(stockholm_file)
    if not rf:
        raise ValueError(f"{stockholm_file} does not contain #=GC RF annotation")
    focus_alignment = _focus_alignment_from_rf(
        rf=rf,
        sequence=row.sequence,
        row_id=row_id,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    lines = [f">{focus_id}", focus_alignment]
    for index, (name, sequence) in enumerate(records.items()):
        if index >= max_records:
            break
        if len(sequence) != len(focus_alignment):
            raise ValueError(
                f"{name} alignment width {len(sequence)} does not match "
                f"focus width {len(focus_alignment)}"
            )
        lines.extend((f">{name}", sequence))
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
    args = parser.parse_args()
    output = build_hmmer_stockholm_focus_fasta_v0(
        benchmark_file=Path(args.benchmark_file),
        row_id=args.row_id,
        stockholm_file=Path(args.stockholm_file),
        output=Path(args.output),
        focus_id=args.focus_id,
        max_records=args.max_records,
    )
    print(output)


if __name__ == "__main__":
    main()
