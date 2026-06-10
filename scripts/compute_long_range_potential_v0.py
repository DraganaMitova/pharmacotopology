#!/usr/bin/env python3
from __future__ import annotations

"""Compute leave-one-out long-range contact potential from locked benchmark rows.

This is a dependency-light stand-in for a larger PDB reference-set builder. It
uses the locked coordinate benchmark rows and excludes --target-row-id when
provided.
"""

import argparse
import json
from pathlib import Path

from pharmacotopology.folding_long_range_calibrated_verifier import build_leave_one_out_long_range_potential
from pharmacotopology.folding_real_coordinate_visual_benchmark import load_real_coordinate_visual_rows

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BENCHMARK_FILE = REPO_ROOT / "data" / "folding_real_coordinate_visual_8.locked.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-file", default=str(DEFAULT_BENCHMARK_FILE))
    parser.add_argument("--target-row-id", default="")
    parser.add_argument("--source-accession", default="", help="Alternative target selector, e.g. 4AKE:A")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    rows = load_real_coordinate_visual_rows(Path(args.benchmark_file))
    target_row_id = args.target_row_id
    if args.source_accession:
        matches = [r for r in rows if r.source_accession == args.source_accession]
        if not matches:
            raise SystemExit(f"No row matched source accession {args.source_accession!r}")
        target_row_id = matches[0].row_id
    if not target_row_id:
        raise SystemExit("Provide --target-row-id or --source-accession")
    potential = build_leave_one_out_long_range_potential(rows, target_row_id)
    Path(args.out).write_text(json.dumps(potential.to_dict(), indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(potential.to_dict(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
