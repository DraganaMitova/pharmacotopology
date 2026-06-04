from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_real_coordinate_visual_benchmark import (  # noqa: E402
    write_real_coordinate_visual_lock,
)


DEFAULT_SOURCE_BENCHMARK_FILE = Path("data/folding_benchmarks_real_10.locked.json")
DEFAULT_OUTPUT_FILE = Path("data/folding_real_coordinate_visual_8.locked.json")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build the locked real-coordinate visual contact benchmark from "
            "local RCSB PDB coordinate files. The generated lock stores "
            "minimal C-alpha coordinate traces, not full PDB files."
        )
    )
    parser.add_argument(
        "--source-benchmark-file",
        default=str(DEFAULT_SOURCE_BENCHMARK_FILE),
    )
    parser.add_argument(
        "--pdb-dir",
        required=True,
        help="Directory containing <PDB_ID>.pdb coordinate files.",
    )
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT_FILE))
    args = parser.parse_args()

    output = write_real_coordinate_visual_lock(
        source_benchmark_file=Path(args.source_benchmark_file),
        pdb_dir=Path(args.pdb_dir),
        output_path=Path(args.output),
    )
    print(output)


if __name__ == "__main__":
    main()
