from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from pharmacotopology.folding_reference_loader import (  # noqa: E402
    load_folding_reference_dataset,
)
from pharmacotopology.folding_structure_benchmark import (  # noqa: E402
    build_structure_evidence_rows,
    write_structure_evidence_file,
)


DEFAULT_BENCHMARK_FILE = Path("data/folding_benchmarks_real_10.locked.json")
DEFAULT_OUTPUT_FILE = Path("data/folding_benchmarks_real_10_structure_evidence.json")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Extract structure-derived folding topology signatures from local "
            "PDB coordinate files plus DisProt disorder-reference rows."
        )
    )
    parser.add_argument(
        "--benchmark-file",
        default=str(DEFAULT_BENCHMARK_FILE),
        help="Locked external benchmark file.",
    )
    parser.add_argument(
        "--pdb-dir",
        required=True,
        help="Directory containing <PDB_ID>.pdb coordinate files.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT_FILE),
        help="Output structure evidence JSON path.",
    )
    args = parser.parse_args()

    benchmark_file = Path(args.benchmark_file)
    dataset = load_folding_reference_dataset(benchmark_file, require_external=True)
    rows = build_structure_evidence_rows(dataset.references, pdb_dir=Path(args.pdb_dir))
    output = write_structure_evidence_file(
        rows,
        Path(args.output),
        source_benchmark_file=benchmark_file,
    )
    print(output)


if __name__ == "__main__":
    main()
